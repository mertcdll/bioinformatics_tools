import pandas as pd
import os
import json
import urllib3
import ssl
import requests

# Custom http Adapter to handle ssl connections
class CustomHttpAdapter(requests.adapters.HTTPAdapter):

    def __init__(self, ssl_context=None, **kwargs):
        self.ssl_context = ssl_context
        super().__init__(**kwargs)

    def init_poolmanager(self, connections, maxsize, block=False):
        self.poolmanager = urllib3.poolmanager.PoolManager(
            num_pools=connections, maxsize=maxsize,
            block=block, ssl_context=self.ssl_context)

# Class for downloading hpo disease data
class HPODiseaseDownloader:
    def __init__(self, phenotype_file, output_folder):
        # Read phenotype file and initialize attributes
        self.df = pd.read_csv(phenotype_file, sep='\t', comment='#')
        self.ptlist = self.df['database_id'].value_counts().index.to_list()
        self.output_folder = output_folder
    # Method to get a legacy session with ssl context
    def get_legacy_session(self):
        ctx = ssl.create_default_context(ssl.Purpose.SERVER_AUTH)
        ctx.options |= 0x4
        session = requests.session()
        session.mount('https://', CustomHttpAdapter(ctx))
        return session
    # Method to download disease data for a specific disease ID
    def download_disease_data(self, disease_id):
        session = self.get_legacy_session()
        url = f"https://hpo.jax.org/api/hpo/disease/{disease_id}"
        response = session.get(url)
        # Create output folder if it doesn't exist
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        # Save the downloaded data as json file
        filename = os.path.join(self.output_folder, f"{disease_id}.json")

        with open(filename, 'w') as file:
            json.dump(response.json(), file)
    # Method to download data for all diseases
    def download_all_diseases(self):
        for disease_id in self.ptlist:
            try:
                self.download_disease_data(disease_id)
                print(disease_id)
            except Exception as e:
                print(f"Error in disease ID {disease_id}: {e}")
                continue

# Class for generating bed files from downloaded HPO disease data
class HPOBedFileGenerator:
    def __init__(self, folder_path, gene_bed_file):
        self.folder_path = folder_path
        self.gene_bed_file = gene_bed_file
        self.data_list = []
        self.df = None
        self.merged_df = None
        self.merged_df_wg = None
        self.merged_df_wog = None
        self.wg_bed = None
        self.wog_bed = None
        self.inh_key = {
            "-": "-",
            "Autosomal recessive inheritance": "AR",
            "Autosomal dominant inheritance": "AD",
            "X-linked recessive inheritance": "XLR",
            "Sporadic": "SP",
            "Typified by incomplete penetrance": "IP",
            "X-linked inheritance": "XL",
            "Typified by somatic mosaicism": "SOM",
            "X-linked dominant inheritance": "XLD",
            "Non-Mendelian inheritance": "NMI",
            "Mitochondrial inheritance": "M",
            "Polygenic inheritance": "PI",
            "Contiguous gene syndrome": "CGS",
            "Typically de novo": "DN",
            "Digenic inheritance": "DI",
            "Genetic anticipation": "GA",
            "Y-linked inheritance": "YL",
            "Typified by age-related disease onset": "ARDO",
            "Sex-limited expression": "SLX",
            "Male-limited expression": "MLE",
            "Oligogenic inheritance": "OI",
            "Genetic anticipation with paternal anticipation bias": "GAPA",
            "Autosomal dominant inheritance with maternal imprinting": "ADMI",
            "Female-limited expression": "FLE",
            "Uniparental isodisomy": "UPiD",
            "Typified by high penetrance": "HP",
            "Uniparental disomy": "UPD"
        }


    # Load data from downloaded json files
    def load_data(self):
        for filename in os.listdir(self.folder_path):
            if filename.endswith('.json'):
                file_path = os.path.join(self.folder_path, filename)

                with open(file_path, 'r') as file:
                    dis_data = json.load(file)

                    try:
                        gene_assoc_list = dis_data['geneAssoc']
                    except (IndexError, KeyError):
                        gene_assoc_list = []

                    if not gene_assoc_list:
                        gene_assoc_list = [{}]

                    for gene_assoc in gene_assoc_list:
                        try:
                            gene_symbol = gene_assoc.get('geneSymbol', '-')
                            gene_id = gene_assoc.get('geneId', '-')
                        except (IndexError, KeyError):
                            gene_symbol = '-'
                            gene_id = '-'

                        try:
                            disease_name = dis_data['disease']['diseaseName']
                        except KeyError:
                            disease_name = '-'

                        try:
                            db_id = dis_data['disease']['dbId']
                        except KeyError:
                            db_id = '-'

                        try:
                            disease_id = dis_data['disease']['diseaseId']
                        except KeyError:
                            disease_id = '-'

                        try:
                            inheritance = '-'
                            for i in range(len(dis_data['catTermsMap'])):
                                if dis_data['catTermsMap'][i]['catLabel'] == 'Inheritance':
                                    names = []
                                    for j in range(len(dis_data['catTermsMap'][i]['terms'])):
                                        names.append(dis_data['catTermsMap'][i]['terms'][j]['name'])

                                    inheritance = "|".join(map(str, names))
                                break

                        except (IndexError, KeyError):
                            inheritance = '-'

                        self.data_list.append({
                            'gene symbol': gene_symbol,
                            'geneid': gene_id,
                            'disease name': disease_name,
                            'dbid': db_id,
                            'diseaseid': disease_id,
                            'inheritance': inheritance
                        })

        self.df = pd.DataFrame(self.data_list)
    # Merge data from json files with gene bed file - gene bed file contains gene coordinates.
    def merge_data(self):
        gene_bed = pd.read_csv(self.gene_bed_file, sep='\t', header=None)
        gene_bed.columns = ['chr', 'start_pos', 'end_pos', 'gene_name']

        self.merged_df = self.df.merge(gene_bed, left_on='gene symbol', right_on='gene_name', how='left')

        self.merged_df.drop('gene_name', axis=1, inplace=True)
        self.merged_df.fillna('-', inplace=True)
    # Create info column based on merged data
    def create_info_column(self):
        self.merged_df['db'] = self.merged_df['diseaseid'].str.split(':', expand=True).iloc[:, 0]

        self.merged_df = self.merged_df.astype(str)

        def to_abbrv(phrase):
            components = phrase.split('|')
            abbreviations = [self.inh_key[component.strip()] for component in components]
            return '|'.join(abbreviations)

        self.merged_df['inheritance_abb'] = self.merged_df['inheritance'].apply(to_abbrv)


        def format_element(element):
            if pd.notna(element):
                element = element.replace(' ', '_')
                if ',' in element:
                    return element.replace(',', '')
            return element


        self.merged_df['info'] = (
            'SYMBOL=' + self.merged_df['gene symbol'] + ';' +
            'GENEID=' + self.merged_df['geneid'] + ';' +
            'DISEASENAME=' + self.merged_df['disease name'].apply(format_element) + ';' +
            'DB=' + self.merged_df['db'] + ';' +
            'DBID=' + self.merged_df['dbid'] + ';' +
            'DISEASEID=' + self.merged_df['diseaseid'] + ';' +
            'INHERITANCE=' + self.merged_df['inheritance_abb']
        )
    # Filter merged data for further processing
    def filter_data(self):
        self.merged_df_wg = self.merged_df[~((self.merged_df['gene symbol'] == '-') & (self.merged_df['geneid'] == '-'))]
        self.merged_df_wog = self.merged_df[((self.merged_df['gene symbol'] == '-') & (self.merged_df['geneid'] == '-'))]
    # Convert data to BED format
    def convert_to_bed_format(self):

       
        self.merged_df_wg['start_pos'] = self.merged_df_wg['start_pos'].str.split('.', expand=True).iloc[:, 0]
        self.merged_df_wg['end_pos'] = self.merged_df_wg['end_pos'].str.split('.', expand=True).iloc[:, 0]
        self.merged_df_wog['start_pos'] = self.merged_df_wog['start_pos'].str.split('.', expand=True).iloc[:, 0]
        self.merged_df_wog['end_pos'] = self.merged_df_wog['end_pos'].str.split('.', expand=True).iloc[:, 0]
    

        self.wg_bed = self.merged_df_wg[['chr', 'start_pos', 'end_pos', 'info']]
        self.wog_bed = self.merged_df_wog[['chr', 'start_pos', 'end_pos', 'info']]
    # Write BED files for the filtered data
    def write_bed_files(self):
        self.wg_bed.to_csv('wg_bed_diseases.bed', index=False, header=None, sep='\t')
        os.system(f"sort -k1,1 -k2,2n {'wg_bed_diseases.bed'} > {'wg_bed_diseases_sorted.bed'}")
        os.system(f"bgzip -c {'wg_bed_diseases_sorted.bed'} > {'wg_bed_diseases_sorted.bed.gz'}")
        os.system(f"tabix -p bed {'wg_bed_diseases_sorted.bed.gz'}")

        self.wog_bed.to_csv('wog_bed_diseases.bed', index=False, header=None, sep='\t')
        os.system(f"sort -k1,1 -k2,2n {'wog_bed_diseases.bed'} > {'wog_bed_diseases_sorted.bed'}")
        os.system(f"bgzip -c {'wog_bed_diseases_sorted.bed'} > {'wog_bed_diseases_sorted.bed.gz'}")
        os.system(f"tabix -p bed {'wog_bed_diseases_sorted.bed.gz'}")
    # Process data through all steps        
    def process_data(self):
        self.load_data()
        self.merge_data()
        self.create_info_column()
        self.filter_data()
        self.convert_to_bed_format()
        self.write_bed_files()

# Path to the phenotype file containing database IDs
phenotype_file = '/home/dell/Documents/7.DiseaseAssociatedGenes/phenotype.hpoa'
# Output folder for downloaded hpo disease data (You'll have to set it)
output_folder = 'Downloads'
downloader = HPODiseaseDownloader(phenotype_file, output_folder)
# Download data for all diseases
downloader.download_all_diseases()

# Path to the gene bed file containing gene coordinates
gene_bed_file = '/home/dell/Documents/NecessaryDocuments/Genelocs/ncbi_ensembl_merged_final.txt'
generator = HPOBedFileGenerator(output_folder, gene_bed_file)
# Process and generate bed files
generator.process_data()