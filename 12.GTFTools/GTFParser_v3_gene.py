import pandas as pd
import os
import gzip
import requests

class ParseGTFAttributes:
    def __init__(self, url, folder_path, chr_key_path):
        # Initialize the class with file paths for GTF and chromosome key files
        self.url = url
        self.file_path = folder_path + '/' + url.split('/')[-1]
        self.db = None
        self.chr_key_path = chr_key_path

    def download_gtf_file(self):
        response = requests.get(self.url, stream=True)
        with open(self.file_path, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
        return self.file_path

    def ReadParseGTF(self):
        # This function parses any GTF file, extracts attributes, and creates a new DataFrame
        with gzip.open(self.file_path, 'rb') as file:
            df = pd.read_csv(file, delimiter='\t', header=None, comment='#', dtype='str')

        # Naming the columns of the DataFrame
        df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

        def parse_attribute(attribute):
            # This function splits the attribute column and creates a dictionary of attributes
            attributes = attribute.split(';')
            parsed_attributes = {}
            for attr in attributes:
                attr = attr.strip()
                if attr:
                    parts = attr.split(' ')
                    if len(parts) >= 2:
                        key, value = parts[0], ' '.join(parts[1:]).strip('"')
                        if key in parsed_attributes:
                            parsed_attributes[key].append(value)
                        else:
                            parsed_attributes[key] = [value]
                    else:
                        parsed_attributes[parts[0]] = None
            for key, value in parsed_attributes.items():
                if isinstance(value, list):
                    parsed_attributes[key] = '|'.join(value)
            return parsed_attributes

        # Applying the attribute parsing function to the 'attribute' column
        attribute_dict = df['attribute'].apply(parse_attribute)
        parsed_attributes_df = pd.DataFrame(attribute_dict.tolist())

        # Combining the parsed attributes DataFrame with the original DataFrame
        parsed_df = pd.concat([df, parsed_attributes_df], axis=1)
        parsed_df.drop(['attribute'], axis=1, inplace=True)
        parsed_df.fillna('-', inplace=True)
        self.db = parsed_df

        # Exporting the parsed GTF to a new file
        output_file_path = self.file_path.rsplit('.', 1)[0] + '_parsed.txt'
        parsed_df.to_csv(output_file_path, sep="\t", index=False)
        return output_file_path
    
    print("finished parsing the gtf file")
    
    def create_bed_file(self, data_from_ncbi = False):
        # This function creates a BED file containing gene names and their locations
        if data_from_ncbi:

            chr_key = pd.read_csv(self.chr_key_path, sep='\t')

            def extract_chromosome_name(row):
                # This function extracts chromosome names from chromosome aliases using the chromosome key
                chromosome_id = row['seqname']

                match = chr_key[chr_key['alias'] == chromosome_id]
                if not match.empty:
                    return match['chr_name'].iloc[0]
                return chromosome_id

            # Applying the extract_chromosome_name function to create a new column 'chr_name'
            self.db['chr_name'] = self.db.apply(extract_chromosome_name, axis=1)

            # Selecting only the 'gene' features from the DataFrame
            bed_df = self.db.loc[self.db['feature'] == 'gene', ['chr_name', 'start', 'end', 'gene', 'gene_id']]
            bed_df.reset_index(drop=True,inplace=True)

            # Grouping genes to remove versions denoted as 'gene_name_{number}'
            gene_groups = bed_df.groupby('gene')
            filtered_groups = [group for _, group in gene_groups if len(group) > 1]

            final_rows = []
            for group in filtered_groups:
                without_underscore_number = group[~group['gene_id'].str.contains(r'_\d+$')]
                if len(without_underscore_number) > 0:
                    final_rows.append(without_underscore_number)

            final_rows.extend([group for _, group in gene_groups if len(group) == 1])
            result_df = pd.concat(final_rows)
            result_df.drop('gene_id', axis=1, inplace=True)
            # Writing the BED file and sorting it
            bed_file_path = self.file_path.rsplit('.', 1)[0] + '.gene_names.bed'
            result_df.to_csv(bed_file_path, sep="\t", index=False, header=False)

            sorted_bed_file_path = bed_file_path.rsplit('.', 1)[0] + '.sorted.bed'
            gzipped_bed_file_path = sorted_bed_file_path + '.gz'

            os.system(f"sort -k1,1 -k2,2n {bed_file_path} > {sorted_bed_file_path}")
            os.system(f"bgzip -c {sorted_bed_file_path} > {gzipped_bed_file_path}")
            os.system(f"tabix -p bed {gzipped_bed_file_path}")

            return gzipped_bed_file_path
        
        # Selecting only the 'gene' features from the DataFrame

        chr_key = pd.read_csv(self.chr_key_path, sep='\t') 
        
        def extract_chromosome_name(row):
            # This function extracts chromosome names from chromosome aliases using the chromosome key
            chromosome_id = row['seqname']

            match = chr_key[chr_key['alias'] == chromosome_id]
            if not match.empty:
                return match['chr_name'].iloc[0]
            return chromosome_id
        
        self.db['chr_name'] = self.db.apply(extract_chromosome_name, axis=1)

        bed_df = self.db.loc[self.db['feature'] == 'gene', ['chr_name', 'start', 'end', 'gene_name', 'gene_id']] #first bed
        bed_df = bed_df.astype(str)
        bed_df.reset_index(drop=True, inplace=True)
        bed_df.replace('-', pd.NA, inplace=True)
        
        # If there is not a gene name for an entry, add gene id.
        bed_df['gene_name'].fillna(bed_df['gene_id'], inplace=True)
        bed_df.drop('gene_id', axis=1, inplace=True)
        bed_df.sort_values(by=['chr_name', 'start'], inplace=True)
        
        # Writing the BED file and sorting it
        bed_file_path = self.file_path.rsplit('.', 1)[0] + '.gene_names.bed'
        bed_df.to_csv(bed_file_path, sep="\t", index=False, header=False)
        
        sorted_bed_file_path = bed_file_path.rsplit('.', 1)[0] + '.sorted.bed'
        gzipped_bed_file_path = sorted_bed_file_path + '.gz'
        
        os.system(f"sort -k1,1 -k2,2n {bed_file_path} > {sorted_bed_file_path}") 
        os.system(f"bgzip -c {sorted_bed_file_path} > {gzipped_bed_file_path}")
        os.system(f"tabix -p bed {gzipped_bed_file_path}")


        return gzipped_bed_file_path


# Usage - GTF files from NCBI

#url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz'
#folder_path = "/home/genwork2/Mert/02.BEDS/NCBI-hg19-genees"
#chr_key_path = "/home/genwork2/Mert/Chralias/final_ultimate_merged_chrkeys_hg19.txt"
#parser = ParseGTFAttributes(url ,folder_path, chr_key_path)
#file_path = parser.download_gtf_file()
#output_file_path = parser.ReadParseGTF()
#bed_file_path = parser.create_bed_file(data_from_ncbi=True)

#print('GTF file downloaded at:', file_path)
#print("Parsed file saved at:", output_file_path)
#print("BED file saved at:", bed_file_path)

#Usage - GTF files from ENSEMBL

url = "https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
file_path = "/home/genwork2/Mert/02.BEDS/ENSEMBL-hg19-genes"
chr_key_path = '/home/genwork2/Mert/Chralias/final_ultimate_merged_chrkeys_hg19.txt'
parser = ParseGTFAttributes(url, file_path, chr_key_path)
file_path = parser.download_gtf_file()
output_file_path = parser.ReadParseGTF()
bed_file_path = parser.create_bed_file(data_from_ncbi=False)

print('GTF file downloaded at:', file_path)
print("Parsed file saved at:", output_file_path)
print("BED file saved at:", bed_file_path)

