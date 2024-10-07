### This one is faster from v2 but v2 is more robust...

import pandas as pd
import gzip
import re
import os

class VCFGenerator:
    def __init__(self, variation_file, chrmap_file):
        self.variation_file = variation_file
        self.chrmap_file = chrmap_file
        self.chrmap_df = None
        self.df = None
        self.vcf_df = None

    ## Reading the variation file downloaded from uniprot
    def read_variation_file(self):
        
        def skip_until_hashtag(file):
            for line in file:
                if line.startswith('#'):
                    break

        
        with gzip.open(self.variation_file, 'rt') as file:
            skip_until_hashtag(file)
            self.df = pd.read_csv(file, delimiter='\t')

        str_cols = []
        for item in self.df.columns:
            strd = item.strip()
            str_cols.append(strd)

        self.df.columns = str_cols
        self.df = self.df.iloc[1:(len(self.df) - 4),]
        self.df = self.df[~((self.df['Clinical Significance'] == '-') & 
                  (self.df['Phenotype/Disease'] == '-') & 
                  (self.df['Phenotype/Disease Source'] == '-'))]
        self.df.reset_index(drop=True, inplace=True)

    ## Extracting chromosome names
    # Import chromosome key

    def read_chrmap_file(self):
        self.chrmap_df = pd.read_csv(self.chrmap_file, sep='\s+', header=None)
        self.chrmap_df.columns = ['ID', 'NAME']
        # Adds CHR_ prefix because this dataset has CHR_ prefix before its aliases.
        def add_prefix(row):
            if row.startswith("HG") or row.startswith("HS"):
                return "CHR_" + row
            return row
        
        self.chrmap_df['ID'] = self.chrmap_df['ID'].apply(add_prefix)
    def process_chromosome_names(self):
        # Extract chromosome name from cytogenetic band and chromosome coordinate column.
        # Function will extract chr name from cytogenetic band if it's present.
        # If cytogenetic band is not present, it extracts chr name from chromosome coordinate column using chr key. 
        def extract_chromosome_name(row):
        

            cytogenetic_band = row['Cytogenetic Band']
            chromosome_coordinate = row['Chromosome Coordinate']

            if pd.notna(cytogenetic_band) and '-' not in cytogenetic_band:
                return 'chr' + re.split(r'(q|p)', cytogenetic_band)[0]

            chromosome_id = re.split(r':g.|:m.', chromosome_coordinate)[0]
            match = self.chrmap_df[self.chrmap_df['ID'] == chromosome_id]

            if not match.empty:
                return match['NAME'].iloc[0]

            return chromosome_id
        # Apply the function and create a new column having chr names
        self.df['Chromosome'] = self.df.apply(extract_chromosome_name, axis=1)

    ## Getting positions

    # Devise a regex pattern to capture position from chromosome coordinates.
    def get_genomic_positions(self):
        pattern = r'[gm]\.(\d+)(?:_(\d+))?'
        genomic_positions = []
        for item in self.df['Chromosome Coordinate']:
            match = re.search(pattern, item)
            if match:
                groups = match.groups()
                genomic_position = '_'.join(gp for gp in groups if gp is not None)
                genomic_positions.append(genomic_position)
            else:
                genomic_positions.append(None)
        self.df['genomic_position']= genomic_positions
        # Get only start position
        self.df['genomic_position'] = self.df['genomic_position'].apply(lambda x: x.split('_')[0] if '_' in str(x) else x)

    ##GET RF AND ALT SEQS

    def parse_genomic_alterations(self):
        #Use regex pattern to extract genomic alterations from chromosome coordinate column.
        self.df['change'] = self.df['Chromosome Coordinate'].str.extract(r'(\[\d+\]|[ACGTRDYH]>[ACGTRDYH]|[^0-9]+)$') #change is ok now.

        # Functions to parse genomic alterations

        def parse_snv(change):
            match = re.match(r'([A-Z])>([A-Z])$', change)
            if match:
                return match.group(1), match.group(2)
            return '-', '-'

        def parse_del(change):
            match = re.match(r'del([A-Z]+)$', change)
            if match:
                return match.group(1), '-'
            return '-', '-'

        def parse_ins(change):
            match = re.match(r'ins([A-Z]+)$', change)
            if match:
                return '-', match.group(1)
            return '-', '-'

        def parse_delins(change):
            match = re.match(r'del([A-Z]+)ins([A-Z]+)$', change)
            if match:
                return match.group(1), match.group(2)
            return '-', '-'

        def parse_other(change):  # this one assigns a hyphen if genomic alteration is delinsNNN, dup, del, inv, ins, [#] or =
            return '-', '-'
        
        self.df['change'] = self.df['change'].astype(str)

        self.df['ref'], self.df['alt'] = zip(*self.df['change'].apply(lambda x: parse_snv(x) if re.match(r'[A-Z]>[A-Z]$', x) else
                                                      parse_del(x) if re.match(r'del[A-Z]+$', x) else
                                                      parse_ins(x) if re.match(r'ins[A-Z]+$', x) else
                                                      parse_delins(x) if re.match(r'del[A-Z]+ins[A-Z]+$', x) else
                                                      parse_other(x)))
    ## Filtering the data and creating the file in vcf format    
    def filter_data_createvcf(self):
        # Filter rows lacking clinical significance, phenotype/disease, phenotype/disease source information
        # Filter df if a row has at least one hyphen in ref and alt column
        # This one filters genomic alterations: delNN, dupNN
        # There will be only SNVs and delNNinsNN type alterations.     
        ref_alt_filtered = self.df[~((self.df['ref'] == '-') | (self.df['alt'] == '-'))]
            
        values_to_filter = ['R', 'Y', 'D', 'H']
        # This one filters ambiguous bases from alt column because vcf format does not accept them
        final_filtered = ref_alt_filtered[~ref_alt_filtered['alt'].isin(values_to_filter)]
        #This one filters '-' chromosome names
        final_filtered_chr = final_filtered.loc[~(final_filtered['Chromosome'] == '-'), :]
        

        # A function to format strings
        def format_element(element):
            if pd.notna(element):
                element = element.replace(' ', '_')
                if ',' in element:
                    return element.replace(',_', '&')
                if ';' in element:
                    return element.replace(';', '')
                return element
            return ''
        # Create info column
        final_filtered_chr['info'] = (
                'SYMBOL=' + final_filtered_chr['Gene Name'].apply(format_element) + ';' +
                'VAC=' + final_filtered_chr['Variant AA Change'].apply(format_element) + ';' +
                'CT=' + final_filtered_chr['Consequence Type'].apply(format_element) + ';' +
                'CLINSIG=' + final_filtered_chr['Clinical Significance'].apply(format_element) + ';' +
                'PD=' + final_filtered_chr['Phenotype/Disease'].apply(format_element) + ';' +
                'PDS=' + final_filtered_chr['Phenotype/Disease Source'].apply(format_element)
        )
        #Create the vcf file
        final_filtered_chr.reset_index(drop=True, inplace=True)
        self.vcf_df = pd.DataFrame({
            'CHROM': final_filtered_chr['Chromosome'],
            'POS': final_filtered_chr['genomic_position'],
            'ID': final_filtered_chr['Source DB ID'],
            'REF': final_filtered_chr['ref'],
            'ALT': final_filtered_chr['alt'],
            'QUAL': pd.Series('30', index=range(len(final_filtered_chr))),
            'FILTER': pd.Series('PASS', index=range(len(final_filtered_chr))),
            'INFO': final_filtered_chr['info'],
            'FORMAT': pd.Series('GT:AD:AF:DP', index=range(len(final_filtered_chr))),
            'UNIPROT': pd.Series('1/1:0,10:1.0000:10', index=range(len(final_filtered_chr)))
        })
            
    ## Exporting vcf

    def export_vcf_file(self):
        header1 = "##fileformat=VCFv4.2"
        header2 = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUNIPROT"

        #vcf_df_final_wt = pd.concat([pd.Series([header1]), pd.Series([header2]), self.vcf_df])

        with open('vcf_df_final2.vcf', 'w') as file:
            file.write(header1 + '\n')
            file.write(header2 + '\n')
            self.vcf_df.to_csv(file, sep='\t', index=False, header=False)

        os.system(f"sort -k1,1 -k2,2n {'vcf_df_final2.vcf'} > {'vcf_df_final2_sorted.vcf'}")
        os.system(f"bgzip -c {'vcf_df_final2_sorted.vcf'} > {'vcf_df_final2_sorted.vcf.gz'}")
        os.system(f"tabix -p vcf {'vcf_df_final2_sorted.vcf.gz'}")


variation_file_path = '/home/genwork2/Mert/uniprot_vcf/homo_sapiens_variation.txt.gz'

chr_id_path = '/home/genwork2/Mert/Chralias/ultimate_chr_keys_merged.txt'

converter = VCFGenerator(variation_file_path, chr_id_path)


converter.read_variation_file()
converter.read_chrmap_file()
converter.process_chromosome_names()
converter.get_genomic_positions()
converter.parse_genomic_alterations()
converter.filter_data_createvcf()
converter.export_vcf_file()


