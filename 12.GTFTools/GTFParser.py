import pandas as pd
import os
import gzip

class ParseGTFAttributes:
    def __init__(self, file_path, chr_key_path):
        # Initialize the class with file paths for GTF and chromosome key files
        self.file_path = file_path
        self.db = None
        self.chr_key_path = chr_key_path

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

    def create_bed_file(self, use_chr_key = False):
        # This function creates a BED file containing gene names and their locations
        if use_chr_key:

            chr_key = pd.read_csv(self.chr_key_path, sep='\t')

            def extract_chromosome_name(row):
                # This function extracts chromosome names from chromosome aliases using the chromosome key
                chromosome_id = row['seqname']

                match = chr_key[chr_key['ALIAS'] == chromosome_id]
                if not match.empty:
                    return match['NAME'].iloc[0]
                return None

            # Applying the extract_chromosome_name function to create a new column 'chr_name'
            self.db['chr_name'] = self.db.apply(extract_chromosome_name, axis=1)

            # Extracting the bare chromosome names by excluding additional characters after 'chrN'
            self.db['chr_name_bare'] = self.db['chr_name'].str.split('_', expand=True).iloc[:, 0]

            # Selecting only the 'gene' features from the DataFrame
            bed_df = self.db.loc[self.db['feature'] == 'gene', ['chr_name_bare', 'start', 'end', 'gene', 'gene_id']]

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
        bed_df = self.db.loc[self.db['feature'] == 'gene',['seqname', 'start', 'end', 'gene_name', 'gene_id']] #first bed
        bed_df = bed_df.astype(str)
        bed_df.replace('-', pd.NA, inplace=True)
        
        # If there is not a gene name for an entry, add gene id.
        bed_df['gene_name'].fillna(bed_df['gene_id'], inplace=True)
        bed_df.loc[:, 'seqname'] = 'chr' + bed_df['seqname'].astype(str)
        bed_df.drop('gene_id', axis=1, inplace=True)

        bed_df.sort_values(by=['seqname', 'start'], inplace=True)
        
        # Writing the BED file and sorting it
        bed_file_path = self.file_path.rsplit('.', 1)[0] + '.gene_names.bed'
        bed_df.to_csv(bed_file_path, sep="\t", index=False, header=False)
        
        sorted_bed_file_path = bed_file_path.rsplit('.', 1)[0] + '.sorted.bed'
        gzipped_bed_file_path = sorted_bed_file_path + '.gz'
        
        os.system(f"sort -k1,1 -k2,2n {bed_file_path} > {sorted_bed_file_path}") 
        os.system(f"bgzip -c {sorted_bed_file_path} > {gzipped_bed_file_path}")
        os.system(f"tabix -p bed {gzipped_bed_file_path}")


        return gzipped_bed_file_path


# Usage with a chromosome key - GTF files from NCBI

file_path = "/home/genwork2/Mert/GRCh38_latest_genomic.gtf.gz"
chr_key_path = "/home/genwork2/Mert/Chralias/final_ultimate_merged_chrkeys.txt"
 
parser = ParseGTFAttributes(file_path, chr_key_path)
output_file_path = parser.ReadParseGTF()
bed_file_path = parser.create_bed_file(use_chr_key=True)

print("Parsed file saved at:", output_file_path)
print("BED file saved at:", bed_file_path)

# Usage without chromosome key - GTF files from ENSEMBL

file_path = "/home/genwork2/Mert/handle_parse_gtfs/Homo_sapiens.GRCh38.110.gtf.gz"
parser = ParseGTFAttributes(file_path, chr_key_path=None)
output_file_path = parser.ReadParseGTF()
bed_file_path = parser.create_bed_file(use_chr_key=False)

print("Parsed file saved at:", output_file_path)
print("BED file saved at:", bed_file_path)

