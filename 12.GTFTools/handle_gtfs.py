import pandas as pd
import gzip
import os

class ParseGTFAttributes:
    def __init__(self, file_path):
        self.file_path = file_path
        self.db = None
    
    def ReadParseGTF(self):
        with gzip.open(self.file_path, 'rb') as file:
            df = pd.read_csv(file, delimiter='\t', header=None, comment='#', dtype='str')

        df.columns = ['seqname','source','feature','start','end','score','strand','frame','attribute']

        def parse_attribute(attribute):
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

        attribute_dict = df['attribute'].apply(parse_attribute)
        parsed_attributes_df = pd.DataFrame(attribute_dict.tolist())
        parsed_df = pd.concat([df, parsed_attributes_df], axis=1)
        parsed_df.drop(['attribute'], axis=1, inplace=True)
        parsed_df.fillna('-', inplace=True)
        self.db = parsed_df
        
        output_file_path = self.file_path.rsplit('.', 1)[0] + '_parsed.txt'
        parsed_df.to_csv(output_file_path, sep="\t", index=False)
        return output_file_path
    
    def create_bed_file(self):

        
        bed_df = self.db.loc[self.db['feature'] == 'gene',['seqname', 'start', 'end', 'gene_name', 'gene_id']] #first bed
        bed_df = bed_df.astype(str)
        bed_df.replace('-', pd.NA, inplace=True)
        bed_df['gene_name'].fillna(bed_df['gene_id'], inplace=True)
        bed_df.loc[:, 'seqname'] = 'chr' + bed_df['seqname'].astype(str)
        bed_df.drop('gene_id', axis=1, inplace=True)

        bed_df.sort_values(by=['seqname', 'start'], inplace=True)
        
        bed_file_path = self.file_path.rsplit('.', 1)[0] + '.gene_names.bed'
        bed_df.to_csv(bed_file_path, sep="\t", index=False, header=False)
        
        sorted_bed_file_path = bed_file_path.rsplit('.', 1)[0] + '.sorted.bed'
        gzipped_bed_file_path = sorted_bed_file_path + '.gz'
        
        os.system(f"sort -k1,1 -k2,2n {bed_file_path} > {sorted_bed_file_path}") 
        os.system(f"bgzip -c {sorted_bed_file_path} > {gzipped_bed_file_path}")
        os.system(f"tabix -p bed {gzipped_bed_file_path}")

        new_bed_df = self.db.loc[self.db['feature'] == 'gene', ['gene_name', 'gene_id', 'start', 'end']] #second bed
        new_bed_df = new_bed_df.astype(str)
        new_bed_df.replace('-', pd.NA, inplace=True)
        new_bed_df['gene_name'].fillna(new_bed_df['gene_id'], inplace=True)
        new_bed_df['metainfo'] = 'START:' + new_bed_df['start'] + ';END:' + new_bed_df['end']
        new_bed_df.drop('gene_id', axis=1, inplace= True)
        
        new_bed_file_path = self.file_path.rsplit('.', 1)[0] + '.gene_strend.bed'
        new_bed_df.to_csv(new_bed_file_path, sep='\t', index=False, header=False)

        new_sorted_bed_file_path = new_bed_file_path.rsplit('.', 1)[0] + '.sorted.bed'
        new_gzipped_bed_file_path = new_sorted_bed_file_path + '.gz'

        os.system(f"sort -k1,1 -k2,2n {new_bed_file_path} > {new_sorted_bed_file_path}") 
        os.system(f"bgzip -c {new_sorted_bed_file_path} > {new_gzipped_bed_file_path}")
        os.system(f"tabix -p bed {new_gzipped_bed_file_path}")
        

        return gzipped_bed_file_path , new_gzipped_bed_file_path


file_path = "/home/genwork2/Mert/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
parser = ParseGTFAttributes(file_path)
output_file_path = parser.ReadParseGTF()
bed_file_path = parser.create_bed_file()

print("Parsed file saved at:", output_file_path)
print("BED files saved at:", bed_file_path[0],',',bed_file_path[1])
