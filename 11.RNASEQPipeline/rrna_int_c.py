import pandas as pd
import os
import gzip
import argparse

class rRNAint:
    def __init__(self, bam_file_directory, gtf_path, output_directory):
        self.bam_file_directory = bam_file_directory
        self.gtf_path = gtf_path
        self.output_directory = output_directory
        self.db = None
        self.gtfparser()
        self.convert_rrna_int_list()
    
    def gtfparser(self):

        with gzip.open(self.gtf_path, 'rb') as file:
            df = pd.read_csv(file, delimiter='\t', header=None, comment='#', dtype='str')

        df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

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

    def convert_rrna_int_list(self):
       
        rrna_df = self.db.loc[self.db['transcript_biotype'].str.contains('rRNA'), :]
    
            
        tr_rrna = rrna_df.loc[rrna_df["feature"] == "transcript", :]

        int_list = pd.DataFrame({
            "chrom": tr_rrna["seqname"],
            "start": tr_rrna["start"],
            "end": tr_rrna["end"],
            "strand": tr_rrna["strand"],
            "id": tr_rrna["transcript_id"]
        }) 

        int_list.reset_index(inplace=True, drop=True) # this is the interval list obtained from gtf

        int_list.to_csv(os.path.join(self.output_directory, "intlist.txt"), sep = "\t", header=None, index=False)
        
        os.system(f"samtools view -H {os.path.join(self.bam_file_directory, 'final.bam')} | grep -P '^@SQ' | cut -f 2,3 > {os.path.join(self.output_directory, 'chromosome_lengths.txt')}")
            
        ch_len = pd.read_csv(os.path.join(self.output_directory, "chromosome_lengths.txt"), sep="\t", header=None)

        intlist_header = pd.concat([pd.DataFrame({"sqs": ["@SQ"]* len(ch_len)}), ch_len], axis=1)

            
        header = "@HD\tVN:1.4\tSO:coordinate"
            
        with open(os.path.join(self.output_directory, "int_list_header.txt"), 'w') as file:
            file.write(header + '\n')
            intlist_header.to_csv(file, sep='\t', index=False, header=None)

        os.system(f"cat {os.path.join(self.output_directory, 'int_list_header.txt')} {os.path.join(self.output_directory, 'intlist.txt')} > {os.path.join(self.output_directory, 'rrna_intlist.txt')}")


        os.remove(os.path.join(self.output_directory, "int_list_header.txt"))
        os.remove(os.path.join(self.output_directory, "intlist.txt"))
        os.remove(os.path.join(self.output_directory, "chromosome_lengths.txt"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Creates an rRNA interval list using bam file and GTF file")
    parser.add_argument("--gtf_path", help="Path to GTF file")
    parser.add_argument("--bam_directory", help ="Directory of BAM file")
    parser.add_argument("--output_directory", help="Directory in which rRNA interval list is written")

    args = parser.parse_args()

    convert = rRNAint(args.bam_directory, args.gtf_path, args.output_directory)


