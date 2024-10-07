import os
import pandas as pd
import argparse

class refflat:
    def __init__(self, gtf_file, output_directory, data_source, gtftogenepred_path):
        self.gtf_file = gtf_file
        self.output_directory = output_directory
        self.data_source = data_source
        self.gtftogenepred_path = gtftogenepred_path
        self.converter()

    def converter(self):
        if self.data_source == "ncbi":
            os.system(f"{self.gtftogenepred_path} {self.gtf_file} {os.path.join(self.output_directory, 'refflat.gp')} -includeVersion -ignoreGroupsWithoutExons -genePredExt")
            gp_df = pd.read_csv(f"{os.path.join(self.output_directory, 'refflat.gp')}", sep='\t', header=None)
            gene_col = gp_df.iloc[:, 11]
            gp_df.insert(0, '', gene_col)
            final_df = gp_df.iloc[:, 0:11]
            final_df.to_csv(f"{os.path.join(self.output_directory, 'refflat.gp')}", sep='\t', header=None, index=False)

        if self.data_source == "ensembl":
            os.system(f"{self.gtftogenepred_path} {self.gtf_file} {os.path.join(self.output_directory, 'refflat.gp')} -ignoreGroupsWithoutExons -genePredExt")
            gp_df = pd.read_csv(f"{os.path.join(self.output_directory, 'refflat.gp')}", sep='\t', header=None)
            gene_col = gp_df.iloc[:, 11]
            gp_df.insert(0, '', gene_col)
            final_df = gp_df.iloc[:, 0:11]
            final_df.to_csv(f"{os.path.join(self.output_directory, 'refflat.gp')}", sep='\t', header=None, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts GTF files from ncbi and emsembl into refFlat format")
    parser.add_argument("--gtf_path", help="Path to GTF file")
    parser.add_argument("--output_path", help ="The path in which refFlat file generated")
    parser.add_argument("--data_source", help="Database that GTF file is obtained (Can be either ensembl or ncbi)")
    parser.add_argument("--program_path", help="Path to gtfToGenePred program")

    args = parser.parse_args()

    convert = refflat(args.gtf_path, args.output_path, args.data_source, args.program_path)