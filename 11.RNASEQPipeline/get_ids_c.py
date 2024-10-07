import re
import pandas as pd
import argparse
import os
import gzip


class GetIDs:
    def __init__ (self, fasta_file, data_source, output_directory):
        self.fasta_file = fasta_file
        self.data_source = data_source
        self.output_directory = output_directory
        self.extractids()
    
    def extractids(self):

        if self.fasta_file.endswith(".gz"):
            with gzip.open(self.fasta_file, 'rt') as file:
                fastaf = file.read()
        else:
            with open (self.fasta_file, 'r') as file:
                fastaf = file.read()

        sp_fasta = fastaf.split("\n")

        headerlist = []

        for item in sp_fasta:
            if item.startswith(">"):
                headerlist.append(item)


        transcript_ids = [line.split()[0][1:] for line in headerlist]


        if self.data_source == "ncbi":
            gene_ids = []
            for item in headerlist:
                matches = re.findall(r'\(([^)]+)\),', item)
    
                if len(matches) == 1:
                    gene_ids.append(matches[0])
                elif len(matches) == 2:
                    second_match = matches[1]
                    if '*' in second_match or ':' in second_match:
                        gene_ids.append(matches[0])
                    else:
                        gene_ids.append(second_match)
                else:
                    gene_ids.append('')

        elif self.data_source == "ensembl":

            gene_ids = [re.search(r'gene:([^ ]+)', line).group(1) for line in headerlist]


        df = pd.DataFrame({
            "TXNAME" : transcript_ids,
            "GENEID": gene_ids
        })


        df.to_csv(os.path.join(self.output_directory, "tx2gene.txt"),sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gets Transcript IDs and their corresponding Gene IDs from transcriptome fasta files")
    parser.add_argument("--fasta_file", help="Path to fasta file")
    parser.add_argument("--data_source", help="Database that GTF file is obtained (Can be either ensembl or ncbi)")
    parser.add_argument("--output_directory", help ="The path in which refFlat file generated")

    args = parser.parse_args()

    convert = GetIDs(args.fasta_file, args.data_source, args.output_directory)

