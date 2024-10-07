import os
import argparse
import pandas as pd

class ABItoFastq:
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.output_folder = output_folder

    def create_pfq_tsv(self):
        
        for filename in os.listdir(self.input_folder):
            if filename.endswith("ab1"):
                ab1_path = os.path.join(self.input_folder, filename)
                ab1_file_name = os.path.splitext(filename)[0]
                out_dir = os.path.join(self.output_folder, "basecalls", ab1_file_name)
                os.makedirs(out_dir, exist_ok=True)
                output_fastq = os.path.join(out_dir, f"{ab1_file_name}.fastq")
                output_tsv = os.path.join(out_dir, f"{ab1_file_name}.tsv")

                command_pfq = f"tracy basecall -f fastq -o {output_fastq} {ab1_path}"
                command_tsv = f"tracy basecall -f tsv -o {output_tsv} {ab1_path}"

                os.system(command_pfq)
                os.system(command_tsv)

    def create_fastqs(self):
        
        for root, dirs, files in os.walk(self.output_folder):
            for filename in files:
                
                if filename.endswith(".tsv"):
                  filename_wo_ext = os.path.splitext(filename)[0]
                  tsv_path  = os.path.join(root, f"{filename_wo_ext}.tsv")
                  fastq_path = os.path.join(root, f"{filename_wo_ext}.fastq")

                  df = pd.read_csv(tsv_path, sep="\t")
                  
                  with open(fastq_path) as file:
                      primary_fastq = file.read()

                  filtered_df = df.dropna(subset=['primary', 'secondary'])

                  secondary_list = list(filtered_df["secondary"])
                  secondary_seq = ''.join(secondary_list)

                  consensus_list = list(filtered_df["consensus"])
                  consensus_seq = ''.join(consensus_list)

                  quality_scores = primary_fastq.split("\n")[3]

                  secondary_fastq = "@secondary" + "\n" + secondary_seq + "\n" + "+" + "\n" + quality_scores + "\n"
                  consensus_fastq = "@consensus" + "\n" + consensus_seq + "\n" + "+" + "\n" + quality_scores + "\n"

                  primary_secondary_fastq = primary_fastq + secondary_fastq

                  with open(os.path.join(root, f"{filename_wo_ext}_ps.fastq"), "w") as output_file:
                      output_file.write(primary_secondary_fastq)

                  with open(os.path.join(root, f"{filename_wo_ext}_cons.fastq"), "w") as output_file:
                      output_file.write(consensus_fastq)
                  
                  filtered_df.to_csv(os.path.join(root, f"{filename_wo_ext}_filtered.tsv"), sep = "\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert ABI files to Fastq format')
    parser.add_argument('-i', '--input_folder', type=str, required=True, help='Path to the ABI file')
    parser.add_argument('-o', '--output_folder', type=str, required=True, help='Output folder path')

    args = parser.parse_args()

    input_folder = args.input_folder
    output_folder = args.output_folder

    ab12fastq = ABItoFastq(input_folder, output_folder)
    ab12fastq.create_pfq_tsv()
    ab12fastq.create_fastqs()