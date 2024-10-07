import pandas as pd
import numpy as np
import os


class CSVtoVCFConverter:
    
    """Converts supplemental_filtered_sequence_mutations_EPC-RUO.csv files into ilyome compatible vcf files."""

    def __init__(self, input_folder, output_folder, ref_fasta):
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.ref_fasta = ref_fasta
        self.process_csv_files()

    def process_csv_files(self):
        for root, dirs, files in os.walk(self.input_folder):
            for file in files:
                if file.endswith("supplemental_filtered_sequence_mutations_EPC-RUO.csv"):
                    file_path = os.path.join(root, file)
                   
                    df = pd.read_csv(file_path, sep=",")
                    df.replace({np.nan: ''}, inplace=True)
                    df.replace(';', '_', regex= True, inplace=True)

                    vcf_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'PGDX']

                    vcf_df = pd.DataFrame(columns=vcf_columns)

                    vcf_df['CHROM'] = df['Nucleotide Position (Genomic hg19)'].str.split(":", expand=True)[0]
                    vcf_df['POS'] = df['Nucleotide Position (Genomic hg19)'].str.split(":", expand=True)[1]
                    vcf_df['ID'] = df['dbSNP RefSNP ID (version 150)'].apply(lambda x: x.split("_")[0] if x else '.')

                    if sum(df["Reference"] == "N") != 0: # FixVcfHeader only accepts this character to fill in, and it checks if any variant has N as a reference allele.
                        raise ValueError("Do not fill 'Reference' for insertion type mutations -- It can create ambiguities")
                    else:
                        pass

                    vcf_df['REF'] = df['Reference'].apply(lambda x: "N" if x =="(null)" else x)
                    vcf_df['ALT'] = df['Mutant'].apply(lambda x: "." if x =="(null)" else x)
                    vcf_df['QUAL'] = 30
                    vcf_df['FILTER'] = df['Overall Sample QC']
                    vcf_df['INFO'] = None
                    vcf_df['FORMAT'] = "GT:AD:AF:DP" 
                    vcf_df['PGDX'] = "./.:0,10:1.0000:10"

                    vcf_columns = ['Nucleotide Position (Genomic hg19)','dbSNP RefSNP ID (version 150)','Reference', 'Mutant', 'Overall Run QC','Overall Sample QC']

                    info_columns = df.columns.difference(vcf_columns)

                    vcf_df['INFO'] = df[info_columns].apply(lambda row: ';'.join(f"{str.lower(col)}={val}" for col, val in row.items()), axis=1)

                    vcf_df['INFO'] = vcf_df['INFO'].str.replace(" ", "_")

                    vcf_df['INFO'] = vcf_df['INFO'].str.replace(">=", "gte_") 

                    vcf_df['INFO'] = vcf_df['INFO'].str.replace(">", "_")

                    vcf_df['INFO'] = vcf_df['INFO'].str.replace("(%)", "pct")

                    vcf_df['INFO'] = vcf_df['INFO'].str.replace("(", "")

                    vcf_df['INFO'] = vcf_df['INFO'].str.replace(")", "")

                    vcf_df['INFO'] = vcf_df['INFO'].str.replace("/", "_")

                    vcf_df['INFO'] = vcf_df['INFO'].str.replace("-", "_")
                    
                    vcf_df['INFO'] = vcf_df['INFO'].str.replace("%", "")
                    
                    vcf_df['CHROM'] = vcf_df['CHROM'].replace({
                        'Chr1': 'chr1',
                        'Chr2': 'chr2',
                        'Chr3': 'chr3',
                        'Chr4': 'chr4',
                        'Chr5': 'chr5',
                        'Chr6': 'chr6',
                        'Chr7': 'chr7',
                        'Chr8': 'chr8',
                        'Chr9': 'chr9',
                        'Chr10': 'chr10',
                        'Chr11': 'chr11',
                        'Chr12': 'chr12',
                        'Chr13': 'chr13',
                        'Chr14': 'chr14',
                        'Chr15': 'chr15',
                        'Chr16': 'chr16',
                        'Chr17': 'chr17',
                        'Chr18': 'chr18',
                        'Chr19': 'chr19',
                        'Chr20': 'chr20',
                        'Chr21': 'chr21',
                        'Chr22': 'chr22',
                        'Chrx': 'chrX',
                        'Chry': 'chrY'
                    })

                    header1 = "##fileformat=VCFv4.2"
                    header2 = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPGDX"

                    filename = os.path.splitext(file)[0]
                    output_name = os.path.join(self.output_folder, filename)
                    output_vcf = output_name + ".vcf"
                    sorted_file = output_name + "_sorted.vcf"
                    zipped_file = output_name + "_sorted.vcf.gz"
                    vcf_wheader = output_name + "_header.vcf.gz"
                    norm_vcf = output_name + "_norm.vcf"
                    final_vcf =  output_name + "_final.vcf"
                    final_vcf_gz = final_vcf + ".gz"                   


                    with open(output_vcf, 'w') as file:
                        file.write(header1 + '\n')
                        file.write(header2 + '\n')
                        vcf_df.to_csv(file, sep='\t', index=False, header=False)

                    os.system(f"sort -k1,1 -k2,2n {output_vcf} > {sorted_file}")

                    with open(sorted_file, 'r') as file:
                        lines = file.readlines()
                        lines.insert(0, header1 + '\n')
                        lines.insert(1, header2 + '\n')

                    with open(sorted_file, 'w') as outfile:
                        outfile.writelines(lines)

                    os.system(f"bgzip -c {sorted_file} > {zipped_file}")
                    os.system(f"tabix -p vcf {zipped_file}")
                    os.system(f"gatk FixVcfHeader -I {zipped_file} -O {vcf_wheader} -R {self.ref_fasta}")
                    os.system(f"bcftools norm {vcf_wheader} -N -f {self.ref_fasta} -c s > {norm_vcf}") # this one changes N's in ref alleles to original nucleotide found in hg19 reference genome (for insertion events)
                    
                    
                    with open(norm_vcf, 'r') as input_file, open(final_vcf, 'w') as output_file:
                        for line in input_file:
                            if line.startswith('#'):
                                output_file.write(line)
                            else:
                                fields = line.strip().split('\t') # Only the mutation types "INS and Homo_INS" in variants has the ref allele as null and normalization has fixed it,
                                                                 #now following code inserts ref allele beside alt allele in order to prevent ambiguities with SNPs.
                                if "mutation_type=Homo_INS" in fields[7] or "mutation_type=INS" in fields[7]:
                                    fields[4] = fields[3] + fields[4]
                                    modified_line = '\t'.join(fields)
                                    output_file.write(modified_line + '\n')
                                else:
                                    output_file.write(line)
                    
                    os.system(f"bgzip -c {final_vcf} > {final_vcf_gz}")
                    os.system(f"tabix -p vcf {final_vcf_gz}")
                    os.system(f"sudo rm {self.output_folder}/*sorted* {self.output_folder}/*header* {self.output_folder}/*norm* {self.output_folder}/*RUO.vcf")


if __name__ == "__main__":
    input_folder = "/home/genwork2/Mert/pgdx"
    output_folder = "/home/genwork2/Mert/pgdx"
    reference_fasta = "/home/genwork2/Mert/pgdx/hg19.fa"
    csv_convert = CSVtoVCFConverter(input_folder, output_folder, reference_fasta)




