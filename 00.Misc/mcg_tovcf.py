"""import pandas as pd

df = pd.read_csv("/home/genwork2/Mert/VariantDatabases/12.MyCancerGenome(ProjectGenie)/data_mutations_extended.txt", sep ="\t")

df.columns
df["NCBI_Build"].value_counts()

df.iloc[123423,]["Tumor_Seq_Allele2"]

df.loc[df["Tumor_Seq_Allele1"] == "-",][["Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]]


df["Tumor_Seq_Allele2"][15035]

sum(df.loc[df["Tumor_Seq_Allele2"] == "-",][["Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]]["Tumor_Seq_Allele1"].isna())



df.iloc[20]["Tumor_Seq_Allele1"]

sum(df["Tumor_Seq_Allele2"] == "nan")


sum(df["Tumor_Seq_Allele1"].isna())

df.isna()


search_value = "ENST00000256078.4:c.34G>T"

"""


#################################################3


import csv
import re
import os

def convert_vcf(row):
    # Extract data from the dictionary
    chrom = "chr" + row["Chromosome"]
    pos = row["Start_Position"]

    id = row["dbSNP_RS"]
    if id is None or id == "":
        id = "."
    else:
        id = row["dbSNP_RS"]

    ref = row.get("Reference_Allele", "") 
    if not ref or ref == "-":
        return None
    
    #a1 = row["Tumor_Seq_Allele1"]
    alt1 = row["Tumor_Seq_Allele1"]
    alt2 = row["Tumor_Seq_Allele2"]

    
    qual = "."
    filter = row["FILTER"]
    
    columns = ["Chromosome", "Start_Position", "dbSNP_RS", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Score", "FILTER"]

    info_cols = [key for key in row.keys() if key not in columns]
    info = ";".join([f"{item}={row[item]}" for item in info_cols])

    vcf_linez = []

    if alt1 and alt1 != "-":
        if ref != alt1:
            vcf_line_alt1 = f"{chrom}\t{pos}\t{id}\t{ref}\t{alt1}\t{qual}\t{filter}\t{info}"
            vcf_lines.append(vcf_line_alt1)

    if alt2 and alt2 != "-":
        if ref != alt2:
            vcf_line_alt2 = f"{chrom}\t{pos}\t{id}\t{ref}\t{alt2}\t{qual}\t{filter}\t{info}"
            vcf_lines.append(vcf_line_alt2)

    return vcf_linez


with open("/home/genwork2/Mert/VariantDatabases/12.MyCancerGenome(ProjectGenie)/data_mutations_extended.txt", "r") as input_file:
    vcf_lines = []
    csv_reader = csv.DictReader(input_file, delimiter='\t')
    for row in csv_reader:
        vcf_line = convert_vcf(row)
        if vcf_line is not None:
            vcf_lines.extend(vcf_line)

vcf_lines_un = list(set(vcf_lines))

with open("/home/genwork2/Mert/VariantDatabases/12.MyCancerGenome(ProjectGenie)/filtered_my_cancer_genome.vcf", "w") as output_file:
    output_file.write("##fileformat=VCFv4.3\n")

    output_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for vcf_line in vcf_lines_un:
        if vcf_line != []:
            output_file.write(vcf_line + "\n")



os.system("sort -k1,1 -k2,2n '/home/genwork2/Mert/VariantDatabases/12.MyCancerGenome(ProjectGenie)/filtered_my_cancer_genome.vcf' > '/home/genwork2/Mert/VariantDatabases/12.MyCancerGenome(ProjectGenie)/filtered_sorted_my_cancer_genome.vcf'")
os.system("bgzip -c '/home/genwork2/Mert/VariantDatabases/12.MyCancerGenome(ProjectGenie)/filtered_sorted_my_cancer_genome.vcf' > '/home/genwork2/Mert/VariantDatabases/12.MyCancerGenome(ProjectGenie)/filtered_sorted_my_cancer_genome.vcf.gz'")
os.system("tabix -p vcf '/home/genwork2/Mert/VariantDatabases/12.MyCancerGenome(ProjectGenie)/filtered_sorted_my_cancer_genome.vcf.gz'")



