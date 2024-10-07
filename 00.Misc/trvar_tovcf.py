import os 

def convert_to_vcf(row):
    chrom = row[0]
    pos = row[2]
    if pos == ".":
        return None #grch38 pos
    id = row[4]
    ref = row[5]
    alt = row[6]
    af = row[7]
    ac = row[8]
    an = row[9]
    hom = row[10]
    het = row[11]
    seq_method = row[12]
    filter = row[13]
    qual = row[14].split("|")[0]
    dp = row[15]
    gene_name = row[16]
    gene_id = row[17]
    feature = row[18]
    feature_id = row[19]
    effect = row[20]
    lof = row[21]
    hgvs_c = row[22]
    hgvs_p = row[23]
    gerp_rs = row[24]
    cadd_phred = row[25]
    sift_pred = row[26]
    polyphen2_pred = row[27]
    af_gnomad_wes = row[28]
    af_gnomad_wgs = row[29]
    gme_af = row[30]
    _1000gp_af = row[31]
    esp_af = row[32]

    info = f"AF={af};AC={ac};AN={an};Hom={hom};Het={het};SequencingMethod={seq_method};DP={dp};GeneName={gene_name};GeneID={gene_id};Feature={feature};FeatureID={feature_id};Effect={effect};LoF={lof};HGVS_C={hgvs_c};HGVS_P={hgvs_p};GERP_RS={gerp_rs};CADD_phred={cadd_phred};SIFT_pred={sift_pred};Polyphen2_HVAR_pred={polyphen2_pred};AF_gnomAD_WES={af_gnomad_wes};AF_gnomAD_WGS={af_gnomad_wgs};GME_AF={gme_af};1000GP_AF={_1000gp_af};ESP_AF={esp_af}"

    vcf_line = f"{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}"
    return vcf_line

with open("/home/genwork2/Mert/VariantDatabases/1.TurkishVariome/30739108.tsv", "r") as input_file:
    vcf_lines = []
    for line in input_file:
        row = line.strip().split("\t")
        vcf_line = convert_to_vcf(row)
        if vcf_line is not None:
            vcf_lines.append(vcf_line)

with open("/home/genwork2/Mert/VariantDatabases/1.TurkishVariome/filtered_turkish_variome_grch38.vcf", "w") as output_file:
    output_file.write("##fileformat=VCFv4.3\n")

    output_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for vcf_line in vcf_lines[1:]:
        output_file.write(vcf_line + "\n")


os.system("sort -k1,1 -k2,2n '/home/genwork2/Mert/VariantDatabases/1.TurkishVariome/filtered_turkish_variome_grch38.vcf' > '/home/genwork2/Mert/VariantDatabases/1.TurkishVariome/filtered_sorted_turkish_variome_grch38.vcf'")
os.system("bgzip -c '/home/genwork2/Mert/VariantDatabases/1.TurkishVariome/filtered_sorted_turkish_variome_grch38.vcf' > '/home/genwork2/Mert/VariantDatabases/1.TurkishVariome/filtered_sorted_turkish_variome_grch38.vcf.gz'")
os.system("tabix -p vcf '/home/genwork2/Mert/VariantDatabases/1.TurkishVariome/filtered_sorted_turkish_variome_grch38.vcf.gz'")




####################

import gzip
def get_nucleotide(fasta, chr_name, pos):


    with gzip.open("/home/genwork2/Mert/VariantDatabases/1.TurkishVariome/hg38.fa.gz", "rt") as file:
        for line in file:
            if line.startswith('>'):
                current_chromosome  = line.strip('>')
            elif current_chromosome == chr_name:
                sequence = line.strip()
                


