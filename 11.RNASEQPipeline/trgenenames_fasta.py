import re
import pandas as pd
##ncbi

with open ("/home/genwork2/Mert/GRCh38_latest_rna.fna", 'r') as file:
    fastaf = file.read()


sp_fasta = fastaf.split("\n")

headerlist = []

for item in sp_fasta:
    if item.startswith(">"):
        headerlist.append(item)


transcript_ids = [line.split()[0][1:] for line in headerlist]

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

df_ncbi = pd.DataFrame({
    "TXNAME" : transcript_ids,
    "GENEID": gene_ids
})


df_ncbi.to_csv("/home/genwork2/Mert/RNAseq/ncbitrgene.txt",sep="\t", index=False)


##ensembl

with open ("/home/genwork2/Mert/RNAseq/human_transcriptome.fa", 'r') as file:
    fastafn = file.read()

sp_fastan = fastafn.split("\n")

headerlistn = []

for item in sp_fastan:
    if item.startswith(">"):
        headerlistn.append(item)

headerlistn[0]


transcript_idsn = [line.split()[0][1:] for line in headerlistn]


gene_idsn = [re.search(r'gene:([^ ]+)', line).group(1) for line in headerlistn]


df_ensembl = pd.DataFrame({
    "TXNAME" : transcript_idsn,
    "GENEID": gene_idsn,
})


df_ensembl.to_csv("/home/genwork2/Mert/RNAseq/ensembltrgene.txt",sep="\t", index=False)

