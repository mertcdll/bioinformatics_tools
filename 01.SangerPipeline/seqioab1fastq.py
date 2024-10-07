from Bio import SeqIO

records = SeqIO.parse("/home/genwork2/Mert/ab1fastq/IFNAR1-T2D-2_001_A01.ab1", "abi")
count = SeqIO.write(records, "/home/genwork2/Mert/ab1fastq/IFNAR1-T2D-2_001_A01.fastq", "fastq")