ParseGTFAttributes - A Phyton Class for Parsing GTF files and Creating BED files

ParseGTFAttributes is a Python class that enables the parsing of GTF (Gene Transfer Format) 
files, extracting attributes, and creating BED (Browser Extensible Data) files containing
gene names and their locations. This class has two main functions:

1. ReadParseGTF: This function parses any GTF file and extracts attribute columns, adding them
as new columns to the DataFrame. The parsed DataFrame is then saved as a new file in 
tab-separated format.

2. create_bed_file: This function creates a BED file from the parsed DataFrame containing gene
names and their genomic locations. The function also handles if chromosome aliases found in 
GTF files.

# Installation

Phyton 3.x and the required libraries ('pandas', 'os' and 'gzip') must be installed. 

# Usage

1. Import necessary libraries:

import pandas as pd
import os
import gzip


2. Create an instance of the ParseGTFAttributes class by providing the url, path to 
the folder at which GTF file is downloaded and the chromosome key file (if necessary), and create a new parsed file:

Usage - GTF files from NCBI

url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz'
folder_path = "/path/to/folder"
chr_key_path = "/path/chr_key_file.txt" 
parser = ParseGTFAttributes(url, folder_path, chr_key_path)
file_path = parser.download_gtf_file()
output_file_path = parser.ReadParseGTF()
print("GTF file downloaded at:", file_path)
print("Parsed file saved at:", output_file_path)

Usage - GTF files from ENSEMBL

url = "https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
folder_path = "/path/to/folder"
chr_key_path = "/path/chr_key_file.txt" 
parser = ParseGTFAttributes(url, folder_path, chr_key_path)
file_path = parser.download_gtf_file()
output_file_path = parser.ReadParseGTF()
print("GTF file downloaded at:", file_path)
print("Parsed file saved at:", output_file_path)

3. Call the create_bed_file function to create a BED file containing gene names 
and their genomic locations:

Usage - GTF files from NCBI

bed_file_path = parser.create_bed_file(data_from_ncbi=True)
print("BED file saved at:", bed_file_path)

Usage - GTF files from ENSEMBL

bed_file_path = parser.create_bed_file(data_from_ncbi=False)
print("BED file saved at:", bed_file_path)

# Input files

1. GTF File (file_path): Provide the path to the GTF file in GZIP format. The GTF file must follow 
the standard GTF format.

2. Chromosome Key File (chr_key_path): Provide the path to the chromosome key file, which maps 
chromosome aliases to their actual names. This is only necessary for handling chromosome aliases that 
some GTF files may use.

Chromosome key file format should be:

alias	chr_name
NC_000001.11	chr1
NC_000002.12	chr2
NC_000003.12	chr3
NC_000004.12	chr4
NC_000005.10	chr5
.
.
.

# Output files

1. Parsed GTF File: The ReadParseGTF function creates a new parsed file in tab-separated format. 
The file name will be the same as the input GTF file but with _parsed.txt appended to it.

2. BED File: The create_bed_file function creates a BED file containing gene names and their genomic 
locations. The file name will be the same as the input GTF file but with .gene_names.bed appended 
to it. The BED file will also be sorted and compressed using bgzip and indexed using tabix. 

# Warnings

If this phyton class will be used with a chromosome key file, process may take quite longer. 
So, be patient. 

# License

This project licensed under GENOKS A.Ş. 

# Acknowledgments

The class ParseGTFAttributes was developed to handle GTF files commonly used in genomics research 
and to extract and process gene-related information from such files efficiently.

