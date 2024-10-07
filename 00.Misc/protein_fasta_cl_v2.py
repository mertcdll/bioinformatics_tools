import pandas as pd
from Bio import SeqIO
import gzip
import os
import csv
import requests

class ProteinFastaConverter:
    def __init__(self, url, folder_path, chr_alias_file):
        self.url = url
        self.fasta_file = folder_path + '/' + url.split('/')[-1]
        self.chr_alias_file = chr_alias_file
        self.prot_info = None
        self.chr_start_end = None
        self.alias_df = None
        self.bed_df_wov = None
        self.bed_df_woid = None
    
    def download_fasta_file(self):
        response = requests.get(self.url, stream=True)
        with open(self.fasta_file, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)

    def read_fasta_convertdf(self):
        # Read protein information from a compressed fasta file and convert to DataFrame
        
        data = []
        with gzip.open(self.fasta_file, "rt") as gz_handle:
            for record in SeqIO.parse(gz_handle, "fasta"):
                parts = record.description.split()
                # Extract relevant information from the fasta record description
                protein_id = parts[0]
                chr_number = parts[2].split(':')[2]
                start_pos = parts[2].split(':')[3]
                end_pos = parts[2].split(':')[4]
                gene_id = parts[3].split(':')[1]
                transcript_id = parts[4].split(':')[1]
                gene_biotype = parts[5].split(':')[1]
                transcript_biotype = parts[6].split(':')[1]
                gene_symbol_part = [part for part in parts if part.startswith("gene_symbol:")]
                if gene_symbol_part:
                    gene_symbol = gene_symbol_part[0].split(':')[1]
                else:
                    gene_symbol = '-'

                sequence = str(record.seq)
                    
                data.append({
                    "Protein_ID": protein_id,
                    "Chr_number": chr_number,
                    "Start Position": start_pos,
                    "End Position": end_pos,
                    "Gene_ID": gene_id,
                    "Transcript_ID": transcript_id,
                    "Gene_Biotype": gene_biotype,
                    "Transcript_Biotype": transcript_biotype,
                    "Gene_Symbol": gene_symbol,
                    "Sequence": sequence
                })

        # Create a DataFrame from the collected protein information
        self.prot_info = pd.DataFrame(data)

        # Create a DataFrame containing chromosome start and end positions
        self.chr_start_end = self.prot_info[['Chr_number', 'Start Position', 'End Position']]
    
    def read_alias_file(self):
        # Read chromosome alias information
        self.alias_df = pd.read_csv(self.chr_alias_file, sep='\t')

    def process_fasta_file(self):
        # Process protein information to create bed format DataFrames
        
        def extract_chromosome_name(row):
            # Helper function to extract chromosome name using alias mapping
            chromosome_id = row['Chr_number']

            match = self.alias_df[self.alias_df['alias'] == chromosome_id]
            if not match.empty:
                return match['chr_name'].iloc[0]
            return chromosome_id   

        # Add a new column for chromosome names using the alias mapping
        self.chr_start_end['chr_name'] = self.chr_start_end.apply(extract_chromosome_name, axis=1)

        # Remove unnecessary columns from the protein information DataFrame
        self.prot_info.drop(['Chr_number', 'Start Position', 'End Position'], inplace=True, axis=1)

        # Create a copy of the protein information DataFrame with modified ids (version information is dropped)
        prot_info_wov = self.prot_info.copy()
        prot_info_wov['Protein_ID'] = prot_info_wov['Protein_ID'].str.split('.').str[0]
        prot_info_wov['Gene_ID'] = prot_info_wov['Gene_ID'].str.split('.').str[0]
        prot_info_wov['Transcript_ID'] = prot_info_wov['Transcript_ID'].str.split('.').str[0]

        # Convert protein information to json format and create a DataFrame
        prot_info_json_wov = prot_info_wov.apply(lambda row: row.to_json(), axis=1)
        prot_info_json_df = pd.DataFrame({'json_string': prot_info_json_wov})

        # Create a bed format DataFrame without protein or gene or transcript versions.
        self.bed_df_wov = pd.DataFrame({
            'chr': self.chr_start_end['chr_name'],
            'start': self.chr_start_end['Start Position'],
            'end': self.chr_start_end['End Position'],
            'info': prot_info_json_df['json_string']
        })
        
        # Create a copy of the protein information DataFrame without protein id
        protinfo_woid = self.prot_info.drop('Protein_ID', axis=1)

        # Convert protein information without protein id to JSON format and create a DataFrame
        prot_info_json_woid = protinfo_woid.apply(lambda row: row.to_json(), axis=1)
        prot_info_json_woid_df = pd.DataFrame({'json_string': prot_info_json_woid})

        # Create a bed format DataFrame having protein id in its first column and a json string without protein id in its fourth column.
        self.bed_df_woid = pd.DataFrame({
            'protein_id': prot_info_wov['Protein_ID'],
            'start': pd.Series('1', index=range(len(prot_info_wov))),
            'end': pd.Series('4', index=range(len(prot_info_wov))),
            'info': prot_info_json_woid_df['json_string']
        })
    
    def create_bed_file(self):
        # Create and process bed files
        
        with open('proteins.bed', 'w') as file:
            self.bed_df_wov.to_csv(file, sep="\t", index=False, header=False,  quoting=csv.QUOTE_NONE, quotechar="", escapechar="\\")

        os.system(f"sort -k1,1 -k2,2n {'proteins.bed'} > {'proteins_sorted.bed'}")
        os.system(f"bgzip -c {'proteins_sorted.bed'} > {'proteins_sorted.bed.gz'}")
        os.system(f"tabix -p bed {'proteins_sorted.bed.gz'}")
        
        with open('proteins_id.bed', 'w') as file:
            self.bed_df_woid.to_csv(file, sep="\t", index=False, header=False,  quoting=csv.QUOTE_NONE, quotechar="", escapechar="\\")

        os.system(f"sort -k1,1 -k2,2n {'proteins_id.bed'} > {'proteins_id_sorted.bed'}")
        os.system(f"bgzip -c {'proteins_id_sorted.bed'} > {'proteins_id_sorted.bed.gz'}")
        os.system(f"tabix -p bed {'proteins_id_sorted.bed.gz'}")


url = "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
fasta_folder_path = '/home/genwork2/Mert/Downloads'
alias_file_path = '/home/genwork2/Mert/Chralias/final_ultimate_merged_chrkeys.txt'

converter = ProteinFastaConverter(url, fasta_folder_path, alias_file_path)

# Run them in order
converter.download_fasta_file()
converter.read_fasta_convertdf()
converter.read_alias_file()
converter.process_fasta_file()
converter.create_bed_file()
