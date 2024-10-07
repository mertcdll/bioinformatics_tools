import requests
import pandas as pd
import os 

class ClingenDataProcessor:
    """A class for processing clingen data."""

    def __init__(self):
        """Initializes ClingenDataProcessor class"""
        self.gd_df = None
        self.idmaps = None

    def download_csv(self, csv_url, save_path):
        """Downloads and writes gene-disease-validity csv file"""
        with open(save_path, 'wb') as f, requests.get(csv_url, stream=True) as r:
            for line in r.iter_lines():
                f.write(line + b'\n')

    def import_data(self, csv_path, idmap_path):
        """Imports downloaded csv file and ip map file which contains hgnc ids and their corresponding ncbi ids"""
        # Import the csv file
        gd_df = pd.read_csv(csv_path, sep=",", header=4)
        gd_df = gd_df.loc[1:,]

        # Import the id map file
        idmaps = pd.read_csv(idmap_path, sep="\t", dtype=str)

        self.gd_df = gd_df
        self.idmaps = idmaps

    def process_data(self):
        """Maps ncbi ids to hgnc ids and makes column names lower case."""
        def map_ids(row):
            hgnc_id = row["GENE ID (HGNC)"]
            match = self.idmaps[self.idmaps["HGNC ID"] == hgnc_id]

            if not match.empty:
                return match["NCBI Gene ID"].iloc[0]
            else:
                return None

        self.gd_df["GENE ID (HGNC)"] = self.gd_df.apply(map_ids, axis=1)
        self.gd_df.rename(columns={"GENE ID (HGNC)": "GENE ID (NCBI)"}, inplace=True)
        self.gd_df.columns = self.gd_df.columns.str.lower()

    def export_csv(self, output_dir):
        """Exports processed csv file"""
        self.gd_df.to_csv(os.path.join(output_dir, "gene-disease-validity_ncbiids.csv"), index=False)

# Create an instance
processor = ClingenDataProcessor()

# URLs and file paths:
csv_url = 'https://search.clinicalgenome.org/kb/gene-validity/download'
csv_path = "/home/genwork2/Mert/clingen/gene-disease_validity.csv"
idmap_path = "/home/genwork2/Mert/clingen/idmap.txt"
output_dir = "/home/genwork2/Mert/clingen"

# Download the csv file
processor.download_csv(csv_url, csv_path)

# Import the csv and id map files
processor.import_data(csv_path, idmap_path)

# Process the data
processor.process_data()

# Export the final DataFrame as a CSV file
processor.export_csv(output_dir)
