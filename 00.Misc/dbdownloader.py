import requests
import pandas as pd
import os

class Downloader:

    def __init__(self, url, output_folder):
        self.url = url
        self.output_folder = output_folder
        self.path = None

    def download_file(self, url, output_folder):
        
        filename = os.path.basename(url)
        file_path = os.path.join(output_folder, filename)
        response = requests.get(url, stream=True)
        with open(file_path, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
        return file_path

    def downloader(self):

        self.path = self.download_file(self.url, self.output_folder)

        return self.path
    
    def tabular(self):

        extension = os.path.splitext(self.path)[-1]

        if extension == ".txt" or extension == ".tsv":
            df = pd.read_csv(self.path, sep="\t")

        elif extension == ".csv":
            df = pd.read_csv(self.path)

        elif extension == ".vcf" or extension == ".gz":
            df = None
        
        tab_dict = {
            "path":self.path,
            "tabular":df,
            }
        
        return tab_dict



if __name__== "__main__":
    url = "https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv"
    output_folder = "/home/genwork2/Mert/dbdownloader"
    down = Downloader(url, output_folder)
    path = down.downloader()
    tabdict = down.tabular()

