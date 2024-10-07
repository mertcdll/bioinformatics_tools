import os
import pandas as pd
import json

class CsvToJsonConverterbTMB:

    """Converts bTMB_EPC_RUO.csv files into ilyome compatible json files."""

    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.output_folder = output_folder

    def process_csv_files(self):
        for root, dirs, files in os.walk(self.input_folder):
            for file in files:
                if file.endswith("bTMB_EPC-RUO.csv"):
                    file_path = os.path.join(root, file)

                    df = pd.read_csv(file_path, sep=",")

                    data_list = df.to_dict(orient='records')

                    json_file_name = os.path.splitext(file)[0] + ".json"
                    output_file_path = os.path.join(self.output_folder, json_file_name)

                    with open(output_file_path, 'w') as json_file:
                        json.dump(data_list, json_file, indent=2)

if __name__ == "__main__":

    input_folder_path = '/home/genwork2/Mert/pgdx'
    output_folder_path = '/home/genwork2/Mert/md5sumchecker'

    csv_processor = CsvToJsonConverterbTMB(input_folder_path, output_folder_path)

    csv_processor.process_csv_files()
