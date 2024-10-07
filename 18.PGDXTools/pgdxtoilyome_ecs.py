import os
import json

class CsvToJsonConverterECS:

    """Converts ECS_EPC_RUO.csv files into ilyome compatible json files."""

    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.output_folder = output_folder

    def convert_csv_to_json(self, csv_filepath):
        with open(csv_filepath, "r") as csv_file:
            csv_content = csv_file.read()

        sections = [section.strip() for section in csv_content.strip().split('\n\n')]

        result = {}
        for section in sections:
            lines = section.split('\n')
            section_name = lines[0].strip('[]')
            section_header = lines[1].split(',')
            csv_lines = lines[2:]

            section_data = {}

            for item in csv_lines:
                metric_dict = {}
                if len(section_header) == 2:
                    metric, value = item.split(',')
                    metric_dict[metric] = value
    
                if len(section_header) > 2:
                    metric_lines = item.split(',')
                    metric = metric_lines[0]
                    metric_dict[metric] = {}
                    for i in range(1, len(section_header)):
                        metric_dict[metric][section_header[i]] = metric_lines[i]
        
                section_data.update(metric_dict)

            result[section_name] = section_data

        json_filename = os.path.splitext(csv_filepath)[0].split("/")[-1] + ".json"
        json_filepath = os.path.join(self.output_folder, json_filename)

        with open(json_filepath, "w") as json_file:
            json.dump(result, json_file, indent=2)

    def convert_all_csv_files(self):
        for root, dirs, files in os.walk(self.input_folder):
            for filename in files:
                if filename.endswith("ECS_EPC_RUO.csv"):
                    csv_filepath = os.path.join(root, filename)
                    self.convert_csv_to_json(csv_filepath)


if __name__ == "__main__":
    
    input_folder_path = "/home/genwork2/Mert/pgdx"
    output_folder_path = "/home/genwork2/Mert/md5sumchecker"
    converter = CsvToJsonConverterECS(input_folder_path, output_folder_path)
    converter.convert_all_csv_files()



