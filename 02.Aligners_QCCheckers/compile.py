import os
import pandas as pd

class MetricsAggregator:
    def __init__(self, input_dir, output_file):
        self.input_dir = input_dir
        self.output_file = output_file

    def extract_sample_id(self, filename):
        pruned = filename.replace("_metrics.txt", "")
        pruned2 = pruned.replace("sorted_filtered_", "")
        return pruned2

    def process_files(self):
        data = []
        metrics_set = set()

        for root, dirs, files in os.walk(self.input_dir):
            for file in files:
                if file.endswith("_metrics.txt"):
                    file_path = os.path.join(root, file)
                    sample_id = self.extract_sample_id(file)

                    with open(file_path, 'r') as f:
                        metrics = {}
                        metrics['sample_id'] = sample_id
                        for line in f:
                            key, value = line.strip().rsplit(' ', 1)
                            metrics[key] = value
                            metrics_set.add(key)
                        data.append(metrics)

        df = pd.DataFrame(data)
        
        df = df[['sample_id'] + sorted(metrics_set)]
        
        df = df.sort_values(by='sample_id')

        df.to_csv(self.output_file, index=False)

input_directory = '/home/genwork2/Mert/NIPT/nipt_metricshg19/mosdepth_metrics'
output_csv_file = '/home/genwork2/Mert/NIPT/nipt_metricshg19/combined_nipt_metrics_mosdepth.csv'

metrics_aggregator = MetricsAggregator(input_directory, output_csv_file)
metrics_aggregator.process_files()
