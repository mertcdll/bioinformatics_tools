import os
import csv
import pandas as pd

directory = '/home/genwork2/Mert/bwa_dedup/genoks-re-test'

dup_files = []

for root, dirs, files in os.walk(directory):
    for file in files:
        if file.endswith('_dup_metrics.txt'):
            dup_files.append(os.path.join(root, file))
 

headers_written = False

with open("/home/genwork2/Mert/bwa_dedup/genoks-re-test/dup_metrics.csv", 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    
    for file in dup_files:
        filename = os.path.basename(file)
        sample_id = filename.replace('_dup_metrics.txt', '')

        with open(file, 'r') as dup_file:
            lines = dup_file.readlines()

        start_idx = None
        for idx, line in enumerate(lines):
            if line.startswith('## METRICS CLASS'):
                start_idx = idx + 1
                break

        if start_idx is None:
            raise ValueError("Metrics class not found in the file.")

        headers = lines[start_idx].strip().split('\t')
        data = lines[start_idx + 1].strip().split('\t')

        # Add sample_id to the headers and data
        if not headers_written:
            headers.insert(0, 'sample_id')
            writer.writerow(headers)
            headers_written = True
        
        data.insert(0, sample_id)
        writer.writerow(data)
        
        
        
df = pd.read_csv("/home/genwork2/Mert/bwa_dedup/genoks-re-test/dup_metrics.csv")


df = df.sort_values(by="sample_id")


df.to_csv("/home/genwork2/Mert/bwa_dedup/genoks-re-test/dup_metrics.csv", index=False)