import pandas as pd
import os

directory = '/home/genwork2/Mert/bwa_dedup/genoks-re-test' 

csv_files = [f for f in os.listdir(directory) if f.endswith('_hs_metrics.csv')]

df_list = []

for file in csv_files:
    file_path = os.path.join(directory, file)
    df = pd.read_csv(file_path)
    
    sample_name = file.replace('_hs_metrics.csv', '')
    
    df.insert(0, 'sample_id', sample_name)
    
    df_list.append(df)

combined_df = pd.concat(df_list, ignore_index=True)

combined_df = combined_df.sort_values(by='sample_id')

combined_df.to_csv(os.path.join(directory, 'combined_metrics.csv'), index=False)

print("Concatenation complete. Combined file saved as 'combined_metrics.csv'.")