import os
from collections import defaultdict
import subprocess
import gzip

def filter_unique_prefixes(paths):

    prefix_count = defaultdict(int)
    
    for path in paths:
        filename = os.path.basename(path)
        prefix = filename.split("_")[0]
        prefix_count[prefix] += 1
    
    filtered_paths = [path for path in paths if prefix_count[os.path.basename(path).split("_")[0]] > 1]
    
    return filtered_paths

def process_illumina_samples(input_folder_paths, name):
    r1_paths = []
    r2_paths = []
    
    for folder in input_folder_paths:
        for root, dirs, files in os.walk(folder):
            for filename in files:
                if filename.startswith(name) and filename.endswith("R1_001.fastq.gz"):
                    r1_path = os.path.join(root, filename)
                    r1_paths.append(r1_path)
                elif filename.startswith(name) and filename.endswith("R2_001.fastq.gz"):
                    r2_path = os.path.join(root, filename)
                    r2_paths.append(r2_path)
    
    return r1_paths, r2_paths

def process_non_illumina_samples(input_folder_paths, item):
    r1_paths = []
    r2_paths = []
    
    for folder in input_folder_paths:
        for root, dirs, files in os.walk(folder):
            for filename in files:
                wofq = filename.split(".")[0]
                machine_code = wofq.split("_")[0]
                sample_code = wofq.split("_")[-2]
                                
                if item[0].split("_")[0] == sample_code and item[1] == machine_code and filename.endswith("1.fq.gz"):
                    r1_path = os.path.join(root, filename)
                    r1_paths.append(r1_path)
                elif item[0].split("_")[0] == sample_code and item[1] == machine_code and filename.endswith("2.fq.gz"):
                    r2_path = os.path.join(root, filename)
                    r2_paths.append(r2_path)
    
    return r1_paths, r2_paths

def merge_fastq_files(input_paths, output_path):
    """Manually merge the gzipped input fastq files into one output file."""
    with gzip.open(output_path, 'wb') as outfile:
        for file_path in input_paths:
            try:
                with gzip.open(file_path, 'rb') as infile:
          
                    outfile.write(infile.read())
            except Exception as e:
                print(f"An error occurred while reading {file_path}: {e}")

def merger(input_folder_paths, sample_ids, output_folder, illumina):
    
    if illumina:
        
        for name in sample_ids:
            print(f"#################{name}###################")
            
            r1_paths, r2_paths = process_illumina_samples(input_folder_paths, name)
            
            fq1_name = f"{name}_MF_R1_001.fastq.gz"
            fq2_name = f"{name}_MF_R2_001.fastq.gz"
            
            filtered_r1_paths = filter_unique_prefixes(r1_paths)
            filtered_r2_paths = filter_unique_prefixes(r2_paths)
            
            print(f"r1 paths are: {r1_paths}")
            print(f"filtered r1 paths are: {filtered_r1_paths}")
            
            print(f"r2 paths are: {r2_paths}")
            print(f"filtered r2 paths are: {filtered_r2_paths}")
            
            fq1_output_path = os.path.join(output_folder, fq1_name)
            try:
                merge_fastq_files(filtered_r1_paths, fq1_output_path)
                print(f"{fq1_name} is created")
            except Exception as e:
                print(f"Failed to create {fq1_name}: {str(e)}")

            fq2_output_path = os.path.join(output_folder, fq2_name)
            try:
                merge_fastq_files(filtered_r2_paths, fq2_output_path)
                print(f"{fq2_name} is created")
            except Exception as e:
                print(f"Failed to create {fq2_name}: {str(e)}")
        
    else:
        for item in sample_ids.items():
            
            print(f"#################{item[0]}###################")
            
            r1_paths, r2_paths = process_non_illumina_samples(input_folder_paths, item)
            
            fq1_name = f"{item[1]}_MF_{item[0].split('_')[0]}_1.fq.gz"
            fq2_name = f"{item[1]}_MF_{item[0].split('_')[0]}_2.fq.gz"
        
            filtered_r1_paths = filter_unique_prefixes(r1_paths)
            filtered_r2_paths = filter_unique_prefixes(r2_paths)
            
            print(f"r1 paths are: {r1_paths}")
            print(f"filtered r1 paths are: {filtered_r1_paths}")
            
            print(f"r2 paths are: {r2_paths}")
            print(f"filtered r2 paths are: {filtered_r2_paths}")
            
            fq1_output_path = os.path.join(output_folder, fq1_name)
            try:
                merge_fastq_files(filtered_r1_paths, fq1_output_path)
                print(f"{fq1_name} is created")
            except Exception as e:
                print(f"Failed to create {fq1_name}: {str(e)}")

            fq2_output_path = os.path.join(output_folder, fq2_name)
            try:
                merge_fastq_files(filtered_r2_paths, fq2_output_path)
                print(f"{fq2_name} is created")
            except Exception as e:
                print(f"Failed to create {fq2_name}: {str(e)}")
                print(f"An unexpected error occurred while processing {fq2_name}: {str(e)}")

folder_paths = ["/mnt/gen100/01.fastq_files/01.others/09.BilkentSehirPediatri/others/mgi"]

# names = [
#     "78", "82", "80",
#     "91", "99", "60", 
#     "74", "67", "91", 
#     "101", "95", "61"
# ]

names_with_machine_codes = {
    "78": "V300090591",
    "82": "V350031245",
    "80": "V350090653",
    "91": "V350128898",
    "99": "V350151974",
    "60": "V350159139",
    "74": "V350159444",
    "67": "V350166743",
    "91_p": "V350166887", 
    "101": "V350166901",
    "95": "V350194579",
    "61": "V350194626"
}


output = "/home/genwork2/Mert/mergedfq/bilkentpediatri"
illumina = False

merger(folder_paths, names_with_machine_codes, output, illumina)


