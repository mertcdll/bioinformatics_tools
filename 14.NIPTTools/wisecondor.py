import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def convert_single_bam(input_file, output_file):
    command = ["WisecondorX", "convert", input_file, output_file]
    
    try:
        subprocess.run(command, check=True)
        return f"Converted {input_file} to {output_file}"
    except subprocess.CalledProcessError as e:
        return f"Error converting {input_file}: {e}"

def convert_bam_to_npz(input_dir, output_dir, max_workers=4):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    bam_files = [
        (os.path.join(input_dir, filename), os.path.join(output_dir, f"{os.path.splitext(filename)[0]}.npz"))
        for filename in os.listdir(input_dir) if filename.endswith(".bam")
    ]
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(convert_single_bam, input_file, output_file): input_file for input_file, output_file in bam_files}
        
        for future in as_completed(futures):
            print(future.result())

input_directory = "/home/genwork2/Mert/NIPT/240902B_NIPT_BAMS"
output_directory = "/home/genwork2/Mert/NIPT/240902B_NIPT_BAMS"

convert_bam_to_npz(input_directory, output_directory, max_workers=40)
