import os
import subprocess

def convert_gz_to_bgz(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".gz"):
            gz_file_path = os.path.join(directory, filename)
            bgz_file_path = gz_file_path.replace(".gz", ".bgz")
            
            command = f"zcat {gz_file_path} | bgzip -c > {bgz_file_path}"
            subprocess.run(command, shell=True, check=True)
            print(f"Converted {gz_file_path} to {bgz_file_path}")

directory = '/path/to/your/directory'

convert_gz_to_bgz(directory)