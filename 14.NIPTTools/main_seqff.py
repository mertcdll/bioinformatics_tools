import os
from concurrent.futures import ThreadPoolExecutor

def run_command(command):
    os.system(command)

def seqff_aut(input_folder, output_folder):
        
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(".bam"):
                bam_name = os.path.splitext(file)[0]
                bam_path = os.path.join(root, file)
                output_file = os.path.join(output_folder, f"{bam_name}.txt")
                command = f"Rscript /home/genwork2/Mert/NIPT/seqff.R -f {bam_path} -o {output_file}"
                run_command(command)
                
                
input_folder = "/home/genwork2/Mert/NIPT/NIPT_BAMS_HG19"
output_folder = "/home/genwork2/Mert/NIPT/nipt_metricshg19"
seqff_aut(input_folder, output_folder)
