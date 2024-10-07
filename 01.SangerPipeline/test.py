import subprocess
import os 


def variant_caller(ab1_path, outprefix_path, reference):
    command = ["tracy", "decompose", "-v", "-a", "homo_sapiens", "-r", reference, "-o", outprefix_path, ab1_path]

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
        if stdout:
            print("Command stdout:", stdout)
        if stderr:
            print("Command stderr:", stderr)
    except subprocess.CalledProcessError as e:
        
        print(f"Error executing command: {e}")

  
variant_caller("/home/genwork2/Mert/ab1fastq/2404_GenoksSanger/AYSEGUL_ERDUR_MBD5_8F_B01.ab1", 
               "/home/genwork2/Mert/ab1fastq/outputtest/nermin_aliyeva" , 
               "/home/genwork2/Mert/01.REFGEN/hg38.fa.gz")


















