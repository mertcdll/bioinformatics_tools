import os

base_directory_path = "/mnt/gen100/01.fastq_files/196.TWIST-ExoV2-DNAPrepWithExomePlus-NovaSeq-RUN186_fastq_files/RUN186-PolyA-DNAPrep-Nextera-Ribozero/Gazi-CFTR"

for root, dirs, files in os.walk(base_directory_path):
    folder_name = os.path.basename(root)
        
    for filename in files:
        if "_S" in filename:

            parts = filename.split("_S", 1)
            new_filename = f"{folder_name}_S{parts[1]}"
            old_filepath = os.path.join(root, filename)
            new_filepath = os.path.join(root, new_filename)
            os.rename(old_filepath, new_filepath)
            print(f"Renamed: {old_filepath} -> {new_filepath}")