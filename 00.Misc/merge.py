import os
import logging
import argparse

logging.basicConfig(format='%(asctime)s: %(levelname)s: %(message)s ', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)

def merge_fastq_files_by_sample_id(input_directory, output_directory, sequence_type):
    
    sample_files = {}
    
    def find_samples(input_directory):
        for filename in os.listdir(input_directory):
            if filename.endswith(".fastq.gz"):
                try:
                    if "L00" in filename:
                        sample_id, row, lane, read, suffix = filename.split("_")
                    else:
                        sample_id, row, read, suffix = filename.split("_")
                        
                    sample_id = sample_id[:-3] if sample_id.endswith("-RE") else sample_id
                except Exception as e:
                    print(filename, "Does not match with fastq name format: SampleID_Row_Lane_Read_001.fastq.gz")
                else:
                    if sample_id not in sample_files:
                        sample_files[sample_id] = {"R1": [], "R2":[]}
                    if filename.endswith("R1_001.fastq.gz"):
                        sample_files[sample_id]["R1"].append(os.path.join(input_directory, filename))
                    if filename.endswith("R2_001.fastq.gz"):
                        sample_files[sample_id]["R2"].append(os.path.join(input_directory, filename))
            
            if filename.endswith(".fq.gz"):
                try:
                    sample_id, lane, barcode, read_suffix = filename.split("_")
                    sample_id = sample_id.strip("-RE") if sample_id.endswith("-RE") else sample_id
                except Exception as e:
                    print(filename, "Does not match with fastq name format: SampleID_Row_Lane_Read_001.fastq.gz")
                else:
                    if sample_id not in sample_files:
                        sample_files[sample_id] = {"R1": [], "R2":[]}
                    if filename.endswith("_1.fq.gz"):
                        sample_files[sample_id]["R1"].append(os.path.join(input_directory, filename))
                    if filename.endswith("_2.fq.gz"):
                        sample_files[sample_id]["R2"].append(os.path.join(input_directory, filename))
    
    find_samples(input_directory)

    if sample_files:
        # Create the output directory for the merged files
        os.makedirs(output_directory, exist_ok=True)
    else:
        print("There is no fastq file in the input directory")
    
    for sample_id, files in sample_files.items():
        r1 = files.get("R1", [])
        r2 = files.get("R2", [])
        r1.sort()
        r2.sort()
        R1, R2 = r1, r2
        if len(R1) > 1:
            if sequence_type == "SE":
                passed = True
            else:
                passed = len(R1) == len(R2)

            if passed:
                try:
                    parts = ' '.join(R1)
                    sample_output_path = os.path.join(output_directory, f'{sample_id}_MF_R1_001.fastq.gz')
                    os.system(f"sudo cat {parts} > {sample_output_path}")
                except Exception as e:
                    print(e)
                else:
                    print(f"{parts} > {sample_output_path}")
                
                if sequence_type == "PE":
                    try:
                        parts2 = ' '.join(R2)
                        sample_output_path2 = os.path.join(output_directory, f"{sample_id}_MF_R2_001.fastq.gz")
                        os.system(f"sudo cat {parts2} > {sample_output_path2}")
                    except Exception as e:
                        print(e)
                    else:
                        print(f"{parts2} > {sample_output_path2}")
            else:
                print(f"Sample {sample_id} has different number of R1 {len(R1)} than R2 ({len(R2)})")
    print("Finished successfully!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ILYOME', description='Ilyome script to merge fastq files')
    parser.add_argument("--version", action="version", version="Ilyome v.3.4")
    parser.add_argument('-dir', metavar="--directory", required=True, type=str, help="The full path to the FASTQ files folder")
    parser.add_argument('-out', metavar="--output", required=True, type=str, help="The output directory")
    parser.add_argument('-st', metavar="--sequence", required=True, type=str, help="Sequencing type : SE for Single end or PE for PE")
    args = parser.parse_args()
    
    if args:
        merge_fastq_files_by_sample_id(args.dir, args.out, args.st)
