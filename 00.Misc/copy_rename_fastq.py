import os
import subprocess

samples = [
  {
    "sample_id": "gfan4697-RE-RUN171-Twist",
    "r1": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/gfan4697-RE_S40_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/gfan4697-RE_S40_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "gfan4737-RE-RUN161-Twist",
    "r1": "/mnt/gen40/01.fastq_files/171.TWIST-ExoV2-DNAPrepWithExomePlus-RUN161_fastq_files/RUN161-TwistExoV2/Arman/gfan4737-RE_S72_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen40/01.fastq_files/171.TWIST-ExoV2-DNAPrepWithExomePlus-RUN161_fastq_files/RUN161-TwistExoV2/Arman/gfan4737-RE_S72_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "gfan4737-RE-RUN171-Twist",
    "r1": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/gfan4737-RE_S43_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/gfan4737-RE_S43_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "gfan4699-RE-RUN171-Twist",
    "r1": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/gfan4699-RE_S48_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/gfan4699-RE_S48_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "gfan4699-RE-RUN169-Twist",
    "r1": "/mnt/gen40/01.fastq_files/179.TWIST-ExoV2-DNAPrepWithExomePlus-RUN169_fastq_files/RUN169-Twist/Arman/AR-61-E-793_MF_R1_001.fastq.gz",
    "r2": "/mnt/gen40/01.fastq_files/179.TWIST-ExoV2-DNAPrepWithExomePlus-RUN169_fastq_files/RUN169-Twist/Arman/AR-61-E-793_MF_R2_001.fastq.gz"
  },
  {
    "sample_id": "AR-61-E-793-RE2-RUN171-Twist",
    "r1": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/AR-61-E-793-RE2_S39_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/AR-61-E-793-RE2_S39_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "AR-61-E-807-RUN169-Twist",
    "r1": "/mnt/gen40/01.fastq_files/179.TWIST-ExoV2-DNAPrepWithExomePlus-RUN169_fastq_files/RUN169-Twist/Arman/AR-61-E-807_S194_L003_R1_001.fastq.gz",
    "r2": "/mnt/gen40/01.fastq_files/179.TWIST-ExoV2-DNAPrepWithExomePlus-RUN169_fastq_files/RUN169-Twist/Arman/AR-61-E-807_S194_L003_R2_001.fastq.gz"
  },
  {
    "sample_id": "AR-61-E-807-RE-RUN169-Twist",
    "r1": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/AR-61-E-807-RE_S24_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/AR-61-E-807-RE_S24_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "AR-61-E-739-RE-RUN171-Twist",
    "r1": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/AR-61-E-739-RE_S47_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/AR-61-E-739-RE_S47_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "Dx48731-RE-RUN169-Twist",
    "r1": "/mnt/gen40/01.fastq_files/179.TWIST-ExoV2-DNAPrepWithExomePlus-RUN169_fastq_files/RUN169-Twist/Najmabadi/Dx48731-RE_S61_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen40/01.fastq_files/179.TWIST-ExoV2-DNAPrepWithExomePlus-RUN169_fastq_files/RUN169-Twist/Najmabadi/Dx48731-RE_S61_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "Dx48731-RE2-RUN171-Twist",
    "r1": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Najmabadi/Dx48731-RE2_S26_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Najmabadi/Dx48731-RE2_S26_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "Dx48273-RE-RUN169-Twist",
    "r1": "/mnt/gen40/01.fastq_files/179.TWIST-ExoV2-DNAPrepWithExomePlus-RUN169_fastq_files/RUN169-Twist/Najmabadi/Dx48273-RE_S68_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen40/01.fastq_files/179.TWIST-ExoV2-DNAPrepWithExomePlus-RUN169_fastq_files/RUN169-Twist/Najmabadi/Dx48273-RE_S68_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "Dx48731-RE2-RUN171-Twist",
    "r1": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Najmabadi/Dx48273-RE2_S28_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Najmabadi/Dx48273-RE2_S28_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "Dx48141-RE-RUN169-Twist",
    "r1": "/mnt/gen40/01.fastq_files/179.TWIST-ExoV2-DNAPrepWithExomePlus-RUN169_fastq_files/RUN169-Twist/Najmabadi/Dx48141-RE_S65_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen40/01.fastq_files/179.TWIST-ExoV2-DNAPrepWithExomePlus-RUN169_fastq_files/RUN169-Twist/Najmabadi/Dx48141-RE_S65_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "Dx48141-RE2-RUN171-Twist",
    "r1": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Najmabadi/Dx48141-RE2_S53_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Najmabadi/Dx48141-RE2_S53_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "AR-Mil86-747-RE-RUN171-Twist",
    "r1": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/AR-Mil86-747-RE_S70_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/Arman/AR-Mil86-747-RE_S70_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "299-24-MAG-RUN171-Twist",
    "r1": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/GenKli/299-24-MAG_S38_L001_R1_001.fastq.gz",
    "r2": "/mnt/gen100/01.fastq_files/181.TWIST-ExoV2-DNAPrepWithExomePlus-RUN171_fastq_files/RUN171-Twist/GenKli/299-24-MAG_S38_L001_R2_001.fastq.gz"
  },
  {
    "sample_id": "1871-23-RKS-RE-RUN164-Twist",
    "r1": "/mnt/gen40/01.fastq_files/174.TWIST-ExoV2-DNAPrepWithExomePlus-RUN164_fastq_files/RUN164-Twist/GenKli/1871-23-RKS-RE_S41_R1_001.fastq.gz",
    "r2": "/mnt/gen40/01.fastq_files/174.TWIST-ExoV2-DNAPrepWithExomePlus-RUN164_fastq_files/RUN164-Twist/GenKli/1871-23-RKS-RE_S41_R2_001.fastq.gz"
  }
]

output_directory = "/mnt/gen100/01.fastq_files/194-ILLUMINA-DEMO-RE-NovaSeq-RUN184_fastq_files/TestTwist"

os.makedirs(output_directory, exist_ok=True)

def copy_and_rename_files(sample):
    run_id = '-'.join(sample["sample_id"].split("-")[-2:])
    r1_source = sample["r1"]
    r2_source = sample["r2"]
    
      
    sample1_name = os.path.basename(r1_source).split("_")[0]
    sample1_remainder =  "_".join(os.path.basename(r1_source).split("_")[1:])
    r1_dest_name = sample1_name + "-" + run_id + "_" + sample1_remainder
    
    
    
    sample2_name = os.path.basename(r2_source).split("_")[0]
    sample2_remainder =  "_".join(os.path.basename(r2_source).split("_")[1:])
    r2_dest_name = sample2_name + "-" + run_id + "_" + sample2_remainder
    
    
    r1_dest = os.path.join(output_directory, r1_dest_name)
    r2_dest = os.path.join(output_directory, r2_dest_name)

    subprocess.run(["sudo", "cp", r1_source, r1_dest], check=True)
    subprocess.run(["sudo", "cp", r2_source, r2_dest], check=True)

    print(f"Copied and renamed: {r1_source} -> {r1_dest}")
    print(f"Copied and renamed: {r2_source} -> {r2_dest}")

for sample in samples:
    copy_and_rename_files(sample)

print("Finished copying and renaming files.")


