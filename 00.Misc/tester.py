import os

def listFastqs(directory):
    r1_files = [] 
    r2_files = []
    Fastqs = []

    for i in os.listdir(directory):
        if i.endswith(('_R1_001.fastq.gz', '1.fq.gz', '_1.fastq.gz')):
            r1_files.append(i)
        elif i.endswith(('_R2_001.fastq.gz','2.fq.gz', '_2.fastq.gz')):
            r2_files.append(i)
            
    for r1_file in r1_files:
        if r1_file.endswith("_R1_001.fastq.gz"):
            r2_file = r1_file.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
        elif r1_file.endswith("1.fq.gz"):
            r2_file = r1_file.replace("1.fq.gz", "2.fq.gz")
        elif r1_file.endswith("_1.fastq.gz"):
            r2_file = r1_file.replace("_1.fastq.gz", "_2.fastq.gz")
        
        if r2_file in r2_files:
            r1_realpath = os.path.join(directory, r1_file)
            r2_realpath = os.path.join(directory, r2_file)
            Fastqs.append([r1_file, r2_file, r1_realpath, r2_realpath])
        else:
            print(f"Read 2 is missing for {r1_file}")
            
    for r2_file in r2_files:
        if r2_file.endswith("_R2_001.fastq.gz"):
            r1_file = r2_file.replace("_R2_001.fastq.gz", "_R1_001.fastq.gz")
        elif r2_file.endswith("2.fq.gz"):
            r1_file = r2_file.replace("2.fq.gz", "1.fq.gz")
        elif r2_file.endswith("_2.fastq.gz"):
            r1_file = r2_file.replace("_2.fastq.gz", "_1.fastq.gz")

        if r1_file not in r1_files:
            print(f"Read 1 is missing for {r2_file}")
        return Fastqs


directory = "/mnt/Gennas/01.Fastq_files/169.TWIST-ExoV2-DNAPrepWithExomePlus-RUN159_fastq_files/RUN159-IlluminaDNAPrepwithExome2.0plusEnrichment/TahirAtik"
listFastqs(directory)