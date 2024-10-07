import os

def renamer(directory):
    v = os.listdir(directory)
    for file in v:
        if file.endswith("vcf.gz"):
            path= os.path.join(directory, file)
            newpath=path.replace("g.vcf.gz","hard-filtered.vcf.gz")
            os.system(f"sudo mv {path} {newpath}")  
            print(newpath)

renamer("/mnt/Gennas/02.Bam_VCF__files/26.TWIST-RUN16_PE150_bam_vcf_files/02.vcf_results/01.GATK4_results")
