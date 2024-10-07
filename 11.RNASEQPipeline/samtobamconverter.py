import os

def convert_sam_to_bam(input_sam_path, output_bam_path, toolname, thread_number):

    os.system(f"samtools fixmate -O bam,level=0 -@{thread_number} -m {input_sam_path} {os.path.join(output_bam_path, f'{toolname}_fixmate.bam')}")
    
    os.system(f"samtools sort -l 1 -@{thread_number} -o {os.path.join(output_bam_path, f'{toolname}_sorted.bam')} {os.path.join(output_bam_path, f'{toolname}_fixmate.bam')}")

    os.system(f"samtools markdup -r -O bam,level=0 -@{thread_number} {os.path.join(output_bam_path, f'{toolname}_sorted.bam')} {os.path.join(output_bam_path, f'{toolname}_markdup.bam')}")

    os.system(f"samtools view -@{thread_number} {os.path.join(output_bam_path, f'{toolname}_markdup.bam')} -o {os.path.join(output_bam_path, f'{toolname}_final.bam')}")

    os.system(f"samtools index -@{thread_number} {os.path.join(output_bam_path, f'{toolname}_final.bam')}")


inputsam = "/home/genwork2/Mert/NIPT/wbowtie2/24B3042918.sam"

outputbam = "/home/genwork2/Mert/NIPT/wbowtie2"

tool = "bowtie2"

threads = 40

convert_sam_to_bam(input_sam_path=inputsam, output_bam_path=outputbam, toolname=tool, thread_number=threads)