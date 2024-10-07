import os
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class Indexer:
    def __init__(self, reference_fasta, algotype):
        self.reference_fasta = reference_fasta
        self.algotype = algotype


    def index(self):
        
        logging.info("bwa - Indexing the reference genome")
        command  = f"bwa index -a {self.algotype} {self.reference_fasta}"

        try:
            os.system(command)
            logging.info("bwa - Successfully indexed the reference genome")
        except Exception as e: 
            logging.error(f"bwa - An error occurred while indexing the reference genome: {e}")


class Mapper:
    def __init__(self, reference_fasta, sample_id, r1, r2, threads, output_path):
        self.reference_fasta = reference_fasta
        self.sample_id = sample_id
        self.r1 = r1
        self.r2 = r2
        self.threads = threads
        self.output_path = output_path

    def map(self):

        sample_bam_path = os.path.join(self.output_path, self.sample_id)
        os.makedirs(sample_bam_path, exist_ok=True)

        logging.info(f"bwa - Mapping {self.r1} and {self.r2} to the reference genome")
        sam_file = self.sample_id+".sam"
        command = f"bwa mem -M -t {self.threads} -R '@RG\\tID:{self.sample_id}\\tSM:{self.sample_id}\\tLB:Mylib\\tPU:Illumina' {self.reference_fasta} {self.r1} {self.r2} > {os.path.join(sample_bam_path, sam_file)}"


        try:
            os.system(command)
            logging.info(f"bwa - Successfully mapped {self.r1} and {self.r2} to the reference genome")
        except Exception as e:
            logging.error(f"bwa - An error occurred while mapping {self.r1} and {self.r2}: {e}")

        return os.path.join(sample_bam_path, sam_file)



class deduplicate:
    def __init__(self, sam_path, sample_id, threads):
        self.sam_path = sam_path
        self.threads = threads
        self.sample_id = sample_id

    def convert_markdedup(self):

        logging.info(f"Sambamba - Creating a bam file for {self.sample_id} from sam file")
        bam_file = self.sample_id+".bam"
        
        out_path = os.path.dirname(self.sam_path)

        output_bam = os.path.join(out_path, bam_file)
        converter_cmd = f"sambamba view -t {self.threads} -S -f bam -o {output_bam} {self.sam_path}"

        try: 
            os.system(converter_cmd)
            logging.info(f"Sambamba - Successfully created bam file for {self.sample_id}")
        except Exception as e:
            logging.error(f"Sambamba - An error occurred while creating bam file for {self.sample_id}")

        logging.info(f"Sambamba - Marking and removing duplicates for bam file of {self.sample_id}")
        
        filtered_bam_file = "filtered_" + self.sample_id + ".bam"
        filtered_bam = os.path.join(out_path, filtered_bam_file)

        marker_cmd = f"sambamba markdup -t {threads} -r {output_bam} {filtered_bam}"

        try:
            os.system(marker_cmd)
            logging.info(f"Sambamba - Successfully marked and removed duplicates for {self.sample_id} > {filtered_bam}")
        except Exception as e:
            logging.error(f"An error occurred while marking and removing duplicates for bam file of {self.sample_id}")

        os.system(f"rm {self.sam_path}")
        os.system(f"rm {output_bam}")

        sorted_bam_file = "sorted_" + filtered_bam_file
        sorted_bam = os.path.join(out_path, sorted_bam_file)

        sort_cmd = f"samtools sort -@ {self.threads} {filtered_bam} -o {sorted_bam}"

        os.system(sort_cmd)

        index_cmd = f"samtools index -@ {self.threads} {sorted_bam}"

        os.system(index_cmd)

        os.system(f"rm {filtered_bam}")


ref_fasta = "/home/genwork2/Mert/01.REFGEN/GRCh38/hg38.fa"
at = "bwtsw"

out_path = "/home/genwork2/Mert/bwa_dedup/brcadual"
r1 = "/home/genwork2/03.Fastq_Output/RUN181-BRCA-DUAL/BRCA-DUAL/771-24-HOE_S1_L004_R1_001.fastq.gz"
r2 = "/home/genwork2/03.Fastq_Output/RUN181-BRCA-DUAL/BRCA-DUAL/771-24-HOE_S1_L004_R2_001.fastq.gz"
sample_id = "771-24-HOE"
threads = 40

#indexit = Indexer(reference_fasta=ref_fasta, algotype=at) # use once
#indexit.index() # use once

mapit = Mapper(reference_fasta=ref_fasta, sample_id=sample_id, r1=r1, r2=r2, output_path=out_path, threads=threads)
sam_file= mapit.map()

dedup = deduplicate(sam_path=sam_file, sample_id=sample_id, threads=threads)
dedup.convert_markdedup()


