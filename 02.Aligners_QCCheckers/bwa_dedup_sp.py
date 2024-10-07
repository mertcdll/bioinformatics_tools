import os
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class Indexer:
    def __init__(self, reference_fasta:str, algotype:str):
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
    def __init__(self, reference_fasta: str,
                 sample_id: str, 
                 r1:str, 
                 r2:str, 
                 threads:int, 
                 output_path:str, 
                 use_dragen:bool=False,
                 hash_table:str=None):
        
        self.reference_fasta = reference_fasta
        self.sample_id = sample_id
        self.r1 = r1
        self.r2 = r2
        self.threads = threads
        self.output_path = output_path
        self.use_dragen = use_dragen
        self.hash_table = hash_table

    def map(self):

        sample_bam_path = os.path.join(self.output_path, self.sample_id)
        os.makedirs(sample_bam_path, exist_ok=True)
        sam_file = self.sample_id+".sam"
        
        if self.use_dragen:
          logging.info(f"dragen - Mapping {self.r1} and {self.r2} to the reference genome")
          command = (
              f"dragen-os -r {self.hash_table} "
              f"-1 {self.r1} -2 {self.r2} "
              f"--RGSM {self.sample_id} --RGID {self.sample_id} "
              f"--num-threads {self.threads} > {os.path.join(sample_bam_path, sam_file)}"
          )
          
          try:
            os.system(command)
            logging.info(f"dragen - Successfully mapped {self.r1} and {self.r2} to the reference genome")
          except Exception as e:
            logging.error(f"dragen - An error occurred while mapping {self.r1} and {self.r2}: {e}")

        
        else:
          logging.info(f"bwa - Mapping {self.r1} and {self.r2} to the reference genome")          
          command = (
              f"bwa mem -M -t {self.threads} "
              f"-R '@RG\\tID:{self.sample_id}\\tSM:{self.sample_id}\\tLB:Mylib\\tPU:Illumina' "
              f"{self.reference_fasta} {self.r1} {self.r2} > {os.path.join(sample_bam_path, sam_file)}"
          )

          try:
              os.system(command)
              logging.info(f"bwa - Successfully mapped {self.r1} and {self.r2} to the reference genome")
          except Exception as e:
              logging.error(f"bwa - An error occurred while mapping {self.r1} and {self.r2}: {e}")

          return os.path.join(sample_bam_path, sam_file)


class BAMConverter:
    def __init__(self, sample_id, sam_path, threads, use_mark_duplicates, remove_all_duplicates, remove_sequencing_duplicates):
        self.sample_id = sample_id
        self.sam_path = sam_path
        self.threads = threads
        self.use_mark_duplicates = use_mark_duplicates
        self.remove_all_duplicates = remove_all_duplicates
        self.remove_sequencing_duplicates = remove_sequencing_duplicates
        self.out_path = os.path.dirname(self.sam_path)
        self.bam_file = self.sample_id + ".bam"
        self.output_bam = os.path.join(self.out_path, self.bam_file)
        self.filtered_bam_file = "filtered_" + self.sample_id + ".bam"
        self.filtered_bam = os.path.join(self.out_path, self.filtered_bam_file)
        self.metrics_file = self.sample_id + "_dup_metrics.txt"
        self.metrics_path = os.path.join(self.out_path, self.metrics_file)
        self.sorted_bam_file = "sorted_" + self.filtered_bam_file
        self.sorted_bam = os.path.join(self.out_path, self.sorted_bam_file)

    def create_bam_file(self):
        logging.info(f"Samtools - Creating a bam file for {self.sample_id} from sam file")
        converter_cmd = f"samtools view -@ {self.threads} -bS -o {self.output_bam} {self.sam_path}"
        
        try:
            os.system(converter_cmd)
            logging.info(f"Samtools - Successfully created bam file for {self.sample_id}")
        except Exception as e:
            logging.error(f"Samtools - An error occurred while creating bam file for {self.sample_id}: {e}")
    
    def mark_duplicates(self):
        if self.use_mark_duplicates:
            logging.info(f"MarkDuplicatesSpark - Marking and removing duplicates for bam file of {self.sample_id}")
            marker_cmd = f"gatk MarkDuplicatesSpark -I {self.output_bam} -O {self.filtered_bam} --metrics-file {self.metrics_path} --spark-master local[{self.threads}]"
            
            if self.remove_sequencing_duplicates:
                marker_cmd += " --remove-sequencing-duplicates"
            
            if self.remove_all_duplicates:
                marker_cmd += " --remove-all-duplicates"
                                
            try:
                os.system(marker_cmd)
                logging.info(f"MarkDuplicatesSpark - Successfully marked and removed duplicates for {self.sample_id} > {self.filtered_bam}")
            except Exception as e:
                logging.error(f"An error occurred while marking and removing duplicates for bam file of {self.sample_id}: {e}")
        else:
            logging.info(f"Sambamba - Marking and removing duplicates for bam file of {self.sample_id}")
            marker_cmd = f"sambamba markdup -t {self.threads} {self.output_bam} {self.filtered_bam}"
            
            if self.remove_all_duplicates:
                marker_cmd += " --remove-duplicates"
            
            try:
                os.system(marker_cmd)
                logging.info(f"Sambamba - Successfully marked and removed duplicates for {self.sample_id} > {self.filtered_bam}")
            except Exception as e:
                logging.error(f"An error occurred while marking and removing duplicates for bam file of {self.sample_id}: {e}")
    
    def sort_bam(self):
        sort_cmd = f"samtools sort -@ {self.threads} {self.filtered_bam} -o {self.sorted_bam}"
        
        try:
            os.system(sort_cmd)
            logging.info(f"Samtools - Successfully sorted bam file for {self.sample_id} > {self.sorted_bam}")
        except Exception as e:
            logging.error(f"Samtools - An error occurred while sorting bam file for {self.sample_id}: {e}")
    
    def index_bam(self):
        index_cmd = f"samtools index -@ {self.threads} {self.sorted_bam}"
        
        try:
            os.system(index_cmd)
            logging.info(f"Samtools - Successfully indexed bam file for {self.sample_id} > {self.sorted_bam}")
        except Exception as e:
            logging.error(f"Samtools - An error occurred while indexing bam file for {self.sample_id}: {e}")

    def clean_up(self):
        try:
            os.remove(self.sam_path)
            os.remove(self.output_bam)
            os.remove(self.filtered_bam)
            logging.info(f"Clean up - Successfully removed intermediate files for {self.sample_id}")
        except Exception as e:
            logging.error(f"Clean up - An error occurred while removing intermediate files for {self.sample_id}: {e}")

    def convert_markdedup(self):
        self.create_bam_file()
        self.mark_duplicates()
        self.sort_bam()
        self.index_bam()
        self.clean_up()


ref_fasta = "/home/genwork2/Mert/01.REFGEN/GRCh37/hg19.fa"
at = "bwtsw"

out_path = "/home/genwork2/Mert/NIPT"
r1 = "/home/genwork2/03.Fastq_Output/240902_NDX550927_0503_AHGJT7BDXY/24B3043751_S2_R1_001.fastq.gz"
r2 = "/home/genwork2/03.Fastq_Output/240902_NDX550927_0503_AHGJT7BDXY/24B3043751_S2_R2_001.fastq.gz"
sample_id = "24B3043751"
threads = 40
remove_dups = True
remove_seq_dups = False
use_md = False

#indexit = Indexer(reference_fasta=ref_fasta, algotype=at) # use once
#indexit.index() # use once

mapit = Mapper(reference_fasta=ref_fasta, sample_id=sample_id, r1=r1, r2=r2, output_path=out_path, threads=threads)
sam_file= mapit.map()

dedup = BAMConverter(sam_path=sam_file, sample_id=sample_id, threads=threads, remove_all_duplicates=remove_dups, remove_sequencing_duplicates=remove_seq_dups, use_mark_duplicates=use_md)
dedup.convert_markdedup()


