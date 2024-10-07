import os
import logging
import csv
import json
import argparse


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
        
        return self.sorted_bam
        
        
class QCMetrics:

    def __init__(self, bam_path, target, bait, ref_fasta, out_path, sample_id):
        self.bam_path = bam_path
        self.target = target
        self.bait = bait
        self.ref_fasta = ref_fasta
        self.out_path = out_path
        self.sample_id = sample_id
        self.qc_path = None
        self.csv_path = None

    def qc(self):
        
        logging.info("gatk CollectHsMetrics - Collecting QC metrics from bam file")
        
        self.qc_path = os.path.join(self.out_path, f"{self.sample_id}_hs_metrics.txt")
        
        command = (
            f"gatk CollectHsMetrics "
            f"--INPUT {self.bam_path} "
            f"--OUTPUT {self.qc_path} "
            f"--REFERENCE_SEQUENCE {self.ref_fasta} "
            f"--BAIT_INTERVALS {self.bait} "
            f"--TARGET_INTERVALS {self.target}"
        )

        logging.info(f"Running command: {command}")
        os.system(command)


    def convert_metrics_to_csv(self):
        with open(self.qc_path, 'r') as metrics_file:
            lines = metrics_file.readlines()
        
        start_idx = None
        for idx, line in enumerate(lines):
            if line.startswith('## METRICS CLASS'):
                start_idx = idx + 1
                break
        
        if start_idx is None:
            raise ValueError("Metrics class not found in the file.")
        
        headers = lines[start_idx].strip().split('\t')
        data = lines[start_idx + 1].strip().split('\t')
        
        self.csv_path = os.path.join(self.out_path, f"{self.sample_id}_hs_metrics.csv")
        
        with open(self.csv_path, 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(headers)
            writer.writerow(data)

def process_samples(config_path):
    with open(config_path, 'r') as config_file:
        config = json.load(config_file)

    reference_fasta = config["reference_fasta"]
    algotype = config["algotype"]
    targets = config["targets"]
    baits = config["baits"]
    output_path = config["output_path"]
    threads = config["aligner_threads"]
    remove_all_duplicates = config["remove_all_duplicates"]
    remove_sequencing_duplicates = config["remove_sequencing_duplicates"]
    use_mark_duplicates = config["use_gatk_mark_duplicates"]
    use_dragen = config["use_dragen"]
    samples = config["samples"]
    index_fasta = config["index_fasta"]
    
    if index_fasta:
    
        indexer = Indexer(reference_fasta=reference_fasta, algotype=algotype)
        indexer.index()

    for sample in samples:
        sample_id = sample["sample_id"]
        r1 = sample["r1"]
        r2 = sample["r2"]

        mapper = Mapper(reference_fasta=reference_fasta, sample_id=sample_id, r1=r1, r2=r2, output_path=output_path, threads=threads, use_dragen=use_dragen)
        sam_file = mapper.map()

        deduplicator = BAMConverter(sam_path=sam_file, sample_id=sample_id, threads=threads, remove_all_duplicates=remove_all_duplicates, remove_sequencing_duplicates=remove_sequencing_duplicates, use_mark_duplicates=use_mark_duplicates)
        bam_file = deduplicator.convert_markdedup()

        qc_metrics = QCMetrics(bam_path=bam_file, target=targets, bait=baits, ref_fasta=reference_fasta, out_path=output_path, sample_id=sample_id)
        qc_metrics.qc()
        qc_metrics.convert_metrics_to_csv()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process genomic data.")
    parser.add_argument("--config", type=str, required=True, help="Path to the configuration file.")
    args = parser.parse_args()

    process_samples(args.config)