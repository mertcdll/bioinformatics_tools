import os
import logging
import subprocess

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')



class PseudoAlign:
    def __init__(self, transcriptome_fasta_file, output_directory, read1, read2, gtf_file, thread_number, tool, bootstrap_value, Rscript_path, data_source, get_ids_path, sample_name):
        self.transcriptome_fasta_file = transcriptome_fasta_file
        self.output_directory = output_directory
        self.read1 = read1
        self.read2 = read2
        self.gtf_file = gtf_file
        self.thread_number = thread_number
        self.tool = tool
        self.bootstrap_value = bootstrap_value
        self.Rscript_path = Rscript_path
        self.data_source = data_source
        self.get_ids_path = get_ids_path
        self.alignment_directory = os.path.join(self.output_directory, "pseudoalignment_output")
        self.index_directory = os.path.join(self.output_directory, "pseudoalignment_index")
        self.sample_name = sample_name

    def pseudoindexer(self):
        if self.transcriptome_fasta_file.endswith(".gz"):
            logging.info("rna_pseudoaligner_counter.py - Unzipping transcriptome fasta file")
            try:
                result = subprocess.run(f"gunzip {self.transcriptome_fasta_file}", shell=True, stderr=subprocess.PIPE, text=True)
                stderr_output = result.stderr.strip()
                if result.returncode != 0:
                    logging.error(f"Gunzip - An error occurred while unzipping the transcriptome fasta file: {stderr_output} ")
                else:
                    logging.info("Gunzip - Transcriptome fasta file is successfully unzipped")
                decompressed_fasta_file = self.transcriptome_fasta_file[:-3]
            except Exception as e: 
                logging.error(f"rna_pseudoaligner_counter.py - An internal error occurred while unzipping the transcriptome fasta file: {e}")
        else:
            decompressed_fasta_file = self.transcriptome_fasta_file

        os.makedirs(self.index_directory, exist_ok=True)
        os.makedirs(self.alignment_directory, exist_ok=True)

        if self.tool == "kallisto":
            logging.info("Kallisto - Indexing the reference transcriptome")
            try:
                result = subprocess.run(f"kallisto index -i {os.path.join(self.index_directory, 'transcriptome_ind.idx')} {decompressed_fasta_file}", shell=True, stderr=subprocess.PIPE, text=True)
                stderr_output = result.stderr.strip()
                if result.returncode != 0:
                    logging.error(f"Kallisto - An error occurred while indexing the reference transcriptome:{stderr_output}")
                else:
                    logging.info(f"Kallisto - Created the index file succesfully in the directory: {self.index_directory}")
            except Exception as e:
                logging.error(f"rna_pseudoaligner_counter.py - An error occurred while Kallisto is indexing the reference transcriptome:{e}")
            logging.info("get_ids.py - Getting Gene IDs of transcripts from transcriptome fasta file")
            try:
                result = subprocess.run(f"python3 {self.get_ids_path} --fasta_file {decompressed_fasta_file} --data_source {self.data_source} --output_directory {self.alignment_directory}", shell=True, stderr=subprocess.PIPE, text=True)
                stderr_output = result.stderr.strip()
                if result.returncode != 0:
                    logging.error(f"get_ids.py - An error occurred while getting Gene IDs of transcripts: {stderr_output}")
                else: 
                    logging.info(f"get_ids.py - Transcript IDs vs Gene IDs were written in the file: {os.path.join(self.alignment_directory, 'tx2gene.txt')}")
            except Exception as e:
                logging.error(f"rna_pseudoaligner_counter.py - An error occurred while get_ids.py is getting Gene IDs of transcripts: {e}")    
    def pseudoaligner(self):
        # Paired-end
        if self.read2:
            if self.read1.endswith('.gz'):
                logging.info("Gunzip - Unzipping the fastq file for read 1")
                try:
                    result = subprocess.run(f"gunzip {self.read1}", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode != 0:
                        logging.error(f"Gunzip - An error occurred while unzipping the fastq file for read 1: {stderr_output}")
                    else: 
                        logging.info("Gunzip - Fastq file of read 1 is successfully unzipped")
                except Exception as e:
                    logging.error(f"rna_pseudoaligner_counter.py - An error occurred while unzipping the fastq file for read 1: {e}")
                decompressedfile1 = self.read1[:-3]
            else:
                decompressedfile1 = self.read1

            if self.read2.endswith('.gz'):
                logging.info("Gunzip - Unzipping the fastq file for read 2")
                try:
                    result = subprocess.run(f"gunzip {self.read2}", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode != 0:
                        logging.error(f"Gunzip - An error occurred while unzipping the fastq file for read 2: {stderr_output}")
                    else:
                        logging.info("Gunzip - Fastq file of read 2 is successfully unzipped") 
                except Exception as e:
                    logging.error(f"rna_pseudoaligner_counter.py - An error occurred while unzipping the fastq file for read 2: {e}")
                decompressedfile2 = self.read2[:-3]
            else:
                decompressedfile2 = self.read2
        # Single-end
        else:
            if self.read1.endswith('.gz'):
                logging.info("Gunzip - Unzipping the fastq file")
                try:
                    result = subprocess.run(f"gunzip {self.read1}", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode != 0:
                        logging.error(f"Gunzip - An error occurred while unzipping the fastq file for read 2: {stderr_output}")
                    else:
                        logging.info("Gunzip - Fastq file is successfully unzipped")
                except Exception as e:
                    logging.error(f"rna_pseudoaligner_counter.py - An error occurred while unzipping the fastq file: {e}")
                decompressedfile1 = self.read1[:-3]
            else:
                decompressedfile1 = self.read1

        if self.tool == "kallisto":
            if self.read2:
                logging.info(f"Kallisto - Getting transcript level abundances for {self.sample_name}")
                try:
                    result = subprocess.run(f"kallisto quant -i {os.path.join(self.index_directory, 'transcriptome_ind.idx')} -o {os.path.join(self.alignment_directory, self.sample_name)} {decompressedfile1} {decompressedfile2} --threads {self.thread_number} --bootstrap-samples {self.bootstrap_value} --plaintext --gtf {self.gtf_file}", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode != 0:
                        logging.error(f"Kallisto - An error occurred while kallisto is getting transcript level abundances for {self.sample_name}: {stderr_output}")
                    else: 
                        logging.info(f"Kallisto  - Successfully calculated transcript level abundances for {self.sample_name} in the file: {os.path.join(self.alignment_directory, self.sample_name, 'abundance.tsv')}")
                except Exception as e:
                    logging.error(f"rna_pseudoaligner_counter.py - An error occurred while Kallisto is getting transcript level abundances for {self.sample_name}: {e}")
                logging.info(f"Tximport - Calculating gene level read counts from transcript level abundances for {self.sample_name}")
                try:
                    result = subprocess.run(f"Rscript {self.Rscript_path} --id_file {os.path.join(self.alignment_directory, 'tx2gene.txt')} --abundance_file {os.path.join(self.alignment_directory, self.sample_name, 'abundance.tsv')} --output_dir {os.path.join(self.alignment_directory, self.sample_name)} --sample_name {self.sample_name}", shell=True, stderr=subprocess.PIPE, text=True) 
                    stderr_output = result.stderr.strip()
                    if result.returncode != 0:
                        logging.error(f"Tximport - An error occurred while calculating gene level counts from transcript level abundances for {self.sample_name}: {stderr_output}")
                    else:
                        logging.info(f"Tximport - Succesfully calculated gene level read counts, gene level abundances and gene lengths for {self.sample_name} in the directory: {os.path.join(self.alignment_directory, self.sample_name)}")    
                except Exception as e:
                    logging.error(f"rna_pseudoaligner_counter.py - An error occurred while Tximport is calculating gene level read counts from transcript level abundances for {self.sample_name}: {e}")
            else:
                logging.info(f"Kallisto - Getting transcript level abundances for {self.sample_name}")
                try:
                    result = subprocess.run(f"kallisto quant -i {os.path.join(self.index_directory, 'transcriptome_ind.idx')} -o {os.path.join(self.alignment_directory, self.sample_name)} {decompressedfile1} --threads {self.thread_number} --bootstrap-samples {self.bootstrap_value} --plaintext --gtf {self.gtf_file} --single", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode != 0:
                        logging.error(f"Kallisto - An error occurred while getting transcript level abundances for {self.sample_name}: {stderr_output}")
                    else: 
                        logging.info(f"Kallisto - Successfully calculated transcript level abundances for {self.sample_name} in the file: {os.path.join(self.alignment_directory, self.sample_name, 'abundance.tsv')}")
                except Exception as e:
                    logging.error(f"rna_pseudoaligner_counter.py - An error occurred while Kallisto is getting transcript level abundances for {self.sample_name}: {e}")                
                logging.info(f"Tximport - Calculating gene level readcounts from transcript level abundances for {self.sample_name}")
                try:
                    result = subprocess.run(f"Rscript {self.Rscript_path} --id_file {os.path.join(self.alignment_directory, 'tx2gene.txt')} --abundance_file {os.path.join(self.alignment_directory, self.sample_name, 'abundance.tsv')} --output_dir {os.path.join(self.alignment_directory, self.sample_name)} --sample_name {self.sample_name}", shell=True, stderr=subprocess.PIPE, text=True) 
                    stderr_output = result.stderr.strip()
                    if result.returncode != 0:
                        logging.error(f"Tximport - An error occurred while calculating gene level read counts from transcript level abundances for {self.sample_name}: {stderr_output}")
                    else:
                        logging.info(f"Tximport - Succesfully calculated gene level read counts, gene level abundances and gene lengths for {self.sample_name} in the directory: {os.path.join(self.alignment_directory, self.sample_name)}")    
                except Exception as e:
                    logging.error(f"rna_pseudoaligner_counter.py - An error occurred while Tximport is calculating gene level read counts from transcript level abundances for {self.sample_name}: {e}")

pseudoalignutil = "kallisto"
transcriptome_fasta_path = "/home/genwork2/Mert/GRCh38_latest_rna.fna" # change according to data source
gtf_path = "/home/genwork2/Mert/GRCh38_latest_genomic.gtf.gz" # change according to data source
output_path = "/home/genwork2/Mert/RNAseq"
thread_num = 40
bootstrap_val = 200
read_1 = "/home/genwork2/Mert/RNAseq/samplernaseqdata/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq"
read_2 = "/home/genwork2/Mert/RNAseq/samplernaseqdata/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq"
rscript_dir = "/home/genwork2/Mert/RNAseq/tximport_f_new.R"
data_origin = "ncbi" # change according to data source
id_ext_path = "/home/genwork2/Mert/RNAseq/get_ids_c.py"
sample_id = "UHR_Rep3_ERCC-Mix1"


# If run type is single read, type read2=None
pseudoaligner = PseudoAlign(
    transcriptome_fasta_file=transcriptome_fasta_path,
    output_directory=output_path,
    read1=read_1,
    read2=read_2,
    gtf_file=gtf_path,
    thread_number=thread_num,
    tool=pseudoalignutil,
    bootstrap_value=bootstrap_val,
    Rscript_path=rscript_dir,
    data_source=data_origin,
    get_ids_path=id_ext_path,
    sample_name=sample_id
)

#pseudoaligner.pseudoindexer() # use it once
pseudoaligner.pseudoaligner()