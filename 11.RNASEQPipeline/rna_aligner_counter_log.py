import os
import logging
import subprocess
import pandas as pd
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class AlignConfig:
    def __init__(self, fasta_file, gtf_file, thread_number, align_tool, count_tool, output_directory, picard_directory, data_source, refflat_script, rrna_intlist_script, gtftogenepred, seq_type, exp_tool, exp_rscript_path, report_rscript_path, template, report_file_name, lfc_cutoff, p_cutoff, interaction):
        self.fasta_file = fasta_file
        self.gtf_file = gtf_file
        self.thread_number = thread_number
        self.align_tool = align_tool
        self.count_tool = count_tool
        self.output_directory = output_directory
        self.picard_directory = picard_directory
        self.data_source = data_source
        self.refflat_script = refflat_script
        self.rrna_intlist_script = rrna_intlist_script
        self.gtftogenepred = gtftogenepred
        self.seq_type = seq_type
        self.exp_tool = exp_tool
        self.exp_rscript_path = exp_rscript_path
        self.report_rscript_path = report_rscript_path
        self.template = template
        self.report_file_name = report_file_name
        self.lfc_cutoff = lfc_cutoff
        self.p_cutoff = p_cutoff
        self.interaction = interaction



class Indexer:
    def __init__(self, fasta_file, gtf_file, output_directory, thread_number, tool):
        self.fasta_file = fasta_file
        self.output_directory = output_directory
        self.thread_number = thread_number
        self.gtf_file = gtf_file 
        self.tool = tool
        self.hisat2_build_path = "/home/genwork2/hisat2/hisat2-build" # Unfortunately hisat2 indexer needs this.

    def indexing(self):
        if self.fasta_file.endswith(".gz"):
            logging.info("rna_aligner_counter.py - Unzipping genome fasta file")
            try:
                result = subprocess.run(f"gunzip {self.fasta_file}", shell=True, stderr=subprocess.PIPE, text=True)
                stderr_output = result.stderr.strip()
                if result.returncode !=0:
                    logging.error(f"Gunzip - An error occurred while unzipping the genome fasta file: {stderr_output} ")
                else:
                     logging.info("Gunzip - Genome fasta file is successfully unzipped")
                   
                decompressed_fasta_file = self.fasta_file[:-3]
            except Exception as e:
                logging.error(f"rna_aligner_counter.py - An internal error occurred while unzipping the genome fasta file: {e}")

        else:
            decompressed_fasta_file = self.fasta_file

        index_directory = os.path.join(self.output_directory, "alignment_index")
        os.makedirs(index_directory, exist_ok=True)


        if self.tool == 'star':
            logging.info("STAR - Indexing the reference genome")
            try:
                if self.gtf_file == None:
                    result = subprocess.run(f"STAR --runThreadN {self.thread_number} --runMode genomeGenerate --genomeDir {index_directory} --genomeFastaFiles {decompressed_fasta_file}", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode !=0:
                        logging.error(f"STAR - An error occurred while indexing the reference genome:{stderr_output}")
                    else:
                        logging.info(f"STAR - Created the index file succesfully in the directory: {index_directory}")
                else:
                    result = subprocess.run(f"STAR --runThreadN {self.thread_number} --runMode genomeGenerate --genomeDir {index_directory} --genomeFastaFiles {decompressed_fasta_file} --sjdbGTFfile {self.gtf_file}", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode !=0:
                        logging.error(f"STAR - An error occurred while indexing the reference genome:{stderr_output}")
                    else:
                        logging.info(f"STAR - Created the index file succesfully in the directory: {index_directory}")
            except Exception as e:
                logging.error(f"rna_aligner_counter.py - An error occurred while Star is indexing the reference genome:{e}")
        elif self.tool == 'hisat2':
            logging.info("HISAT2 - Indexing the reference genome")
            try:
                result = subprocess.run(f"python3 {self.hisat2_build_path} --threads {self.thread_number} {decompressed_fasta_file} {os.path.join(index_directory, 'genome_ind')}", shell=True, stderr=subprocess.PIPE, text=True)
                stderr_output = result.stderr.strip()
                if result.returncode !=0:
                    logging.error(f"HISAT2 - An error occurred while indexing the reference genome:{stderr_output}")
                else:
                    logging.info(f"HISAT2 - Created the index file succesfully in the directory: {index_directory}")
            except Exception as e:
                logging.error(f"rna_aligner_counter.py - An error occurred while Hisat2 is indexing the reference genome:{e}")
        elif self.tool == 'magicblast':
            logging.info("Magicblast - Indexing the reference genome")
            try:
                result = subprocess.run(f"makeblastdb -in {decompressed_fasta_file} -out {os.path.join(index_directory, 'genome_ind')} -parse_seqids -dbtype nucl", shell=True, stderr=subprocess.PIPE, text=True )
                stderr_output = result.stderr.strip()
                if result.returncode !=0:
                    logging.error(f"Magicblast - An error occurred while indexing the reference genome:{stderr_output}")
                else:
                    logging.info(f"Magicblast - Created the index file succesfully in the directory: {index_directory}")
            except Exception as e:
                logging.error(f"rna_aligner_counter.py - An error occurred while Magicblast is indexing the reference genome:{e}")


        return index_directory


class Aligner:
    def __init__(self, sample_name, read1, read2, index_directory, output_directory, thread_number, tool):
        self.sample_name = sample_name
        self.read1 = read1
        self.read2 = read2
        self.index_directory = index_directory
        self.output_directory = output_directory
        self.thread_number = thread_number
        self.tool = tool

    def align(self):
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
                    logging.error(f"rna_aligner_counter.py - An error occurred while unzipping the fastq file for read 1: {e}")
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
                    logging.error(f"rna_aligner_counter.py - An error occurred while unzipping the fastq file for read 2: {e}")
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
                    logging.error(f"rna_aligner_counter.py - An error occurred while unzipping the fastq file: {e}")
                decompressedfile1 = self.read1[:-3]
            else:
                decompressedfile1 = self.read1

        alignment_directory = os.path.join(self.output_directory, "alignment_output")
        os.makedirs(alignment_directory, exist_ok=True)

        sample_directory = os.path.join(alignment_directory, self.sample_name)
        os.makedirs(sample_directory, exist_ok=True)

        if self.tool == 'star':
            if self.read2:
                logging.info(f"STAR - Aligning {self.sample_name} to reference genome")
                try:
                    result = subprocess.run(f"STAR --runThreadN {self.thread_number} --genomeDir {self.index_directory} --readFilesIn {decompressedfile1} {decompressedfile2} --outFileNamePrefix {os.path.join(sample_directory, '')}", shell=True, stderr=subprocess.PIPE, text=True)    
                    stderr_output = result.stderr.strip()
                    if result.returncode !=0:
                        logging.error(f"STAR - An error occurred while aligning {self.sample_name} to reference genome:{stderr_output}")
                    else: 
                        logging.info(f"STAR - Successfully aligned {self.sample_name} to reference genome in the file: {os.path.join(sample_directory, 'Aligned.out.sam')}")
                except Exception as e:
                    logging.error(f"rna_aligner_counter.py - An error occurred while star is aligning {self.sample_name} to reference genome: {e}")
            else:
                logging.info(f"STAR - Aligning {self.sample_name} to reference genome")
                try:
                    result = subprocess.run(f"STAR --runThreadN {self.thread_number} --genomeDir {self.index_directory} --readFilesIn {decompressedfile1} --outFileNamePrefix {os.path.join(sample_directory, '')}", shell=True, stderr=subprocess.PIPE, text=True)    
                    stderr_output = result.stderr.strip()
                    if result.returncode !=0:
                        logging.error(f"STAR - An error occurred while aligning {self.sample_name} to reference genome:{stderr_output}")
                    else: 
                        logging.info(f"STAR - Successfully aligned {self.sample_name} to reference genome in the file: {os.path.join(sample_directory, 'Aligned.out.sam')}")
                except Exception as e:
                    logging.error(f"rna_aligner_counter.py - An error occurred while star is aligning {self.sample_name} to reference genome: {e}")     
        elif self.tool == 'hisat2':
            if self.read2:
                logging.info(f"HISAT2 - Aligning {self.sample_name} to reference genome")
                try:
                    result = subprocess.run(f"hisat2 -p {self.thread_number} -x {os.path.join(self.index_directory, 'genome_ind')} -1 {decompressedfile1} -2 {decompressedfile2} -S {os.path.join(sample_directory, 'Aligned.out.sam')} > {os.path.join(sample_directory, 'initial_metrics.txt')}", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode !=0:
                        logging.error(f"HISAT2 - An error occurred while aligning {self.sample_name} to reference genome:{stderr_output}")
                    else: 
                        logging.info(f"HISAT2 - Successfully aligned {self.sample_name} to reference genome in the file: {os.path.join(sample_directory, 'Aligned.out.sam')}")
                except Exception as e:
                    logging.error(f"rna_aligner_counter.py - An error occurred while hisat2 is aligning {self.sample_name} to reference genome: {e}")
            else:
                logging.info(f"HISAT2 - Aligning {self.sample_name} to reference genome")
                try:
                    result = subprocess.run(f"hisat2 -p {self.thread_number} -x {os.path.join(self.index_directory, 'genome_ind')} -U {decompressedfile1} -S {os.path.join(sample_directory, 'Aligned.out.sam')} > {os.path.join(sample_directory, 'initial_metrics.txt')}", shell=True, stderr=subprocess.PIPE, text=True)       
                    stderr_output = result.stderr.strip()
                    if result.returncode !=0:
                        logging.error(f"HISAT2 - An error occurred while aligning {self.sample_name} to reference genome:{stderr_output}")
                    else: 
                        logging.info(f"HISAT2 - Successfully aligned {self.sample_name} to reference genome in the file: {os.path.join(sample_directory, 'Aligned.out.sam')}")
                except Exception as e:
                    logging.error(f"rna_aligner_counter.py - An error occurred while hisat2 is aligning {self.sample_name} to reference genome: {e}")
        elif self.tool == 'magicblast':
            if self.read2:
                logging.info(f"Magicblast - Aligning {self.sample_name} to reference genome")
                try:
                    result = subprocess.run(f"magicblast -query {decompressedfile1} -query_mate {decompressedfile2} -db {os.path.join(self.index_directory, 'genome_ind')} -num_threads {self.thread_number} -out {os.path.join(sample_directory, 'Aligned.out.sam')} -infmt fastq -outfmt sam", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode !=0:
                        logging.error(f"Magicblast - An error occurred while aligning {self.sample_name} to reference genome:{stderr_output}")
                    else: 
                        logging.info(f"HISAT2 - Successfully aligned {self.sample_name} to reference genome in the file: {os.path.join(sample_directory, 'Aligned.out.sam')}")
                except Exception as e:
                    logging.error(f"rna_aligner_counter.py - An error occurred while magicblast is aligning {self.sample_name} to reference genome: {e}")                    
            else:
                logging.info(f"Magicblast - Aligning {self.sample_name} to reference genome")
                try:
                    result = subprocess.run(f"magicblast -query {decompressedfile1} -db {os.path.join(self.index_directory, 'genome_ind')} -num_threads {self.thread_number} -out {os.path.join(sample_directory, 'Aligned.out.sam')} -infmt fastq -outfmt sam", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode !=0:
                        logging.error(f"Magicblast - An error occurred while aligning {self.sample_name} to reference genome:{stderr_output}")
                    else: 
                        logging.info(f"HISAT2 - Successfully aligned {self.sample_name} to reference genome in the file: {os.path.join(sample_directory, 'Aligned.out.sam')}")
                except Exception as e:
                    logging.error(f"rna_aligner_counter.py - An error occurred while magicblast is aligning {self.sample_name} to reference genome: {e}")  

        return sample_directory
        
class SamtoBamConversion:
    def __init__(self, conversion_directory, thread_number):
        self.conversion_directory = conversion_directory
        self.thread_number = thread_number
        self.converter()
        
    def converter(self):
        logging.info("Samtools - Started to process sam file to create a bam file")
        logging.info("Samtools fixmate - Correcting missing mate information by adjusting the alignment coordinates and flags of the reads")
        try:
            result = subprocess.run(f"samtools fixmate -O bam,level=0 -@{self.thread_number} -m {os.path.join(self.conversion_directory, 'Aligned.out.sam')} {os.path.join(self.conversion_directory, 'fixmate.bam')}", shell=True, stderr=subprocess.PIPE, text=True)
            stderr_output = result.stderr.strip()
            if result.returncode !=0:
                logging.error(f"Samtools fixmate - An error occurred while performing fixmate module: {stderr_output}")
            else:
                logging.info("Samtools fixmate - Successfully corrected missing mate information")
        except Exception as e:
            logging.error(f"rna_aligner_counter.py - An error occurred while samtools fixmate correcting missing mate information: {e}")
        
        logging.info("Samtools sort - Sorting the alignment in records based on their alignment positions")
        try:
            result = subprocess.run(f"samtools sort -l 1 -@{self.thread_number} -o {os.path.join(self.conversion_directory, 'sorted.bam')} {os.path.join(self.conversion_directory, 'fixmate.bam')}", shell=True, stderr=subprocess.PIPE, text=True)
            stderr_output = result.stderr.strip()
            if result.returncode !=0:
                logging.error(f"Samtools sort - An error occurred while performing sort module: {stderr_output}")
            else:
                logging.info("Samtools sort - Successfully corrected missing mate information")
        except Exception as e:
            logging.error(f"rna_aligner_counter.py - An error occurred while samtools sorting alignment records: {e}")                    
        
        logging.info("Samtools markdup - Marking duplicates")
        try:
            result = subprocess.run(f"samtools markdup -O bam,level=0 -@{self.thread_number} {os.path.join(self.conversion_directory, 'sorted.bam')} {os.path.join(self.conversion_directory, 'markdup.bam')}", shell=True, stderr=subprocess.PIPE, text=True)
            stderr_output = result.stderr.strip()
            if result.returncode !=0:
                logging.error(f"Samtools markdup - An error occurred while performing markdup module: {stderr_output}")
            else:
                logging.info("Samtools markdup - Successfully marked duplicates")
        except Exception as e:
            logging.error(f"rna_aligner_counter.py - An error occurred while samtools markdup marking duplicates: {e}")  
        
        logging.info("Samtools view - Creating final bam file")
        try:
            result = subprocess.run(f"samtools view -@{self.thread_number} {os.path.join(self.conversion_directory, 'markdup.bam')} -o {os.path.join(self.conversion_directory, 'final.bam')}", shell=True, stderr=subprocess.PIPE, text=True)
            stderr_output = result.stderr.strip()
            if result.returncode !=0:
                logging.error(f"Samtools view - An error occurred while performing view module: {stderr_output}")
            else:
                logging.info(f"Samtools view - Successfully created final bam file in: {os.path.join(self.conversion_directory, 'final.bam')}")
        except Exception as e:
            logging.error(f"rna_aligner_counter.py - An error occurred while samtools view creating final bam file: {e}")          
        
        logging.info("Samtools index - Creating index file")
        try: 
            result = subprocess.run(f"samtools index -@{self.thread_number} {os.path.join(self.conversion_directory, 'final.bam')}", shell=True, stderr=subprocess.PIPE, text=True)
            stderr_output = result.stderr.strip()
            if result.returncode !=0:
                logging.error(f"Samtools index - An error occurred while performing index module: {stderr_output}")
            else:
                logging.info(f"Samtools index - Successfully created final index file in: {self.conversion_directory}")
        except Exception as e:
            logging.error(f"rna_aligner_counter.py - An error occurred while samtools index creating final index file: {e}")                  
        
        os.remove(os.path.join(self.conversion_directory, 'fixmate.bam'))
        os.remove(os.path.join(self.conversion_directory, 'sorted.bam'))
        os.remove(os.path.join(self.conversion_directory, 'markdup.bam'))

class RefflatRNAint:
   
    def __init__(self, gtf_file, output_directory, gtftogenepred_path, data_source, ref_flat_sc_path, rrna_int_list_sc_path):
        self.gtf_file = gtf_file
        self.output_directory = output_directory
        self.gtftogenepred_path = gtftogenepred_path
        self.data_source = data_source
        self.ref_flat_sc_path = ref_flat_sc_path
        self.rrna_int_list_sc_path = rrna_int_list_sc_path
        self.createfiles()

    def createfiles(self):
        logging.info("rna_aligner_counter.py - Creating refflat and rRNA interval list")
        logging.info("refflat_c.py - Creating refflat file")
        try:
            result = subprocess.run(f"python3 {self.ref_flat_sc_path} --gtf_path {self.gtf_file} --output_path {self.output_directory} --data_source {self.data_source} --program_path {self.gtftogenepred_path}", shell=True, stderr=subprocess.PIPE, text=True)
            stderr_output = result.stderr.strip()
            if result.returncode !=0:
                logging.error(f"refflat_c.py - An error occurred while creating refflat file: {stderr_output}")
            else:
                logging.info(f"refflat_c.py - Successfully created refflat file in: {self.output_directory}")
        except Exception as e:
            logging.error(f"rna_aligner_counter.py - An error occurred while refflat_c.py creating refflat file: {e}")   
        
        logging.info("rrna_int_c.py - Creating rRNA interval list")
        try:
            result = subprocess.run(f"python3 {self.rrna_int_list_sc_path} --gtf_path {self.gtf_file} --bam_directory {self.output_directory} --output_directory {self.output_directory}", shell=True, stderr=subprocess.PIPE, text=True)
            stderr_output = result.stderr.strip()
            if result.returncode !=0:
                logging.error(f"rrna_int_c.py - An error occurred while creating rRNA interval list file: {stderr_output}")
            else:
                logging.info(f"rrna_int_c.py - Successfully created rRNA interval list file in: {self.output_directory}")
        except Exception as e:
            logging.error(f"rna_aligner_counter.py - An error occurred while rrna_int_c.py creating rRNA interval list file: {e}")   

class QCMetrics:
    def __init__(self, qc_directory, picard_path, thread_number):
        self.qc_directory = qc_directory
        self.picard_path = picard_path
        self.thread_number = thread_number
        self.qc()
    
    def qc(self): 
        logging.info("rna_aligner_counter.py - Creating qc metrics")
        logging.info("Samtools flagstat - Calculating qc metrics")
        try:
            result = subprocess.run(f"samtools flagstat -@{self.thread_number} {os.path.join(self.qc_directory, 'final.bam')} > {os.path.join(self.qc_directory, 'flagstat_output.txt')}", shell = True, stderr = subprocess.PIPE, text = True) 
            stderr_output = result.stderr.strip()
            if result.returncode !=0:
                logging.error(f"Samtools flagstat - An error occurred while calculating qc metrics: {stderr_output}")
            else:
                logging.info(f"Samtools flagstat - Successfully calculated qc metrics in the directory: {self.qc_directory}")
        except Exception as e:
            logging.error(f"rna_aligner_counter.py - An error occurred while Samtools flagstat calculating qc metrics: {e}")

        logging.info("Picard CollectRnaSeqMetrics - Collecting RNASeq metrics")
        try:
            result = subprocess.run(f"java -jar {self.picard_path} CollectRnaSeqMetrics -I {os.path.join(self.qc_directory, 'final.bam')} -O {os.path.join(self.qc_directory, 'output.collect_RNA_Metrics')} -REF_FLAT {os.path.join(self.qc_directory, 'refflat.gp')} -STRAND SECOND_READ_TRANSCRIPTION_STRAND -RIBOSOMAL_INTERVALS {os.path.join(self.qc_directory, 'rrna_intlist.txt')} -CHART {os.path.join(self.qc_directory, 'PositionvsCoverage.pdf')}", shell=True, stderr=subprocess.PIPE, text=True)
            stderr_output = result.stderr.strip()
            if result.returncode !=0:
                logging.error(f"Picard CollectRnaSeqMetrics - An error occurred while collecting RNASeq metrics: {stderr_output}")
            else:
                logging.info(f"Picard CollectRnaSeqMetrics - Successfully collected RNASeq metrics in the directory: {self.qc_directory}")
        except Exception as e:
            logging.error(f"rna_aligner_counter.py - An error occurred while Picard CollectRnaSeqMetrics collecting RNASeq metrics: {e}")

class Quant:
    def __init__(self, sample_name, output_directory, alignment_directory, gtf_file, thread_number, tool, read_type):
        self.sample_name = sample_name
        self.output_directory = output_directory
        self.gtf_file = gtf_file
        self.tool = tool
        self.alignment_directory = alignment_directory
        self.thread_number = thread_number
        self.read_type = read_type


    def quantification(self):

        read_count_dir = os.path.join(self.output_directory, "read_counts")
        os.makedirs(read_count_dir, exist_ok=True)

        sample_directory = os.path.join(read_count_dir, self.sample_name)
        os.makedirs(sample_directory, exist_ok=True)

        if self.tool == "featurecounts":
            logging.info("Featurecounts - Calculating gene level counts")
            if self.read_type == "paired-end":
                try: 
                    result = subprocess.run(f"featureCounts -T {self.thread_number} -t exon -g gene_id -a {self.gtf_file} -o {os.path.join(sample_directory, 'featurecounts_counts.txt')} {os.path.join(self.alignment_directory, 'final.bam')} -p --countReadPairs", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode !=0:
                        logging.error(f"Featurecounts - An error occurred while calculating gene level counts: {stderr_output}")
                    else:
                        logging.info("Featurecounts - Successfully calculated gene level counts")
                except Exception as e:
                    logging.error(f"rna_aligner_counter.py - An error occurred while Featurecounts is calculating gene level counts: {e}")
            elif self.read_type == "single-end": 
                try:
                    result = subprocess.run(f"featureCounts -T {self.thread_number} -t exon -g gene_id -a {self.gtf_file} -o {os.path.join(sample_directory, 'featurecounts_counts.txt')} {os.path.join(self.alignment_directory, 'final.bam')}", shell=True, stderr=subprocess.PIPE, text=True)               
                    stderr_output = result.stderr.strip()
                    if result.returncode !=0:
                        logging.error(f"Featurecounts - An error occurred while calculating gene level counts: {stderr_output}")
                    else:
                        logging.info("Featurecounts - Successfully calculated gene level counts")
                except Exception as e:
                    logging.error(f"rna_aligner_counter.py - An error occurred while Featurecounts is calculating gene level counts: {e}")
            else:
                raise ValueError("Invalid read_type. Supported vales are 'paired-end' and 'single-end'")
            
            counts = pd.read_csv(os.path.join(sample_directory, 'featurecounts_counts.txt'), sep="\t", comment="#")

            counts = counts.iloc[:, [0,-1]]

            with open(os.path.join(sample_directory, 'gene_level_counts.txt'), "w") as file:
                file.write(self.sample_name + "\n")

                counts.to_csv(file, sep="\t", header=False, index=False)

        if self.tool == "htseq":
            logging.info("HTSeq - Calculating gene level counts")
            if self.read_type == "paired-end":
                try:
                    result = subprocess.run(f"htseq-count -f bam -r pos -i gene_id -m union -s no {os.path.join(self.alignment_directory, 'final.bam')} {self.gtf_file} > {os.path.join(sample_directory, 'htseq_counts.txt')}", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode !=0:
                        logging.error(f"HTSeq - An error occurred while calculating gene level counts: {stderr_output}")
                    else:
                        logging.info("HTSeq - Successfully calculated gene level counts")
                except Exception as e:
                    logging.error(f"rna_aligner_counter.py - An error occurred while HTSeq is calculating gene level counts: {e}")
            elif self.read_type == "single-end":
                try:
                    result = subprocess.run(f"htseq-count -f bam -r pos -i gene_id -m union -s no {os.path.join(self.alignment_directory, 'final.bam')} {self.gtf_file} > {os.path.join(sample_directory, 'htseq_counts.txt')}", shell=True, stderr=subprocess.PIPE, text=True)
                    stderr_output = result.stderr.strip()
                    if result.returncode !=0:
                        logging.error(f"HTSeq - An error occurred while calculating gene level counts: {stderr_output}")
                    else:
                        logging.info("HTSeq - Successfully calculated gene level counts")
                except Exception as e:
                    logging.error(f"rna_aligner_counter.py - An error occurred while HTSeq is calculating gene level counts: {e}")                    

            else:
                raise ValueError("Invalid read_type. Supported values are 'paired-end' and 'single-end'")
            
            sample_name = f"{self.sample_name}\n"

            with open(os.path.join(sample_directory, 'htseq_counts.txt'), "r") as file:
                lines = file.readlines()
        
            lines = [line for line in lines if not line.startswith("__")]

            lines.insert(0, sample_name)

            with open(os.path.join(sample_directory, 'gene_level_counts.txt'), "w") as file:
                file.writelines(lines)


        return os.path.join(sample_directory, 'gene_level_counts.txt')
        
class DiffExp:
    def __init__ (self, sample_info_file, output_dir, exp_tool, exp_rscript_path, data_source, report_rscript_path, template, report_file_name, lfc_cutoff, p_cutoff, interaction):
        self.sample_info_file = sample_info_file
        self.output_dir = output_dir
        self.exp_tool = exp_tool
        self.exp_rscript_path = exp_rscript_path
        self.data_source = data_source
        self.report_rscript_path = report_rscript_path
        self.template = template
        self.report_file_name = report_file_name
        self.lfc_cutoff = lfc_cutoff
        self.p_cutoff = p_cutoff
        self.interaction = interaction
        self.expression_dir = os.path.join(output_dir, "alignment_expression")

    def differential_exp(self):
        
        os.makedirs(self.expression_dir, exist_ok=True)

        if self.exp_tool == "edgeR":
            logging.info(f"edgeR - Performing differential expression analysis")        
            try:
                result = subprocess.run(f"Rscript {self.exp_rscript_path} --sample_info_file {self.sample_info_file} --output_dir {self.expression_dir} --data_source {self.data_source} --interaction {self.interaction} --lfc_cutoff {self.lfc_cutoff} --p_cutoff {self.p_cutoff}", shell=True, stderr=subprocess.PIPE, text=True)    
                stderr_output = result.stderr.strip()
                if result.returncode !=0:
                    logging.error(f"edgeR - An error occurred while performing differential expression analysis: {stderr_output}")
                else:
                    logging.info(f"edgeR - Successfully performed differential expression analysis. Output files are in: {self.expression_dir}")
            except Exception as e:
                logging.error("rna_pseudoaligner_counter.py - An error occurred while edgeR is performing differential expression analysis")
            
            logging.info(f"rnaseq_report_automator.R - Creating a differential expression analysis report")        
            try:
                result = subprocess.run(f"Rscript {self.report_rscript_path} --template {self.template} --output_directory {self.expression_dir} --output_file {self.report_file_name} --lfc_cutoff {self.lfc_cutoff} --p_cutoff {self.p_cutoff} --interaction {self.interaction} --sample_info_file {self.sample_info_file} --data_source {self.data_source}", shell=True, stderr=subprocess.PIPE, text=True)    
                stderr_output = result.stderr.strip()
                if result.returncode !=0:
                    logging.error(f"edgeR - An error occurred while performing differential expression analysis: {stderr_output}")
                else:
                    logging.info(f"edgeR - Successfully performed differential expression analysis. Output files are in: {self.expression_dir}")
            except Exception as e:
                logging.error("rna_pseudoaligner_counter.py - An error occurred while edgeR is performing differential expression analysis")
        return self.expression_dir  


class ProcessSamples:
    def __init__(self, config, sample_info):
        self.config = config
        self.sample_info = sample_info
        self.index_created = False
    
    def process_samples(self):
        if not self.index_created:
            indexer = Indexer(
                fasta_file=self.config.fasta_file,
                gtf_file=None,
                output_directory=self.config.output_directory,
                thread_number=self.config.thread_number,
                tool=self.config.align_tool
            )
            index_dir = indexer.indexing()
            self.index_created = True

            for sample in self.sample_info:
                aligner = Aligner(
                    sample_name=sample["sample_name"],
                    read1=sample["read1"],
                    read2=sample["read2"],
                    index_directory=index_dir,
                    output_directory=self.config.output_directory,
                    thread_number=self.config.thread_number,
                    tool=self.config.align_tool
                )

                align_dir = aligner.align()

                converter = SamtoBamConversion(
                    conversion_directory=align_dir,
                    thread_number=self.config.thread_number
                )

                refflat_rnaint = RefflatRNAint(
                    gtf_file=self.config.gtf_file,
                    output_directory=align_dir,
                    gtftogenepred_path=self.config.gtftogenepred,
                    data_source=self.config.data_source,
                    ref_flat_sc_path=self.config.refflat_script,
                    rrna_int_list_sc_path=self.config.rrna_intlist_script,                
                )

                qc_controller = QCMetrics(
                    qc_directory=align_dir,
                    picard_path=self.config.picard_directory,
                    thread_number=self.config.thread_number
                )


                counts = Quant(
                    sample_name=sample["sample_name"],
                    output_directory=self.config.output_directory,
                    gtf_file=self.config.gtf_file,
                    alignment_directory=align_dir,
                    thread_number=self.config.thread_number,
                    tool=self.config.count_tool,
                    read_type=self.config.seq_type

                )

                sample["counts_file"] = counts.quantification()

        names = [item["sample_name"] for item in self.sample_info]
        subjects = [item["subject"] for item in self.sample_info]
        count_paths = [item["counts_file"] for item in self.sample_info]

        met_df = pd.DataFrame({"names" : names, "subject": subjects, "count_dirs": count_paths})
        
        pattern = r'^f\d+'
        factors = []
        values = []

        for item in self.sample_info:
            for key in item:
                if re.match(pattern, key):
                    factors.append(key)
                    values.append(item[key])

        un_factors = factors[0:int(len(factors)/len(self.sample_info))]

        num_factors = len(un_factors)

        factor_dict = {factor: [] for factor in un_factors}

        for i in range(0, len(values), num_factors):
            for j, factor in enumerate(un_factors):
                factor_dict[factor].append(values[i + j])

        factor_df = pd.DataFrame(factor_dict)

        fac_names = [item.split("_")[1] for item in un_factors]

        factor_df.columns = fac_names
        
        sample_df = pd.concat([met_df, factor_df], axis=1)
        
        
        sample_df.to_csv(os.path.join(self.config.output_directory, "read_counts", "sample_info.txt"), sep="\t", index=False)


        differential_expression = DiffExp(
            sample_info_file= os.path.join(self.config.output_directory,"read_counts", "sample_info.txt"),
            output_dir=self.config.output_directory,
            exp_tool=self.config.exp_tool,
            exp_rscript_path=self.config.exp_rscript_path,
            data_source=self.config.data_source,
            report_rscript_path=self.config.report_rscript_path,
            template=self.config.template,
            report_file_name=self.config.report_file_name,
            lfc_cutoff=self.config.lfc_cutoff,
            p_cutoff=self.config.p_cutoff,
            interaction=self.config.interaction
        )

        differential_expression.differential_exp()


common_config = AlignConfig(
    fasta_file= "/home/genwork2/Mert/GRCh38_latest_genomic.fna.gz",
    gtf_file= "/home/genwork2/Mert/GRCh38_latest_genomic.gtf.gz",
    thread_number= 40,
    align_tool= "magicblast",
    count_tool= "featurecounts",
    output_directory= "/home/genwork2/Mert/RNAseq/ALIGNEROUTPUTFE",
    picard_directory= "/home/genwork2/Mert/apps/picard.jar",
    data_source= "ncbi",
    refflat_script= "/home/genwork2/Mert/gtftorefflat/refflat_c.py",
    rrna_intlist_script= "/home/genwork2/Mert/gtftorefflat/rrna_int_c.py",
    gtftogenepred= "/home/genwork2/Mert/gtftorefflat/gtfToGenePred",
    seq_type= "paired-end",
    exp_tool="edgeR",
    exp_rscript_path= "/home/genwork2/Mert/RNAseq/edgeR_f_multi.R",
    report_rscript_path="/home/genwork2/Mert/RNAseq/rnaseq_report_automator.R",
    template="/home/genwork2/Mert/RNAseq/DifferentialExpWithedgeRTemplate.Rmd",
    report_file_name="DEReport.html",
    lfc_cutoff="1.2",
    p_cutoff="0.05",
    interaction="True"
)



sample_info = [
    {
        "sample_name": "HBR_Rep1",
        "subject" : "s1",
        "read1": "/home/genwork2/Mert/RNAseq/samplernaseqdata/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq",
        "read2": "/home/genwork2/Mert/RNAseq/samplernaseqdata/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq",
        "f1_treatment": "HBR",
        "f2_drug":"A",
        "organism": "HomoSapiens"
    },
    {
        "sample_name": "HBR_Rep2",
        "subject":"s1",
        "read1": "/home/genwork2/Mert/RNAseq/samplernaseqdata/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq",
        "read2": "/home/genwork2/Mert/RNAseq/samplernaseqdata/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq", 
        "f1_treatment": "HBR",
        "f2_drug": "A",
        "organism" : "HomoSapiens" 
    },
    {
        "sample_name": "HBR_Rep3",
        "subject": "s1",
        "read1": "/home/genwork2/Mert/RNAseq/samplernaseqdata/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq",
        "read2": "/home/genwork2/Mert/RNAseq/samplernaseqdata/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq",
        "f1_treatment": "HBR",
        "f2_drug": "B",
        "organism": "HomoSapiens"  
    },
    {
        "sample_name": "UHR_Rep1",
        "subject": "s2",
        "read1": "/home/genwork2/Mert/RNAseq/samplernaseqdata/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq",
        "read2": "/home/genwork2/Mert/RNAseq/samplernaseqdata/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq",
        "f1_treatment": "UHR",
        "f2_drug": "A",
        "organism": "HomoSapiens"  
    },
    {
        "sample_name": "UHR_Rep2",
        "subject": "s2",
        "read1": "/home/genwork2/Mert/RNAseq/samplernaseqdata/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq",
        "read2": "/home/genwork2/Mert/RNAseq/samplernaseqdata/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq",
        "f1_treatment": "UHR",
        "f2_drug": "A",
        "organism": "HomoSapiens"  
    },
    {
        "sample_name": "UHR_Rep3",
        "subject": "s2",
        "read1": "/home/genwork2/Mert/RNAseq/samplernaseqdata/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq",
        "read2": "/home/genwork2/Mert/RNAseq/samplernaseqdata/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq",
        "f1_treatment": "UHR",
        "f2_drug": "B",
        "organism": "HomoSapiens"  
    }
]


processor = ProcessSamples(common_config, sample_info)

processor.process_samples()
