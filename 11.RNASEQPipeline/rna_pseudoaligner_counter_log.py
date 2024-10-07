import os
import logging
import subprocess
import pandas as pd
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class PseudoAlignConfig:
    def __init__(self, transcriptome_fasta_file, gtf_file, thread_number, align_tool, bootstrap_value, Rscript_path, data_source, get_ids_path, output_directory, exp_tool, exp_rscript_path, report_rscript_path, template, report_file_name, lfc_cutoff, p_cutoff, interaction):
        self.transcriptome_fasta_file = transcriptome_fasta_file
        self.gtf_file = gtf_file
        self.thread_number = thread_number
        self.align_tool = align_tool
        self.bootstrap_value = bootstrap_value
        self.Rscript_path = Rscript_path
        self.data_source = data_source
        self.get_ids_path = get_ids_path
        self.output_directory = output_directory
        self.exp_tool = exp_tool
        self.exp_rscript_path = exp_rscript_path
        self.report_rscript_path = report_rscript_path
        self.template = template
        self.report_file_name = report_file_name
        self.lfc_cutoff = lfc_cutoff
        self.p_cutoff = p_cutoff
        self.interaction = interaction


class PseudoAlign:
    def __init__(self, transcriptome_fasta_file, output_directory, read1, read2, gtf_file, thread_number, align_tool, bootstrap_value, Rscript_path, data_source, get_ids_path, sample_name):
        self.transcriptome_fasta_file = transcriptome_fasta_file
        self.output_directory = output_directory
        self.read1 = read1
        self.read2 = read2
        self.gtf_file = gtf_file
        self.thread_number = thread_number
        self.align_tool = align_tool
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

        if self.align_tool == "kallisto":
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

        if self.align_tool == "kallisto":
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

        return os.path.join(self.alignment_directory, self.sample_name, "gene_level_counts.txt")



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
        self.expression_dir = os.path.join(output_dir, "pseudoalignment_expression")

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
            
            logging.info(f"rnaseq_report_automator.py - Creating a differential expression analysis report")        
            try:
                result = subprocess.run(f"Rscript {self.report_rscript_path} --template {self.template} --output_directory {self.expression_dir} --output_file {self.report_file_name} --lfc_cutoff {self.lfc_cutoff} --p_cutoff {self.p_cutoff} --interaction {self.interaction} --sample_info_file {self.sample_info_file} --data_source {self.data_source}", shell=True, stderr=subprocess.PIPE, text=True)    
                stderr_output = result.stderr.strip()
                if result.returncode !=0:
                    logging.error(f"edgeR - An error occurred while creating differential expression analysis report: {stderr_output}")
                else:
                    logging.info(f"edgeR - Successfully created differential expression analysis report. Output html file is in: {self.expression_dir}")
            except Exception as e:
                logging.error("rna_pseudoaligner_counter.py - An error occurred while edgeR is creating differential expression analysis report")
        return self.expression_dir        




class ProcessSamples:
    def __init__(self, config, sample_info):
        self.config = config
        self.sample_info = sample_info
        self.index_created = False

    def process_samples(self):
        if not self.index_created:
            pseudoaligner = PseudoAlign(
                transcriptome_fasta_file=self.config.transcriptome_fasta_file,
                gtf_file=self.config.gtf_file,
                output_directory=self.config.output_directory,
                read1=None, 
                read2=None,
                thread_number=self.config.thread_number,
                align_tool=self.config.align_tool,
                bootstrap_value=self.config.bootstrap_value,
                Rscript_path=self.config.Rscript_path,
                data_source=self.config.data_source,
                get_ids_path=self.config.get_ids_path,
                sample_name=None 
            )
            #pseudoaligner.pseudoindexer()
            self.index_created = True 

            for sample in self.sample_info:
                pseudoaligner = PseudoAlign(
                    transcriptome_fasta_file=self.config.transcriptome_fasta_file,
                    gtf_file=self.config.gtf_file,
                    output_directory=self.config.output_directory,
                    read1=sample["read1"],
                    read2=sample["read2"],
                    thread_number=self.config.thread_number,
                    align_tool=self.config.align_tool,
                    bootstrap_value=self.config.bootstrap_value,
                    Rscript_path=self.config.Rscript_path,
                    data_source=self.config.data_source,
                    get_ids_path=self.config.get_ids_path,
                    sample_name=sample["sample_name"]
                )
            
                sample["counts_file"] = pseudoaligner.pseudoaligner()

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
        
        
        sample_df.to_csv(os.path.join(self.config.output_directory,"pseudoalignment_output", "sample_info.txt"), sep="\t", index=False)

            
        differential_expression = DiffExp(
            sample_info_file= os.path.join(self.config.output_directory,"pseudoalignment_output", "sample_info.txt"),
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

common_config = PseudoAlignConfig(
    transcriptome_fasta_file="/home/genwork2/Mert/GRCh38_latest_rna.fna",
    gtf_file="/home/genwork2/Mert/GRCh38_latest_genomic.gtf.gz", 
    thread_number=40,
    align_tool="kallisto",
    bootstrap_value=200,
    Rscript_path="/home/genwork2/Mert/RNAseq/tximport_f_new.R",
    data_source="ncbi",
    get_ids_path="/home/genwork2/Mert/RNAseq/get_ids_c.py",
    output_directory="/home/genwork2/Mert/RNASEQOUTPUT",
    exp_tool= "edgeR",
    exp_rscript_path="/home/genwork2/Mert/RNAseq/edgeR_f_multi.R",
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




