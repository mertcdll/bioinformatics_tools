import os

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
            os.system(f"gunzip {self.fasta_file}") #with open?
            decompressed_fasta_file = self.fasta_file[:-3]
        else:
            decompressed_fasta_file = self.fasta_file

        index_directory = os.path.join(self.output_directory, "genome_index")
        os.makedirs(index_directory, exist_ok=True)


        if self.tool == 'star':
            if self.gtf_file == None:
                os.system(f"STAR --runThreadN {self.thread_number} --runMode genomeGenerate --genomeDir {index_directory} --genomeFastaFiles {decompressed_fasta_file}")
            else:
                os.system(f"STAR --runThreadN {self.thread_number} --runMode genomeGenerate --genomeDir {index_directory} --genomeFastaFiles {decompressed_fasta_file} --sjdbGTFfile {self.gtf_file}")
        elif self.tool == 'hisat2':
            os.system(f"python3 {self.hisat2_build_path} --threads {self.thread_number} {decompressed_fasta_file} {os.path.join(index_directory, 'genome_ind')}")
        elif self.tool == 'magicblast':
            os.system(f"makeblastdb -in {decompressed_fasta_file} -out {os.path.join(index_directory, 'genome_ind')} -parse_seqids -dbtype nucl")
        
        return index_directory


class Aligner:
    def __init__(self, read1, read2, index_directory, output_directory, thread_number, tool):
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
                os.system(f"gunzip {self.read1}") 
                decompressedfile1 = self.read1[:-3]
            else:
                decompressedfile1 = self.read1

            if self.read2.endswith('.gz'):
                os.system(f"gunzip {self.read2}")
                decompressedfile2 = self.read2[:-3]
            else:
                decompressedfile2 = self.read2
        # Single-end
        else:

            if self.read1.endswith('.gz'):
                os.system(f"gunzip {self.read1}") 
                decompressedfile1 = self.read1[:-3]
            else:
                decompressedfile1 = self.read1

        alignment_directory = os.path.join(self.output_directory, "alignment_output")
        os.makedirs(alignment_directory, exist_ok=True)

        if self.tool == 'star':
            if self.read2:
                os.system(f"STAR --runThreadN {self.thread_number} --genomeDir {self.index_directory} --readFilesIn {decompressedfile1} {decompressedfile2} --outFileNamePrefix {os.path.join(alignment_directory, '')}")    
            else:
                os.system(f"STAR --runThreadN {self.thread_number} --genomeDir {self.index_directory} --readFilesIn {decompressedfile1} --outFileNamePrefix {os.path.join(alignment_directory, '')}")    
        elif self.tool == 'hisat2':
            if self.read2:
                os.system(f"hisat2 -p {self.thread_number} -x {os.path.join(self.index_directory, 'genome_ind')} -1 {decompressedfile1} -2 {decompressedfile2} -S {os.path.join(alignment_directory, 'Aligned.out.sam')} > {os.path.join(alignment_directory, 'initial_metrics.txt')}")
            else:
                os.system(f"hisat2 -p {self.thread_number} -x {os.path.join(self.index_directory, 'genome_ind')} -U {decompressedfile1} -S {os.path.join(alignment_directory, 'Aligned.out.sam')} > {os.path.join(alignment_directory, 'initial_metrics.txt')}")       
        elif self.tool == 'magicblast':
            if self.read2:
                os.system(f"magicblast -query {decompressedfile1} -query_mate {decompressedfile2} -db {os.path.join(self.index_directory, 'genome_ind')} -num_threads {self.thread_number} -out {os.path.join(alignment_directory, 'Aligned.out.sam')} -infmt fastq -outfmt sam")
            else:
                os.system(f"magicblast -query {decompressedfile1} -db {os.path.join(self.index_directory, 'genome_ind')} -num_threads {self.thread_number} -out {os.path.join(alignment_directory, 'Aligned.out.sam')} -infmt fastq -outfmt sam")
        
        return alignment_directory
        
class SamtoBamConversion:
    def __init__(self, conversion_directory, thread_number):
        self.conversion_directory = conversion_directory
        self.thread_number = thread_number
        self.converter()
        
    def converter(self):
        os.system(f"samtools fixmate -O bam,level=0 -@{self.thread_number} -m {os.path.join(self.conversion_directory, 'Aligned.out.sam')} {os.path.join(self.conversion_directory, 'fixmate.bam')}")
        os.system(f"samtools sort -l 1 -@{self.thread_number} -o {os.path.join(self.conversion_directory, 'sorted.bam')} {os.path.join(self.conversion_directory, 'fixmate.bam')}")
        os.system(f"samtools markdup -O bam,level=0 -@{self.thread_number} {os.path.join(self.conversion_directory, 'sorted.bam')} {os.path.join(self.conversion_directory, 'markdup.bam')}")
        os.system(f"samtools view -@{self.thread_number} {os.path.join(self.conversion_directory, 'markdup.bam')} -o {os.path.join(self.conversion_directory, 'final.bam')}")
        os.system(f"samtools index -@{self.thread_number} {os.path.join(self.conversion_directory, 'final.bam')}")
                
        
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
        os.system(f"python3 {self.ref_flat_sc_path} --gtf_path {self.gtf_file} --output_path {self.output_directory} --data_source {self.data_source} --program_path {self.gtftogenepred_path}")
        os.system(f"python3 {self.rrna_int_list_sc_path} --gtf_path {self.gtf_file} --bam_directory {self.output_directory} --output_directory {self.output_directory}")


class QCMetrics:
    def __init__(self, qc_directory, picard_path, thread_number):
        self.qc_directory = qc_directory
        self.picard_path = picard_path
        self.thread_number = thread_number
        self.qc_control()
    
    def qc_control(self): 
        os.system(f"samtools flagstat -@{self.thread_number} {os.path.join(self.qc_directory, 'final.bam')} > {os.path.join(self.qc_directory, 'flagstat_output.txt')}") 
        os.system(f"java -jar {self.picard_path} CollectRnaSeqMetrics -I {os.path.join(self.qc_directory, 'final.bam')} -O {os.path.join(self.qc_directory, 'output.collect_RNA_Metrics')} -REF_FLAT {os.path.join(self.qc_directory, 'refflat.gp')} -STRAND SECOND_READ_TRANSCRIPTION_STRAND -RIBOSOMAL_INTERVALS {os.path.join(self.qc_directory, 'rrna_intlist.txt')}")


class Quant:
    def __init__(self, output_directory, alignment_directory, gtf_file, thread_number, tool, read_type):
        self.output_directory = output_directory
        self.gtf_file = gtf_file
        self.tool = tool
        self.alignment_directory = alignment_directory
        self.thread_number = thread_number
        self.read_type = read_type
        self.quantification()


    def quantification(self):

        read_count_dir = os.path.join(self.output_directory, "read_counts")
        os.makedirs(read_count_dir, exist_ok=True)

        if self.tool == "featurecounts":
            if self.read_type == "paired-end":
                os.system(f"featureCounts -T {self.thread_number} -t exon -g gene_id -a {self.gtf_file} -o {os.path.join(read_count_dir, 'gene_level_counts.txt')} {os.path.join(self.alignment_directory, 'final.bam')} -p --countReadPairs")
            elif self.read_type == "single-end": 
                os.system(f"featureCounts -T {self.thread_number} -t exon -g gene_id -a {self.gtf_file} -o {os.path.join(read_count_dir, 'gene_level_counts.txt')} {os.path.join(self.alignment_directory, 'final.bam')}")               
            else:
                raise ValueError("Invalid read_type. Supported vales are 'paired-end' and 'single-end'")
        if self.tool == "htseq":
            if self.read_type == "paired-end":
                os.system(f"htseq-count -f bam -r pos -i gene_id -m union -s no {os.path.join(self.alignment_directory, 'final.bam')} {self.gtf_file} > {os.path.join(read_count_dir, 'gene_level_counts.txt')}")
            elif self.read_type == "single-end":
                os.system(f"htseq-count -f bam -r pos -i gene_id -m union -s no {os.path.join(self.alignment_directory, 'final.bam')} {self.gtf_file} > {os.path.join(read_count_dir, 'gene_level_counts.txt')}")
            else:
                raise ValueError("Invalid read_type. Supported values are 'paired-end' and 'single-end'")

        

alignutil = "hisat2" # change if you want to use a different aligner (other options are star and magicblast)
countutil = "featurecounts" # change if you want to use a different counter (other option is htseq)
genome_fasta_path = "/home/genwork2/Mert/RNAseq/Homo_sapiens.GRCh38.dna.toplevel.fa" # change according to data source
transcriptome_fasta_path =  "/home/genwork2/Mert/RNAseq/human_transcriptome.fa" #change according to data source
gtf_path = "/home/genwork2/Mert/RNAseq/Homo_sapiens.GRCh38.110.gtf.gz" #change according to data source
output_path = "/home/genwork2/Mert/RNAseq"
thread_num = 40
read_1 = "/home/genwork2/Mert/RNAseq/67-ZEYNEP_AKKUS_TH_S59_L003_R1_001.fastq"
read_2 = "/home/genwork2/Mert/RNAseq/67-ZEYNEP_AKKUS_TH_S59_L003_R2_001.fastq"
picard_dir = "/home/genwork2/Mert/apps/picard.jar"
data_origin = "ensembl" #change according to data source
refflat_script_path = "/home/genwork2/Mert/gtftorefflat/refflat_c.py"
rrna_intlist_script_path = "/home/genwork2/Mert/gtftorefflat/rrna_int_c.py"
gtftogenepred_pr_path = "/home/genwork2/Mert/gtftorefflat/gtfToGenePred"
seqtype = "paired-end" # change according to read type

align_dir = "/home/genwork2/Mert/RNAseq/alignment_output"

# Indexer can be used once
indexer = Indexer(fasta_file=genome_fasta_path, gtf_file=None, output_directory=output_path, thread_number=thread_num, tool=alignutil)
index_dir = indexer.indexing()

# Then align fastq files (If run type is single read, type read2=None)
aligner = Aligner(read1=read_1, read2=read_2, index_directory=index_dir, output_directory= output_path, thread_number=thread_num, tool=alignutil)
align_dir = aligner.align()

# Convert them into sorted bam
converter = SamtoBamConversion(conversion_directory=align_dir, thread_number=thread_num)

# Create refflat and rna interval list files, can be used once like indexer
refflat_rnaint = RefflatRNAint(gtf_file=gtf_path, output_directory=align_dir, gtftogenepred_path=gtftogenepred_pr_path, data_source=data_origin,ref_flat_sc_path=refflat_script_path, rrna_int_list_sc_path=rrna_intlist_script_path)

# Check quality of each alignment, and collect rnaseq metrics
qc_controller = QCMetrics(qc_directory=align_dir, picard_path=picard_dir,thread_number=thread_num)

# Count aligned reads
counter = Quant(output_directory=output_path, gtf_file= gtf_path, alignment_directory= align_dir, thread_number=thread_num, tool=countutil, read_type=seqtype)










































