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
            os.system(f"gunzip {self.fasta_file}")
            decompressed_fasta_file = self.fasta_file
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
        elif self.tool == 'bowtie2':
            os.system(f"bowtie2-build --large-index {self.fasta_file} {os.path.join(index_directory, 'genome_ind')} --threads {self.thread_number}")
        elif self.tool == 'segemehl':
            os.system(f"segemehl.x -x {os.path.join(index_directory, 'genome_ind.idx')} -d {decompressed_fasta_file} --threads {self.thread_number}")
        elif self.tool == 'gem':
            os.system(f"gem-indexer --threads {self.thread_number} -i {decompressed_fasta_file} -o {os.path.join(index_directory, 'genome_ind')}")
        elif self.tool == 'minimap2':
            os.system(f"minimap2 -d {os.path.join(index_directory, 'genome_ind.mmi')} {decompressed_fasta_file} -t {self.thread_number}")
        elif self.tool == 'magicblast':
            os.system(f"makeblastdb -in {decompressed_fasta_file} -out {os.path.join(index_directory, 'genome_ind')} -parse_seqids -dbtype nucl")
        
        return index_directory


class Aligner:
    def __init__(self, read1, read2, index_directory, fasta_file, output_directory, thread_number, tool):
        self.read1 = read1
        self.read2 = read2
        self.index_directory = index_directory
        self.fasta_file = fasta_file
        self.output_directory = output_directory
        self.thread_number = thread_number
        self.tool = tool

    def align(self):
        if self.read1.endswith('.gz'):
            os.system(f"gunzip {self.read1}")
            decompressedfile1 = self.read1
        else:
            decompressedfile1 = self.read1

        if self.read2.endswith('.gz'):
            os.system(f"gunzip {self.read2}")
            decompressedfile2 = self.read2
        else:
            decompressedfile2 = self.read2

        alignment_directory = os.path.join(self.output_directory, "alignment_output")
        os.makedirs(alignment_directory, exist_ok=True)

        if self.tool == 'star':
            os.system(f"STAR --runThreadN {self.thread_number} --genomeDir {self.index_directory} --readFilesIn {decompressedfile1} {decompressedfile2} --outFileNamePrefix {os.path.join(alignment_directory, '')}")    
        elif self.tool == 'hisat2':
            os.system(f"hisat2 -p {self.thread_number} -x {os.path.join(self.index_directory, 'genome_ind')} -1 {decompressedfile1} -2 {decompressedfile2} -S {os.path.join(alignment_directory, 'Aligned.out.sam')} > {os.path.join(alignment_directory, 'initial_metrics.txt')}")
        elif self.tool == 'bowtie2':
            os.system(f"bowtie2 --threads {self.thread_number} -x {os.path.join(self.index_directory, 'genome_ind')} -1 {decompressedfile1} -2 {decompressedfile2} -S {os.path.join(alignment_directory, 'Aligned.out.sam')} > {os.path.join(alignment_directory, 'initial_metrics.txt')}")
        elif self.tool == 'segemehl':
            os.system(f"segemehl.x --threads {self.thread_number} -i {os.path.join(self.index_directory, 'genome_ind.idx')} -d {self.fasta_file} -q {decompressedfile1} -p {decompressedfile2} --progressbar >  {os.path.join(alignment_directory, 'Aligned.out.sam')} ")
        elif self.tool == 'gem':
            os.system(f"gem-mapper -threads {self.thread_number} -I {os.path.join(self.index_directory, 'genome_ind.gem')} -1 {decompressedfile1} -2 {decompressedfile2} -o {os.path.join(alignment_directory, 'Aligned.out.sam')}")
        elif self.tool == 'minimap2':
            os.system(f"minimap2 -ax sr {os.path.join(self.index_directory, 'genome_ind.mmi')} {decompressedfile1} {decompressedfile2} -t {self.thread_number} > {os.path.join(alignment_directory, 'Aligned.out.sam')}")
        elif self.tool == 'magicblast':
            os.system(f"magicblast -query {decompressedfile1} -query_mate {decompressedfile2} -db {os.path.join(self.index_directory, 'genome_ind')} -num_threads {self.thread_number} -out {os.path.join(alignment_directory, 'Aligned.out.sam')} -infmt fastq -outfmt sam")
        
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

        

class QCMetrics:
    def __init__(self, qc_directory, ref_flat_path, rrna_int_list_path, picard_path, thread_number):
        self.qc_directory = qc_directory
        self.ref_flat_path = ref_flat_path
        self.rrna_int_list_path = rrna_int_list_path
        self.picard_path = picard_path
        self.thread_number = thread_number
        self.qc_control()
    
    def qc_control(self):

        os.system(f"samtools flagstat -@{self.thread_number} {os.path.join(self.qc_directory, 'final.bam')} > {os.path.join(self.qc_directory, 'flagstat_output.txt')}") 
        os.system(f"java -jar {self.picard_path} CollectRnaSeqMetrics -I {os.path.join(self.qc_directory, 'final.bam')} -O {os.path.join(self.qc_directory, 'output.collect_RNA_Metrics')} -REF_FLAT {self.ref_flat_path} -STRAND SECOND_READ_TRANSCRIPTION_STRAND -RIBOSOMAL_INTERVALS {self.rrna_int_list_path}")


class Quant:
    def __init__(self, output_directory, alignment_directory, gtf_file, thread_number, tool):
        self.output_directory = output_directory
        self.gtf_file = gtf_file
        self.tool = tool
        self.alignment_directory = alignment_directory
        self.thread_number = thread_number


    def quantification(self):

        read_count_dir = os.path.join(self.output_directory, "read_counts")
        os.makedirs(read_count_dir, exist_ok=True)

        if self.tool == "featurecounts":
            os.system(f"featureCounts -T {self.thread_number} -t exon -g gene_id -a {self.gtf_file} -o {os.path.join(read_count_dir, 'gene_level_counts.txt')} {os.path.join(self.alignment_directory, 'final.bam')} -p --countReadPairs")

        if self.tool == "htseq":
            os.system(f"htseq-count -f bam -r pos -i gene_id -m union -s no {os.path.join(self.alignment_directory, 'final.bam')} {self.gtf_file} > {os.path.join(read_count_dir, 'gene_level_counts.txt')}")    

        
class PseudoAlign:
    def __init__ (self, transcriptome_fasta_file, output_directory, read1, read2, gtf_file, thread_number, tool, bootstrap_value):
        self.transcriptome_fasta_file = transcriptome_fasta_file
        self.output_directory = output_directory
        self.read1 = read1
        self.read2 = read2
        self.gtf_file = gtf_file
        self.thread_number = thread_number
        self.tool = tool
        self.bootstrap_value = bootstrap_value
        self.pseudoaligner()

    def pseudoaligner(self):

        if self.transcriptome_fasta_file.endswith(".gz"):
            os.system(f"gunzip {self.transcriptome_fasta_file}")
            decompressed_fasta_file = self.transcriptome_fasta_file
        else:
            decompressed_fasta_file = self.transcriptome_fasta_file

        if self.read1.endswith('.gz'):
            os.system(f"gunzip {self.read1}")
            decompressedfile1 = self.read1
        else:
            decompressedfile1 = self.read1

        if self.read2.endswith('.gz'):
            os.system(f"gunzip {self.read2}")
            decompressedfile2 = self.read2
        else:
            decompressedfile2 = self.read2
 
        index_directory = os.path.join(self.output_directory, "transcriptome_index")
        os.makedirs(index_directory, exist_ok=True)

        alignment_directory = os.path.join(self.output_directory, "alignment_output")
        os.makedirs(alignment_directory, exist_ok=True)

        if self.tool == "kallisto": # gives transcript level abundances >> use tximport to obtain gene level abundances.
            os.system(f"kallisto index -i {os.path.join(index_directory, 'transcriptome_ind.idx')} {decompressed_fasta_file} --threads {self.thread_number}")
            os.system(f"kallisto quant -i {os.path.join(index_directory, 'transcriptome_ind.idx')} -o {alignment_directory} {decompressedfile1} {decompressedfile2} --threads {self.thread_number} --bootstrap-samples {self.bootstrap_value} --plaintext --genomebam --gtf {self.gtf_file}")
        

alignutil = "hisat2"
countutil = "featurecounts"
pseudoalignutil = "kallisto"
genome_fasta_path = "/home/genwork2/Mert/RNAseq/Homo_sapiens.GRCh38.dna.toplevel.fa" #segemehl needs in alignment also.
transcriptome_fasta_path =  "/home/genwork2/Mert/RNAseq/human_transcriptome.fa" 
gtf_path = "/home/genwork2/Mert/RNAseq/Homo_sapiens.GRCh38.110.gtf" #star needs in indexing stage.
output_path = "/home/genwork2/Mert/RNAseq"
thread_num = 40
bootstrap_val = 200
read_1 = "/home/genwork2/Mert/RNAseq/67-ZEYNEP_AKKUS_TH_S59_L003_R1_001.fastq"
read_2 = "/home/genwork2/Mert/RNAseq/67-ZEYNEP_AKKUS_TH_S59_L003_R2_001.fastq"
ref_flat_dir = "/home/genwork2/Mert/gtftorefflat/reFflat_ensembl_nohyp.txt"
rrna_intlist_dir = "/home/genwork2/Mert/gtftorefflat/rrna_int_list_ensembl.txt"
picard_dir = "/home/genwork2/Mert/apps/picard.jar"


indexer = Indexer(fasta_file=genome_fasta_path, gtf_file=None, output_directory=output_path, thread_number=thread_num, tool=alignutil)
index_dir = indexer.indexing()
aligner = Aligner(read1=read_1, read2=read_2, index_directory=index_dir, fasta_file=None, output_directory= output_path, thread_number=thread_num, tool=alignutil)
align_dir = aligner.align()
converter = SamtoBamConversion(conversion_directory=align_dir, thread_number=thread_num)
qc_controller = QCMetrics(qc_directory=align_dir, ref_flat_path=ref_flat_dir,rrna_int_list_path=rrna_intlist_dir,picard_path=picard_dir,thread_number=thread_num)
counter = Quant(output_directory=output_path, gtf_file= gtf_path, alignment_directory= align_dir, thread_number=thread_num, tool=countutil)



pseudoaligner = PseudoAlign(transcriptome_fasta_file=transcriptome_fasta_path, output_directory=output_path, read1=read_1, read2=read_2, gtf_file=gtf_path, thread_number=thread_num, tool=pseudoalignutil, bootstrap_value=bootstrap_val)











































