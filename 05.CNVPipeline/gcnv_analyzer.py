import os
from multiprocessing import Pool



class GCNV:
    def __init__(self, bam_dir, exome_loc, genome_fasta, ref_dict, contig_pp, num_processes, class_coherence_length, cnv_coherence_length, min_contig_length):
        self.bam_dir = bam_dir
        self.exome_loc = exome_loc
        self.genome_fasta = genome_fasta
        self.ref_dict = ref_dict
        self.contig_pp = contig_pp
        self.num_processes = num_processes
        self.class_coherence_length = class_coherence_length
        self.cnv_coherence_length = cnv_coherence_length
        self.min_contig_length = min_contig_length
        self.input_string = None
        self.filtered_il = None
        self.interval_list = None
        self.ann_interval_list = None
        self.bam_list = None
        self.sample_list = None

    def preprocess_intervals(self):
        gatkt = "PreprocessIntervals"
        imr = "OVERLAPPING_ONLY"

        self.interval_list = os.path.join(self.bam_dir, "hg38_exome.preprocessed.interval_list")
        command = f"gatk {gatkt} -R {self.genome_fasta} -L {self.exome_loc} --bin-length 0 -imr {imr} -O {self.interval_list}"
        os.system(command)

    def annotate_intervals(self):
        gatkt = "AnnotateIntervals"
        imr = "OVERLAPPING_ONLY"

        self.ann_interval_list = os.path.join(self.bam_dir, "hg38_exome.annotated.tsv")
        command = f"gatk {gatkt} -L {self.interval_list} -R {self.genome_fasta} -imr {imr} -O {self.ann_interval_list}"    
        os.system(command)


    def run_gatk_collect_rc(self, args):
        bam_file, interval_list, genome_fasta = args
        sample_name = os.path.splitext(bam_file)[0]
        output_name = sample_name + "_counts.tsv"
        gatkt = "CollectReadCounts"
        imr = "OVERLAPPING_ONLY"

        command = f"gatk {gatkt} -L {interval_list} -R {genome_fasta} -imr {imr} -I {os.path.join(bam_dir, bam_file)} -O {os.path.join(bam_dir, output_name)}" 
        os.system(command)


    def count_reads(self):

        self.bam_list = [item for item in os.listdir(self.bam_dir) if item.endswith(".bam") or item.endswith(".cram")] #1

        arguments = [(bam_file, self.interval_list, self.genome_fasta) for bam_file in self.bam_list]


        try: 
            with Pool(self.num_processes) as pool:
                pool.map(self.run_gatk_collect_rc, arguments)

        except Exception as e:
            print(f"error: {e}")


    def get_input_string(self):

        input_string = ""

        self.clist = [os.path.join(self.bam_dir,item) for item in os.listdir(self.bam_dir) if item.endswith("_counts.tsv")] #1

        for path in self.clist:
            input_string += f"-I {path} "

        self.input_string = input_string.strip()

        return self.input_string

    def filterintervals(self):

        self.filtered_il = os.path.join(self.bam_dir, "cohort_filtered.interval_list")
        gatkt = "FilterIntervals"
        imr = "OVERLAPPING_ONLY"
        command = f"gatk {gatkt} -L {self.interval_list} --annotated-intervals {self.ann_interval_list} {self.input_string} -imr {imr} -O {self.filtered_il}"
        os.system(command)

    def determinecontig(self):

        gatkt = "DetermineGermlineContigPloidy"
        imr = "OVERLAPPING_ONLY"
        command = f"gatk {gatkt} -L {self.filtered_il} -imr {imr} {self.input_string} --contig-ploidy-priors {self.contig_pp} --output {self.bam_dir} --output-prefix ploidy --verbosity DEBUG"
        os.system(command)
    
    def CNVCaller(self):

        gatkt = "GermlineCNVCaller"
        imr = "OVERLAPPING_ONLY"
        ploidy_calls = os.path.join(self.bam_dir, "ploidy-calls")
        output = os.path.join(self.bam_dir, "cnv-cohort")
        coherence_lengths = f"--class-coherence-length {self.class_coherence_length} --cnv-coherence-length {self.cnv_coherence_length}"
        command = f"gatk {gatkt} --run-mode COHORT -L {self.filtered_il} -imr {imr} {self.input_string} --contig-ploidy-calls {ploidy_calls} --annotated-intervals {self.ann_interval_list} --output {output} --output-prefix cnv-cohort --verbosity DEBUG {coherence_lengths}"                   
        os.system(command)


    def postprocess_sample(self):
        
        gatkt = "PostprocessGermlineCNVCalls"
        imr = "OVERLAPPING_ONLY"
        ploidy_calls = os.path.join(self.bam_dir, "ploidy-calls")
        model = os.path.join(self.bam_dir, "cnv-cohort", "cnv-cohort-model")
        calls = os.path.join(self.bam_dir, "cnv-cohort", "cnv-cohort-calls")
        

        for i in range(len(self.clist)): 
            with open(os.path.join(calls, f"SAMPLE_{i}", "sample_name.txt")) as file:
                name = file.read()
                name = name.strip()
            
            interval_output = os.path.join(self.bam_dir, f"{name}_intervals.vcf.gz")
            segment_output = os.path.join(self.bam_dir, f"{name}_segments.vcf.gz")
            denoised = os.path.join(self.bam_dir, f"{name}_denoised_copy_ratios.tsv")

            command = f"gatk {gatkt} -imr {imr} --model-shard-path {model} --calls-shard-path {calls} --allosomal-contig chrX --allosomal-contig chrY --contig-ploidy-calls {ploidy_calls} --sample-index {i} --output-genotyped-intervals {interval_output} --output-genotyped-segments {segment_output} --output-denoised-copy-ratios {denoised} --sequence-dictionary {self.ref_dict}"
            os.system(command)

    def panel_of_normals(self):

        gatkt = "CreateReadCountPanelOfNormals"
        output = os.path.join(self.bam_dir, "cnvpon.hdf5")
        opts = "--minimum-interval-median-percentile 5.0 --spark-master local[40]"

        command = f"gatk {gatkt} {self.input_string} {opts} -O {output}"

        os.system(command)

        return(output)

    def stardardize_denoise(self, pon):

        gatkt = "DenoiseReadCounts"
        counts_list = [item for item in os.listdir(self.bam_dir) if item.endswith("_counts.tsv")]

        cr_list = []
        for item in counts_list:
            count_path = os.path.join(self.bam_dir, item)
            item_prefix = item.split("_counts.tsv")[0]
            sdcr_output = os.path.join(self.bam_dir, f"{item_prefix}_standardizedCR.tsv")
            dncr_output = os.path.join(self.bam_dir, f"{item_prefix}_denoisedCR.tsv")
            cr_output = [item_prefix, sdcr_output, dncr_output]

            cpon = f"--count-panel-of-normals {pon}"
            sdcr = f"--standardized-copy-ratios {sdcr_output}"
            dncr = f"--denoised-copy-ratios {dncr_output}"

            command = f"gatk {gatkt} -I {count_path} {cpon} {sdcr} {dncr}"

            cr_list.append(cr_output)

            os.system(command)
        
        return(cr_list)
    
    def plot_pt1(self, copy_ratio_paths):

        gatkt = "PlotDenoisedCopyRatios"
        min_contig = f"--minimum-contig-length {self.min_contig_length}"
        output = f"--output {os.path.join(self.bam_dir, 'plots')}"
        seq_dict = f"--sequence-dictionary {self.ref_dict}"
        
        for item in copy_ratio_paths:
            out_prefix = f"--output-prefix {item[0]}"
            sd_cr = f"--standardized-copy-ratios {item[1]}"
            dn_cr = f"--denoised-copy-ratios {item[2]}"
            
            command = f"gatk {gatkt} {sd_cr} {dn_cr} {seq_dict} {min_contig} {output} {out_prefix}"

            os.system(command)

    def run_gatk_collect_ac(self, args):
        bam_file, interval_list, genome_fasta = args
        sample_name = os.path.splitext(bam_file)[0]
        output_name = sample_name + "_allelic_counts.tsv"
        gatkt = "CollectAllelicCounts"
        imr = "OVERLAPPING_ONLY"
        command = f"gatk {gatkt} -L {interval_list} -R {genome_fasta} -imr {imr} -I {os.path.join(self.bam_dir, bam_file)} -O {os.path.join(self.bam_dir, output_name)}" 
        os.system(command)


    def allelic_counts(self):

        arguments = [(bam_file, self.filtered_il, self.genome_fasta) for bam_file in self.bam_list]

        try: 
            with Pool(self.num_processes) as pool:
                pool.map(self.run_gatk_collect_ac, arguments)

        except Exception as e:
            print(f"error: {e}")

    def run_gatk_model_segments(self, args):
        sample_name, sample_path, bam_dir = args
        gatkt = "ModelSegments"
        output_dir = os.path.join(bam_dir, 'segmentation_results')

        output_prefix = f"--output-prefix {sample_name}"
        denoised = f"--denoised-copy-ratios {sample_path}_denoisedCR.tsv"
        ac = f"--allelic-counts {sample_path}_allelic_counts.tsv"

        command = f"gatk {gatkt} {denoised} {ac} --output {output_dir} {output_prefix}"
        os.system(command)


    def model_segments(self):

        self.sample_list = [os.path.splitext(item)[0] for item in os.listdir(self.bam_dir) if item.endswith(".bam") or item.endswith(".cram")]

        arguments = [(sample_name, os.path.join(self.bam_dir, sample_name), self.bam_dir) for sample_name in self.sample_list]
        
        try:
            with Pool(self.num_processes) as pool:
                pool.map(self.run_gatk_model_segments, arguments)

        except Exception as e:
            print(f"error: {e}")

        return os.path.join(self.bam_dir, 'segmentation_results')
       
    def copy_ratio_segments(self, segment_folder):

        gatkt = "CallCopyRatioSegments"
        input_list = [item for item in os.listdir(segment_folder) if item.endswith(".cr.seg")]
        
        for item in input_list:
            input = f"--input {os.path.join(segment_folder, item)}"
            sample_name = item.split(".cr.seg")[0]
            output= f"--output {os.path.join(segment_folder, f'{sample_name}.called.seg')}"

            command = f"gatk {gatkt} {input} {output}"
            os.system(command)


    def plot_segments(self, segment_folder):

        gatkt = "PlotModeledSegments"
        output = os.path.join(self.bam_dir, "plots_seg")

        for item in self.sample_list:
            denoised = f"--denoised-copy-ratios {os.path.join(self.bam_dir, f'{item}_denoisedCR.tsv')}"
            ac = f"--allelic-counts {os.path.join(segment_folder, f'{item}.hets.tsv')}"
            segments = f"--segments {os.path.join(segment_folder, f'{item}.modelFinal.seg')}"
            seq_dict= f"--sequence-dictionary {self.ref_dict}"
            min_contig_l = f"--minimum-contig-length {self.min_contig_length}"
            output_folder = f"--output {output}"
            output_prefix = f"--output-prefix {item}"

            command = f"gatk {gatkt} {denoised} {ac} {segments} {seq_dict} {min_contig_l} {output_folder} {output_prefix}"
            os.system(command)



if __name__ == "__main__":
    ref_dict = "/home/genwork2/Mert/hg38.dict"
    bam_dir = "/home/genwork2/Mert/CNVpip/newtrial"
    exome_loc = "/home/genwork2/Mert/CNVpip/hg38_exome_v2.0.2_targets_sorted_validated.annotated_0.bed"
    genome_fasta = "/home/genwork2/Mert/hg38.fa"
    contig_pp = "/home/genwork2/Mert/CNVpip/nesdoc/contig_ploidy_prior_full.tsv"
    class_coherence_length = 10000000
    cnv_coherence_length = 10000000
    num_processes = 40
    min_contig_length = 46709983

    counts = GCNV(bam_dir, exome_loc, genome_fasta, ref_dict, contig_pp, num_processes, class_coherence_length, cnv_coherence_length, min_contig_length)

    counts.preprocess_intervals()
    counts.annotate_intervals()
    counts.count_reads()
    counts.get_input_string()
    counts.filterintervals()
    counts.determinecontig()
    counts.CNVCaller()
    counts.postprocess_sample()
    pon = counts.panel_of_normals()
    cr_list = counts.stardardize_denoise(pon)
    counts.plot_pt1(cr_list)
    counts.allelic_counts()
    seg_dir = counts.model_segments()
    counts.copy_ratio_segments(seg_dir)
    counts.plot_segments(seg_dir)