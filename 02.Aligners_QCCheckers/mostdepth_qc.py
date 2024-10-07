import os
import pandas as pd
import subprocess

class CoverageAnalysis:
    def __init__(self, reference, output_dir, max_worker=40, bed=None):
        self.reference = reference
        self.bed = bed
        self.output_dir = output_dir
        self.max_worker = max_worker
    
    def run_mosdepth(self, bam, output_prefix, wgs=False):
        output = os.path.join(self.output_dir, output_prefix)
        
        if wgs:
            command = [
                "mosdepth", output, f"-x -Q 1 -t {self.max_worker}",
                f"-f {self.reference}",
                f"-n",
                bam
            ]
            subprocess.run(" ".join(command), shell=True)
        else:
            command = [
                "mosdepth", output, f"-x -Q 1 -t {self.max_worker}",
                "--by", self.bed, f"-f {self.reference}",
                f"-n -T 1,5,10,15,20,30,50,100,500,1500",
                bam
            ]
            subprocess.run(" ".join(command), shell=True)
    
    def calculate_coverages(self, qc_thresholds_file):
        result = []
        df = pd.read_csv(qc_thresholds_file, compression='gzip', header=None, sep='\t')
        columns = df.iloc[0].to_list()
        df = df.drop(df.index[0])
        df.columns = columns
        df['size'] = df['end'].astype(int) - df['start'].astype(int)
        
        total_size = df['size'].sum()
        print(total_size)
        
        # Coverage columns
        coverage_columns = [col for col in df.columns if col.endswith('X')]
        
        # Calculate the sum of each coverage column and then calculate coverage
        coverage_dict = {}
        for col in coverage_columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').astype(int)
            coverage_sum = df[col].sum()
            coverage = coverage_sum / total_size
            coverage_dict[col] = coverage
        
        # Display the coverage values
        for col, cov in coverage_dict.items():
            result.append({"key": f"PCT of QC coverage region with coverage [{col.lower()}: inf):", "value": round(cov * 100, 2)})
        
        coverage_keys = list(coverage_dict.keys())
        for i in range(len(coverage_keys) - 1):
            diff = df[coverage_keys[i]] - df[coverage_keys[i + 1]]
            result.append({"key": f"PCT of QC coverage region with coverage [{coverage_keys[i].lower()}:{coverage_keys[i + 1].lower()}):", "value": round((diff.sum() / total_size) * 100, 2)})
        
        return result

    def OnTargetContigs(self, depth_summary):
        df = pd.read_csv(depth_summary, sep="\t")
        df = df.iloc[:50]
        # Select values without _region. all chromosomes
        df1 = df[~df['chrom'].str.endswith('_region')]
        total_bases = df1['bases'].astype(int).sum()

        # Select target regions only
        df = df[df['chrom'].str.endswith('_region')]
        df['chrom'] = df['chrom'].str.replace('_region', '')
        target_total_bases = df['bases'].astype(int).sum()
        target_length = df['length'].astype(int).sum()
        mean_coverage = round(target_total_bases / target_length, 2)

        df = df.drop(['length', 'bases', 'min', 'max'], axis=1)

        df.columns = ['key', 'value']
        # uniformity_2 = sum(1 for x in total if x > 0.2 * mean_coverage) / total_bases * 100
        # uniformity_4 = sum(1 for x in total if x > 0.4 * mean_coverage) / total_bases * 100
        final_result = [
            {"key": "Aligned bases", "value": total_bases},
            {"key":"Aligned bases in QC coverage region", "value":target_total_bases},
            {"key": "Average alignment coverage over QC coverage region", "value": mean_coverage},
            # {"key": "Uniformity of coverage (PCT > 0.2*mean)", "value": uniformity_2},
            # {"key": "Uniformity of coverage (PCT > 0.4*mean)", "value": uniformity_4},
            {"key": "PCT of QC coverage region with coverage [0x: inf):", "value": 100}
        ]
        return df.to_dict(orient="records"), final_result


    def process_directory(self, input_dir, wgs=False):
        for root, dirs, files in os.walk(input_dir):
            for file in files:
                if file.endswith(".bam"):
                    bam_path = os.path.join(root, file)
                    output_prefix = os.path.splitext(file)[0]
                    self.run_mosdepth(bam_path, output_prefix, wgs)
                    
                    qc_thresholds_file = os.path.join(self.output_dir, f"{output_prefix}.thresholds.bed.gz")
                    depth_summary_file = os.path.join(self.output_dir, f"{output_prefix}.mosdepth.summary.txt")
                    
                    if os.path.exists(qc_thresholds_file):
                        metrics = self.calculate_coverages(qc_thresholds_file)
                        
                        if os.path.exists(depth_summary_file):
                            _, final_results = self.OnTargetContigs(depth_summary_file)
                            metrics.extend(final_results)

                        output_metrics_file = os.path.join(self.output_dir, f"{output_prefix}_metrics.txt")
                        with open(output_metrics_file, 'w') as f:
                            for metric in metrics:
                                f.write(f"{metric['key']} {metric['value']}\n")


reference_genome = '/home/genwork2/Mert/01.REFGEN/GRCh37/hg19.fa'
bed_file = '/home/genwork2/Mert/02.BEDS/hg19_chrbed.bed'
output_directory = '/home/genwork2/Mert/NIPT/nipt_metricshg19/mosdepth_metrics'
input_directory = '/home/genwork2/Mert/NIPT/NIPT_BAMS_HG19'
wgs = False

coverage_analysis = CoverageAnalysis(reference=reference_genome, output_dir=output_directory, bed = bed_file)
coverage_analysis.process_directory(input_dir=input_directory, wgs=wgs)


