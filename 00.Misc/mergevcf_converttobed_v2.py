import os
import subprocess
import pandas as pd
import re

class mergeVCF:
    def __init__(self, directory, output_file):
        self.directory = directory
        self.output_file = output_file

    def compress_vcf_files(self):
        
        vcf_files = [f for f in os.listdir(self.directory) if f.endswith('.vcf')] # get the list of vcfs

        
        for vcf_file in vcf_files: #compress vcf files if they exist
            input_file = os.path.join(self.directory, vcf_file)
            output_file = os.path.join(self.directory, vcf_file + '.gz')

            if os.path.exists(output_file): #if only there are vcf.gz, continue without  compressing
                continue  

            subprocess.run(['bgzip', '-c', input_file], stdout=open(output_file, 'wb'))

        self.create_tabix_index()

    def create_tabix_index(self): # tabix index for each vcf.gz
        
        vcf_files = [f for f in os.listdir(self.directory) if f.endswith('.vcf.gz')]

        for vcf_file in vcf_files:
            input_file = os.path.join(self.directory, vcf_file)
            subprocess.run(['tabix', '-p', 'vcf', input_file])

    def merge_vcf_files(self): #merge vcfs

        vcf_files = [f for f in os.listdir(self.directory) if f.endswith('.vcf.gz')]

        command = ['bcftools', 'merge', '-Ov', '-o', self.output_file] # create the command
        command.extend([os.path.join(self.directory, f) for f in vcf_files])

        subprocess.run(command)


class VCFConverterCountsToBed:
    def __init__(self, vcf_file, output_file):
        self.vcf_file = vcf_file
        self.output_file = output_file

    def CNVCountsToBed(self):
        cnv_dict = {}

        with open(self.vcf_file, "r") as vcf:
            for line in vcf:
                if not line.startswith("#"):
                    fields = line.strip().split("\t")
                    chrom = fields[0]
                    start_pos = int(fields[1])
                    end_pos = int(fields[7].split(";")[2].split("=")[1])

                    samples = fields[9:]
                    total_samples = len(samples)
                    wild_type_count = 0
                    deletion_count = 0
                    duplication_count = 0
                    # Calculate deletion and duplication counts
                    for sample in samples:
                        genotype_data = sample.split(':')[0]
                        if genotype_data == "./." or genotype_data == ".":
                            wild_type_count += 1
                        else:
                            genotype_values = genotype_data.split("/")
                            genotype_sum = sum(int(x) for x in genotype_values if x != '.')
                            copy_number = int(sample.split(':')[2])
                            if genotype_sum < copy_number:
                                duplication_count += 1
                            if genotype_sum >= copy_number:
                                deletion_count += 1

                    cnv_key = (chrom, start_pos, end_pos)
                    if cnv_key in cnv_dict:
                        cnv_dict[cnv_key]["Number of Deletions"] += deletion_count
                        cnv_dict[cnv_key]["Number of Duplications"] += duplication_count
                    else:
                        cnv_dict[cnv_key] = {
                            "Chromosome": chrom,
                            "Start Position": start_pos,
                            "End Position": end_pos,
                            "Number of Wild Types": wild_type_count,
                            "Number of Deletions": deletion_count,
                            "Number of Duplications": duplication_count,
                            "Total Number of Samples": total_samples
                        }

        cnv_data = list(cnv_dict.values())

        def custom_sorting_key(chromosome):  # sorting key for chromosomal order
            if isinstance(chromosome, str) and chromosome.startswith('chr'):
                number_match = re.match(r'chr(\d+)$', chromosome)
                if number_match:
                    return (int(number_match.group(1)), chromosome)
                else:
                    return (ord(chromosome[3:]), chromosome)
            return (0, chromosome)

        sorted_data = sorted(cnv_data, key=lambda x: (custom_sorting_key(x['Chromosome']), x['Start Position']))
        
        bed_lines = []
        for cnv in sorted_data:
            chrom = cnv['Chromosome']
            start = cnv['Start Position']
            end = cnv['End Position']
            deletions = cnv['Number of Deletions']
            duplications = cnv['Number of Duplications']
            total_samples = cnv['Total Number of Samples']

            merged_wt = total_samples - (deletions + duplications)
            bed_line = f"{chrom}\t{start}\t{end}\tWT:{merged_wt};DEL:{deletions};DUP:{duplications};TNS:{total_samples}"
            bed_lines.append(bed_line)

        with open(self.output_file, 'w') as file:
            file.write('\n'.join(bed_lines))

        compressed_file = self.output_file + '.gz'
        bgzip_command = ['bgzip', '-c', self.output_file]
        with open(compressed_file, 'wb') as f_out:
            subprocess.run(bgzip_command, stdout=f_out)


merger = mergeVCF('/home/dell/Documents', '/home/dell/Documents/merged.vcf') #merge first
merger.compress_vcf_files()
merger.merge_vcf_files()

converter = VCFConverterCountsToBed("/home/dell/Documents/merged.vcf", "/home/dell/Documents/output.bed")#convert to bed after
converter.CNVCountsToBed()