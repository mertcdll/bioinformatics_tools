import vcf
import gzip
import os
import multiprocessing

class FurtherTrim:
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.output_folder = output_folder

    def trimmer(self, vcf_file):
        input_path = os.path.join(self.input_folder, vcf_file)
        output_path = os.path.join(self.output_folder, f'{os.path.splitext(vcf_file)[0].replace(".vcf", "")}.ftrm.vcf.gz')

        with gzip.open(output_path, "wt") as f_writer:
            vcf_reader = vcf.Reader(filename=input_path)
            vcf_writer = vcf.Writer(f_writer, vcf_reader)
            
            for record in vcf_reader:
                fields_to_keep = {}
                for key in record.INFO:
                    if key in ("AC_joint", "AF_joint", "AN_joint", "nhomalt_joint", "AC_exomes", "AF_exomes", "AN_exomes", "nhomalt_exomes", "AC_genomes", "AF_genomes", "AN_genomes", "nhomalt_genomes"):
                        fields_to_keep[key] = record.INFO[key]
                record.INFO.clear()
                record.INFO.update(fields_to_keep)
                vcf_writer.write_record(record)
                
        print(f"{output_path} created successfully")

    def process_files(self, vcf_filesbwa):
        with multiprocessing.Pool() as pool:
            pool.map(self.trimmer, vcf_files)

if __name__ == "__main__":
    input_folder = "/mnt/gen100/gnomad4/gnomad4genomes/deciphered"
    output_folder = "/home/genwork2/Mert/gnomad4/genomesftrm"
    trimmer = FurtherTrim(input_folder, output_folder)
    files = os.listdir(input_folder)
    vcf_files = [vcf_file for vcf_file in files if vcf_file.endswith(".vcf.gz")]

    trimmer.process_files(vcf_files)
