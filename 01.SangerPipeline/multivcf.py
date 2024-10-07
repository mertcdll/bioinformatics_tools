import vcf
import os
import gzip
import argparse

class MultiVcf:
    def __init__(self, input_folder, output_merged, merged_filename):
        self.input_folder = input_folder
        self.output_merged = output_merged
        self.merged_filename = merged_filename
        self.IlyomeComp()
        self.normalize_vcf()
        self.merge_vcfs()

    def IlyomeComp(self):
        for root, dirs, files in os.walk(self.input_folder):
            for filename in files:
                if filename.endswith(".vcf.gz"):
                    vcf_path = os.path.join(root, filename)
                    ab1vcf = vcf.Reader(filename=vcf_path)
                    
                    new_formats = [
                        vcf.parser._Format('AD', 1, 'String', 'Allelic depths for the ref and alt alleles in the order listed'),
                        vcf.parser._Format('AF', 1, 'Float', 'Allele fraction'),
                        vcf.parser._Format('DP', 1, 'Integer', 'Read depth'),
                        vcf.parser._Format('FT', 1, 'String', 'Filter Information')
                    ]

                    ab1vcf.formats.update({fmt.id: fmt for fmt in new_formats})

                    ab1vcf.formats.pop("GQ", None)

                    new_header_lines = []
                    for line in ab1vcf._header_lines:
                        if line.startswith("##FORMAT=<ID=GQ"):
                            continue
                        new_header_lines.append(line)

                    for fmt in new_formats:
                        new_header_lines.append('##FORMAT=<ID={},Number={},Type={},Description="{}">'.format(fmt.id, fmt.num, fmt.type, fmt.desc))

                    ab1vcf._header_lines = new_header_lines

                    ab1vcf.samples = [filename.replace(".vcf.gz", "")] #.split("_")[0]

                    output_file = os.path.join(root, filename.replace(".vcf.gz", "_ilyome.vcf.gz"))

                    file = gzip.open(output_file, "wt")
                    output_vcf = vcf.Writer(file, ab1vcf)
        
                    for record in ab1vcf:
                        record.FORMAT = "GT:AD:AF:DP:FT"
                        filter_info = "PASS" if record.FILTER == [] else record.FILTER
                        for sample in record.samples:
                            updated_data = sample.data._asdict()
                            updated_data["AD"] = f"{updated_data['GQ']},{updated_data['GQ']}"
                            if sample.data.GT == "0/1":
                                updated_data["AF"] = 0.5
                            else:
                                updated_data["AF"] = 1
                            updated_data["DP"] = 50
                            updated_data["FT"] = filter_info
                            del updated_data['GQ']
                            new_data = vcf.model.make_calldata_tuple(updated_data.keys())(**updated_data)
                            sample.data = new_data
                        
                        output_vcf.write_record(record)

                    output_vcf.close()

                    os.system(f"tabix -f -p vcf {output_file}")
                    
    def normalize_vcf(self):
        for root, dirs, files in os.walk(self.input_folder):
            for filename in files:
                if filename.endswith("_ilyome.vcf.gz"):
                    in_vcf_path = os.path.join(root, filename)
                    out_name = filename.replace("_ilyome.vcf.gz", "_norm.vcf.gz")
                    out_vcf_path = os.path.join(root, out_name)

                    os.system(f"bcftools norm -m- -O z -o {out_vcf_path} --threads 20 {in_vcf_path}")
                    os.system(f"tabix -f -p vcf {out_vcf_path}")

    def merge_vcfs(self):
        paths=[]

        for root, dirs, files in os.walk(self.input_folder):
            for filename in files:
                if filename.endswith("_norm.vcf.gz"):
                    vcf_path = os.path.join(root, filename)
                    paths.append(vcf_path)
        paths_str = " ".join(paths)
        output_merged_path = os.path.join(self.output_merged, f"{self.merged_filename}.vcf.gz")
        out_norm_merged_path = os.path.join(self.output_merged, f"{self.merged_filename}_mernorm.vcf.gz")

        os.system(f'bcftools merge {paths_str} -O z -o {output_merged_path}')
        os.system(f'bcftools norm -m- -O z -o {out_norm_merged_path} --threads 20 {output_merged_path}')
        os.system(f"tabix -f -p vcf {output_merged_path}")
        os.system(f"tabix -f -p vcf {out_norm_merged_path}")

    


def main():
    parser = argparse.ArgumentParser(description="Makes Sanger VCFs Ilyome Compatible and merges them")
    parser.add_argument("-i", "--input_folder", required=True, help="Input folder containing AB1 files")
    parser.add_argument("-o", "--output_merged", required=True, help="Output folder for variant call results")
    parser.add_argument("-fn", "--merged_filename", required=True, help="Path to the reference genome file")
    args = parser.parse_args()

    MultiVcf(args.input_folder, args.output_merged, args.merged_filename)

if __name__ == "__main__":
    main()



# import vcf


# vcffile = vcf.Reader(filename="/home/genwork2/Mert/ab1fastq/outputtest/ALPARSLAN_ALIK_REZTEP_SMARCD1_EX12_13F_D09/ALPARSLAN_ALIK_REZTEP_SMARCD1_EX12_13F_D09.vcf.gz")


# records = []
# for record in vcffile:
#     records.append(record)


# records[0].FILTER == []
    

# filter_info = "PASS" if records[0].FILTER == [] else records[0].FILTER

# filter_info