import os

def liftover(pathin, pathout):
    for file in os.listdir(pathin):
        if file.endswith(".vcf.gz") or file.endswith(".vcf"):

            rgsm = file.split('.')[0]
            inn = os.path.join(pathin, file)
            os.makedirs(pathout, exist_ok=True)
            output = os.path.join(pathout, f'{rgsm}_hg38.vcf.gz')
            rejected = os.path.join(pathout, f'{rgsm}_rejected.vcf.gz')
            cmd = f"gatk LiftoverVcf -R /home/genwork2/Mert/hg38.fa -I {inn} -C /home/genwork2/Mert/VariantDatabases/1.TurkishVariome/hg19ToHg38.over.chain.gz  --REJECT {rejected} -O {output}"
            os.system(cmd)

liftover("/home/genwork2/Mert/00.liftover", "/home/genwork2/Mert/00.liftover/hg38_liftover")


import vcf

input_vcf = "/home/genwork2/Mert/00.liftover/HE4188A_hg38.hard-filtered_fixed5.vcf"

input_vcf_path = "/home/genwork2/Mert/00.liftover/HE4188A_hg38.hard-filtered_fixed5.vcf"
output_vcf_path = "/home/genwork2/Mert/00.liftover/HE4188A_hg38.hard-filtered_fixed9.vcf"

with open(input_vcf_path, 'r') as vcf_file:
    vcf_reader = vcf.Reader(vcf_file)
    
    with open(output_vcf_path, 'w') as output_vcf_file:
        vcf_writer = vcf.Writer(output_vcf_file, vcf_reader)

        for record in vcf_reader:
            new_info = {key.strip('"'): value.strip('"') if isinstance(value, str) else value for key, value in record.INFO.items()}
            record.INFO = new_info

            vcf_writer.write_record(record)
