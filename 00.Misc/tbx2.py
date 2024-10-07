import os

def decomptab(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".vcf.gz"):
                vcf_gz_name = os.path.join(root, file)
                vcf_name = file.replace(".vcf.gz", "")
                vcf_bgz_name = vcf_name + ".bg.vcf.gz"

                os.system(f"zcat {vcf_gz_name} | bgzip -c > {os.path.join(root, vcf_bgz_name)}")

                os.system(f"tabix -p vcf {os.path.join(root, vcf_bgz_name)}")


decomptab("/mnt/gen100/gnomad4_rep")
