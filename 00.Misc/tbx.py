import os
from multiprocessing import Pool

def process_file(file):
    os.system(f"sudo gunzip {file}")
    #os.system(f"bgzip -@ 40 {vcf_file}")
    #os.system(f"tabix -p vcf {vcf_gz_file}")


def process_vcf_files(directory):
    file_list = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith('.vcf.gz')]

    with Pool() as pool:
        pool.map(process_file, file_list)

if __name__ == "__main__":
    directory_to_search = '/mnt/gen100/gnomad4/gnomad4exomes/deciphered'
    process_vcf_files(directory_to_search)




    