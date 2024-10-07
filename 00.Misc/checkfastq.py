import os
import gzip
import zlib
from multiprocessing import Pool


class FastqFileChecker:
    def __init__(self, fastq_info):
        self.fastq_info = fastq_info
        self.r1_path = fastq_info.get("filefastq").get("r1")
        self.r2_path = fastq_info.get("filefastq").get("r2")

    def check_existence(self, file_path):
        return os.path.isfile(file_path)

    def check_reads(self):
        r1_name = self.r1_path.split("/")[-1]
        r2_name = self.r2_path.split("/")[-1]
        
        if "_R1_" in r1_name or "1.fq" in r1_name:
            pass
        else:
            raise AssertionError(f"Invalid R1 file name:{r1_name}")
        
        if "_R2_" in r2_name or "2.fq" in r2_name:
            pass
        else:
            raise AssertionError(f"Invalid R2 file name:{r2_name}")
        

    def check_truncation(self, file_path):
        try:
            with gzip.open(file_path, 'rb') as file:
                file.read()
            return False
        except (gzip.BadGzipFile, zlib.error):
            return True
        except IOError:
            return True

    def check_files(self):

        
        r1_existence = self.check_existence(self.r1_path)
        r1_truncation = self.check_truncation(self.r1_path)

        r2_existence = self.check_existence(self.r2_path)
        r2_truncation = self.check_truncation(self.r2_path)

        self.fastq_info["filefastq"]["r1_existence"] = r1_existence
        self.fastq_info["filefastq"]["r1_truncation"] = r1_truncation
        self.fastq_info["filefastq"]["r2_existence"] = r2_existence
        self.fastq_info["filefastq"]["r2_truncation"] = r2_truncation

        return self.fastq_info


def process_fastq_info(fastq_info):
    fastq_checker = FastqFileChecker(fastq_info)
    fastq_checker.check_reads()
    result_dict = fastq_checker.check_files()
    return result_dict

if __name__ == "__main__":
    fastq_info1 = {'id': 8913,
                'sampleid': 'Sample 8071',
                'filefastq': {'id': 8071, 
                'r1': '/home/genwork2/Mert/fastqfix/215-22-PPK-RE_S205_L004_R1_001.fastq.gz',
                'r2': '/home/genwork2/Mert/fastqfix/215-22-PPK-RE_S205_L004_R2_001.fastq.gz'}}



    
    fastq_info_list = [fastq_info1]
    
    num_processes = len(fastq_info_list)

    with Pool(num_processes) as pool:
        results = pool.map(process_fastq_info, fastq_info_list)

    for result in results:
        print(result)


