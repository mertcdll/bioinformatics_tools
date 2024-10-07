import os
import argparse

class MD5SumCalculator:
    def __init__(self, directory, extension):
        self.directory = directory
        self.extension = extension

    def calculate_md5sum(self, filename):
        os.system(f"md5sum {filename} > {filename}.md5sum")

    def process_directory(self):
        for filename in os.listdir(self.directory):
            if filename.endswith(self.extension):
                file_path = os.path.join(self.directory, filename)
                if os.path.isfile(file_path):
                    print(f"Processing {file_path} ...")
                    self.calculate_md5sum(file_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='VCFMD5CHECKER', description='Script to generate md5sum of vcfs')
    parser.add_argument('--dir', metavar="--directory", required=True, type=str, help="The full path to the target folder")
    parser.add_argument('--ext', metavar="--extension", required=True, type=str, help="Extension of the file")
    args = parser.parse_args()
    
    if args:
        md5sum_calculator = MD5SumCalculator(args.dir, args.ext)
        md5sum_calculator.process_directory()