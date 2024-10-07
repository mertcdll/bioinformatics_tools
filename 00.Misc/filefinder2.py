import os
import time

class FileFinder:
    def __init__(self, directory, extension):
        self.directory = directory
        self.extension = extension

    def finder(self):
        start_time = time.time()  

        cmd = f"find {self.directory} -type f -name *{self.extension}"
        os.system(cmd)

        end_time = time.time() 
        elapsed_time = end_time - start_time
        print(f"Execution time: {elapsed_time:.2f} seconds")

directory = "/mnt/Gennas"
extension = "fastq.gz"

findall = FileFinder(directory, extension)
findall.finder()