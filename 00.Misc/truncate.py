
import gzip
import zlib

file3 = "/home/genwork2/Mert/fastqfix/215-22-PPK-RE/215-22-PPK-RE_S205_L004_R1_001.fastq.gz"


try:
    with gzip.open(file3, 'rb') as file:
        trunk_d = file.read()
    print("File read successfully.")
except gzip.BadGzipFile:
    print("The file is truncated or invalid.")
except zlib.error:
    print("The file is truncated or invalid.")

