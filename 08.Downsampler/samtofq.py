import os


class Samtofq:

  def __init__(self, input_folder, refseq):
    self.input_folder = input_folder
    self.refseq = refseq
    

  def stfq(self):

    files = os.listdir(self.input_folder)
    ds = [file for file in files if file.startswith("DS_")]

    for downsampled in ds:

      ds_path = os.path.join(self.input_folder, downsampled)
      file_name = os.path.splitext(downsampled)[0].split("_")[-1]
      r1_name = f"{file_name}_20DS_R1_001.fastq.gz"
      r2_name = f"{file_name}_20DS_R2_001.fastq.gz"
      r1_path = os.path.join(self.input_folder, r1_name)
      r2_path = os.path.join(self.input_folder, r2_name)


      command = f"gatk SamToFastq -I {ds_path} -F {r1_path} -F2 {r2_path} -R {self.refseq}"
      os.system(command)


if __name__ == "__main__":
    input_folder = "/home/genwork2/Mert/downsamplebams/run174bams"  
    refseq = "/home/genwork2/Mert/hg38.fa"
    samtofq_converter = Samtofq(input_folder, refseq)
    samtofq_converter.stfq()