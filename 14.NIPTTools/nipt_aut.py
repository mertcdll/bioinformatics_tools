import os
import argparse
import subprocess

class NIPTBcl2Fastq:

    def __init__(self, input_folder, output_folder, num_readers, num_writers, num_processors, compression_level):
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.num_readers = num_readers
        self.num_writers = num_writers
        self.num_processors = num_processors
        self.compression_level = compression_level
        self.bclfolders = os.listdir(self.input_folder)
        self.create_folders()
        self.bcl_to_fastq()

    def create_folders(self):

        for folder in self.bclfolders:
            fastq_folder = os.path.join(self.output_folder, folder)
            os.makedirs(fastq_folder, exist_ok=True)

    def bcl_to_fastq(self):

        for bcl_folder in self.bclfolders:
            bcl_abs_path = os.path.join(self.input_folder, bcl_folder)
            fastq_abs_path = os.path.join(self.output_folder, bcl_folder)
            samplesheet_path = os.path.join(bcl_abs_path, "SampleSheet.csv")
                        
            run_command = (
                f"nohup bcl2fastq -R {bcl_abs_path} -o {fastq_abs_path} "
                f"--sample-sheet {samplesheet_path} -r {self.num_readers} "
                f"-w {self.num_writers} -p {self.num_processors} "
                f"--fastq-compression-level {self.compression_level} "
                "--no-lane-splitting &"
            )

            subprocess.run(run_command, shell=True)


def parse_args():
    parser = argparse.ArgumentParser(description="Run bcl2fastq conversion for multiple BCL folders -- To run it, please make sure that only the bcl files which are going to be converted are inside the input folder")
    parser.add_argument("-i", "--input-folder", required=True, help="Path to the input BCL folder")
    parser.add_argument("-o", "--output-folder", required=True, help="Path to the output Fastq folder")
    parser.add_argument("-r", "--num-readers", type=int, default=5, help="Number of readers")
    parser.add_argument("-w","--num-writers", type=int, default=5, help="Number of writers")
    parser.add_argument("-p","--num-processors", type=int, default=20, help="Number of processors")
    parser.add_argument("-cl","--compression-level", type=int, default=8, help="Compression level")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    nipt = NIPTBcl2Fastq(
        args.input_folder, args.output_folder,
        args.num_readers, args.num_writers,
        args.num_processors, args.compression_level
    )




