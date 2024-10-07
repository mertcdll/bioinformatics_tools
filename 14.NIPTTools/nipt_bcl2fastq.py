import os
import argparse
import subprocess
import ftplib

class NIPTBcl2Fastq:

    def __init__(self, output_folder,
                 num_readers, num_writers, num_processors, 
                 compression_level, fastq_file_names, nipt_part, nipt_date, additional_folder=None):
        
        self.output_folder = output_folder
        self.num_readers = num_readers
        self.num_writers = num_writers
        self.num_processors = num_processors
        self.compression_level = compression_level
        self.nipt_part = nipt_part
        self.fastq_file_names = fastq_file_names
        self.nipt_date = nipt_date
        self.additional_folder = additional_folder
        self.bclfolders =  []
        self.fastq_folder_paths = []
        self.fastq_file_paths = []
        self.copybclfromnas()
        self.create_folders()
        self.bcl_to_fastq()
        self.take_sample_paths()
        self.upload_file_secure()
        
    def copybclfromnas(self):
        
        print("Copying bcl files to server")
        self.bclfolders = [folder for folder in os.listdir("/mnt/nextseq/") if folder.startswith(str(self.nipt_date))]
        
        if self.bclfolders:
            source_paths = " ".join([f"/mnt/nextseq/{folder}" for folder in self.bclfolders])
            command = f"cp -r {source_paths} /home/genwork2/04.BCL_NIPT"
            
            subprocess.run(command, shell=True)
        else:
            print("No BCL folders found to copy.")
    
    def create_folders(self):
        
        self.bclfolders = [folder for folder in os.listdir("/home/genwork2/04.BCL_NIPT") if folder.startswith(str(self.nipt_date))]

        print("Creating fastq folders")
        for folder in self.bclfolders:
            fastq_folder = os.path.join(self.output_folder, folder)
            os.makedirs(fastq_folder, exist_ok=True)
            self.fastq_folder_paths.append(fastq_folder)
        print("Fastq folders created")
        
    def bcl_to_fastq(self):
        
        print("Starting bcl2fastq")
        
        self.fastq_folder_paths = []
        
        for bcl_folder in self.bclfolders:
            bcl_abs_path = os.path.join("/home/genwork2/04.BCL_NIPT", bcl_folder)
            fastq_abs_path = os.path.join(self.output_folder, bcl_folder)
            samplesheet_path = os.path.join(bcl_abs_path, "SampleSheet.csv")
            
            self.fastq_folder_paths.append(fastq_abs_path)
                        
            run_command = (
                f"bcl2fastq -R {bcl_abs_path} -o {fastq_abs_path} "
                f"--sample-sheet {samplesheet_path} -r {self.num_readers} "
                f"-w {self.num_writers} -p {self.num_processors} "
                f"--fastq-compression-level {self.compression_level} "
                "--no-lane-splitting"
            )

            subprocess.run(run_command, shell=True)

    def take_sample_paths(self):

        print("Taking sample paths")
        self.bclfolders = [folder for folder in os.listdir("/home/genwork2/04.BCL_NIPT") if folder.startswith(str(self.nipt_date))]
        if self.additional_folder:
            all_folders = self.bclfolders + [self.additional_folder]
        else:
            all_folders = self.bclfolders
        
        all_folder_paths = [os.path.join(self.output_folder, folder) for folder in all_folders]
        print(all_folder_paths)
        self.fastq_file_paths = []
        
        for folder_path in all_folder_paths:
            for root, _, files in os.walk(folder_path):
                for file in files:
                    for file_name in self.fastq_file_names:
                        if file.startswith(file_name):
                            self.fastq_file_paths.append(os.path.join(root, file))

        print("Took sample paths")
        
    def upload_file_secure(self, ftp_server="", ftp_user="", ftp_password=""):
        
        print("Started uploading the files")
        
        try:
            ftps = ftplib.FTP_TLS()
            ftps.connect(ftp_server, 21)
            ftps.login(user=ftp_user, passwd=ftp_password)
            ftps.prot_p()
            ftps.set_pasv(True) 
            counter = 0
            for file_path in self.fastq_file_paths:
                folder_path = file_path.split("/")[-2]
                fastq_name = os.path.basename(file_path)
                remote_path = os.path.join(f"NIPT-Part{str(self.nipt_part)}", folder_path, fastq_name)
                with open(file_path, 'rb') as file:
                    ftps.storbinary(f'STOR {remote_path}', file)
                    print(f"File uploaded successfully: {remote_path}")
                    counter += 1
            print(f"In total {counter} Fastq files uploaded to the server")
            ftps.quit()
        except ftplib.all_errors as e:
            print(f"FTP error: {e}")


def parse_args():
    parser = argparse.ArgumentParser(description="Run bcl2fastq conversion for multiple BCL folders, upload resulting fastqs to the server.")
    parser.add_argument("-o", "--output-folder", required=True, help="Path to the output Fastq folder")
    parser.add_argument("-r", "--num-readers", type=int, default=10, help="Number of readers")
    parser.add_argument("-w","--num-writers", type=int, default=10, help="Number of writers")
    parser.add_argument("-p","--num-processors", type=int, default=40, help="Number of processors")
    parser.add_argument("-cl","--compression-level", type=int, default=8, help="Compression level")
    parser.add_argument("-names","--fastq-names", nargs="+", help="Fastq sample names - For example 24B3043312", required=True)
    parser.add_argument("-part", "--nipt-part", type=int, help="Part Number", required=True)
    parser.add_argument("-date", "--nipt-date", type=int, help= "Nipt Date", required=True)
    parser.add_argument("-add-folder", "--additional-folder", type=str, help="Additional folder name if present")
    
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    nipt = NIPTBcl2Fastq(
        args.output_folder,
        args.num_readers, args.num_writers,
        args.num_processors, args.compression_level, 
        args.fastq_names, args.nipt_part, args.nipt_date, args.additional_folder
    )


