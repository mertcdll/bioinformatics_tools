import os

def check_files(directory):

    for filename in os.listdir(directory):
        if filename.endswith("fastq.gz") or filename.endswith("fq.gz"):
            if os.path.isfile(os.path.join(directory, filename)):
                print("Checking file:", os.path.join(directory, filename))
                print("Performing validation for file:", filename)
                command = f"fastQValidator --file {os.path.join(directory, filename)} --noeof"
                print("Command:", command)
                os.system(command)
        else:
            pass

directory_path = "/home/genwork2/Mert/fastqfix/DS2/2773"
check_files(directory_path)


