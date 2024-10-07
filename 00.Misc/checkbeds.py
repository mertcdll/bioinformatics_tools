import subprocess

class BEDChecker:
    def __init__(self, input_bed):
        self.input_bed = input_bed

    def bed_file_checker(self):
        try:
            result = subprocess.run(f"bedtools sort -i {self.input_bed}", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
            stderr_output = result.stderr.strip()
            if result.returncode !=0:
                print(f"Bed file is not valid: {stderr_output}")
            else: 
                print("Bed file is valid")
        except Exception as e:
            print(f"Error occurred in the validator script")


bed_file = "/mnt/Gennas/11.Bed_Files/48-PSMB5_gene/PSMB5_gene2.bed"

checker = BEDChecker(bed_file)

checker.bed_file_checker()