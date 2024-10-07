import os
import argparse
import subprocess

class AB1VariantCall:
    def __init__(self, input_folder, output_folder, reference_genome, trim_left, trim_right):
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.reference_genome = reference_genome
        self.trim_left = trim_left
        self.trim_right = trim_right
        self.variant_caller()
        self.bcf2vcf()

    def variant_caller(self):
        sample_names = []
        
        for filename in os.listdir(self.input_folder):
                if filename.endswith("ab1"):
                    ab1_path = os.path.join(self.input_folder, filename)
                    outprefix = os.path.splitext(filename)[0]
                    sample_names.append(outprefix)
                    out_dir = os.path.join(self.output_folder, "variantcalls", outprefix)
                    os.makedirs(out_dir, exist_ok=True)
                    outprefix_path = os.path.join(out_dir, outprefix)
                    command = ["tracy", "decompose", "-v", "-a", "homo_sapiens", "-r", self.reference_genome, "-o", outprefix_path]

                    if self.trim_left is not None:
                        command.extend(["-q", str(self.trim_left)])
                    if self.trim_right is not None:
                        command.extend(["-u", str(self.trim_right)])

                    command.append(ab1_path)

                    try:
                        result = subprocess.run(command, capture_output=True, text=True, check=True)
                        stdout = result.stdout
                        stderr = result.stderr
                        if stdout:
                            print("Command stdout:", stdout)
                        if stderr:
                            print("Command stderr:", stderr)
                    except subprocess.CalledProcessError as e:
        
                        print(f"Error executing command: {e}")
                        
                        if e.stderr:
                            print("Error Message:", e.stderr) # AB1 file either lacks basecalls or Sanger trace could not be anchored to reference genome. You can add fail or success here, while adding a new variable. 

        
        return sample_names      
    
    def bcf2vcf(self):
        for root, dirs, files in os.walk(self.output_folder):
            for filename in files:
                if filename.endswith("bcf"):
                    fn = os.path.splitext(filename)[0]
                    vcf_name = fn + ".vcf.gz"
                    bcf_path = os.path.join(root, filename)
                    vcf_path = os.path.join(root, vcf_name)

                    command = f"bcftools convert -O z -o {vcf_path} {bcf_path}"

                    print(f"Converting {filename} to compressed vcf")
                    os.system(command)

                    os.system(f"tabix -p vcf {vcf_path}")

def main():
    parser = argparse.ArgumentParser(description="AB1 Variant Caller")
    parser.add_argument("-i", "--input_folder", required=True, help="Input folder containing AB1 files")
    parser.add_argument("-o", "--output_folder", required=True, help="Output folder for variant call results")
    parser.add_argument("-r", "--reference_genome", required=True, help="Path to the reference genome file")
    parser.add_argument("-q", "--left_trim", required=False, help= "Trim size left")
    parser.add_argument("-u", "--right_trim", required=False, help= "Trim size right")
    args = parser.parse_args()

    AB1VariantCall(args.input_folder, args.output_folder, args.reference_genome, args.left_trim, args.right_trim)


if __name__ == "__main__":
    main()