import json
import pandas as pd

def getBaseINT(bases:str) -> str:

    ints = {
        "A":"1",
        "T":"2",
        "C":"3",
        "G":"4"
    }
    res = ""
    for base in bases:
        try:
            res += ints[base.upper()]
        except Exception as e:
            pass # handel the exception later
    return res

def getChromINT(chrom: str) -> str:
    with open('/home/genwork2/Mert/DecipherVCF/ints.json', 'r') as file:
        ints = json.load(file)
    ints = {key: str(value) for key, value in ints.items()}    

    try:
        return ints[chrom]
    except Exception as e:
        pass # handel the exception later

def decipherVariant(chrom:str, pos:int, ref:str, alt:str) -> int:
    if chrom.startswith("chr"):
        chrom =  chrom.strip("chr")
    chrom = getChromINT(chrom=chrom)
    ref= getBaseINT(ref)
    alt = getBaseINT(alt)
    try:
        return int(chrom + str(pos) + ref + alt)
    except Exception as e:
        pass # handel the exception later


#######################

class FastaVCF:
    def __init__(self, fasta, key_file):
        self.fasta = fasta
        self.key_file = self.load_key_file(key_file)
        self.readFasta()
        

    def load_key_file(self, key_file):
        df = pd.read_csv(key_file, sep ="\t")
        return df
    
    def get_chromosome_name(self, line):

        alias = line.strip(">").split()[0]

        match = self.key_file[self.key_file["alias"] == alias]

        if not match.empty:
            return match["chr_name"].iloc[0]


    def readFasta(self):

        Alternatives = {
            "A": ["C", "G", "T"],
            "C": ["A", "G", "T"],
            "G": ["A", "C", "T"],
            "T": ["A", "C", "G"]
        }
        NCtCommon = ["A", "C", "G", "T"]
        NCtids = {
            "R": ["G", "A"],
            "Y": ["C", "T"],
            "K": ["G", "T"],
            "M": ["A", "C"],
            "S": ["G", "C"],
            "W":["A", "T"],
        }
        output = None

        # output.write(self.vcfHeader)
        with open(self.fasta) as T:
            location = 0
            chro = None

            for line in T:
                if line.startswith(">"):
                    chro = self.get_chromosome_name(line)
                    
                    if chro:

                        outputLocations = f"Mert/DecipherVCF/All/Chr{chro}_GenomicVariants.vcf"
                        output = open(outputLocations, "w")
                        output.write("##fileformat=VCFv4.2\n")
                        output.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description='Approximate read depth (reads with MQ=255 or with bad mates are filtered)'>\n")
                        output.write("##INFO=<ID=ID,Number=A,Type=Integer,Description='Allele count in genotypes, for each ALT allele, in the same order as listed'>\n")
                        output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGENOMICVARIANTS\n")
                    
                    
                    location = 0
                              
                else:
                    for base in line.strip(""):
                        location += 1
                        if base in NCtCommon:
                                
                            for Alt in Alternatives[base]:
                                try:
                                    id = decipherVariant(chro, str(location), base, Alt)
                                    if id != 1:
                                        output.write(f"{'chr'+chro}\t{location}\t.\t{base}\t{Alt}\t1\t.\tID={id}\tDP\t1\n")
                                    else:
                                        print("ERROR: One variant was missed!")
                                        exit()
                                except KeyError:
                                    pass

                        elif base in NCtids:
                            for base in NCtids[base]:
                                for Alt in Alternatives[base]:
                                    try:
                                        id = decipherVariant(chro, str(location), base, Alt)
                                        if id != 1:
                                            output.write(f"{'chr'+chro}\t{location}\t.\t{base}\t{Alt}\t1\t.\tID={id}\tDP\t1\n")
                                        else:
                                            print("ERROR: One variant was missed!")
                                            exit()
                                    except KeyError:
                                        pass
        if output:
            output.close()     


if __name__ == "__main__":
    key_file = "/home/genwork2/Mert/DecipherVCF/ncbi_keys.txt"
    fasta = "/home/genwork2/Mert/DecipherVCF/GRCh38_latest_genomic.fna"
    t = FastaVCF(fasta, key_file)





