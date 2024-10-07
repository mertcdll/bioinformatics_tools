import re

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

def getChromINT(chrom:str) -> str:

    ints = {
        "1":"111",
        "2":"112",
        "3":"113",
        "4":"114",
        "5":"115",
        "6":"116",
        "7":"117",
        "8":"118",
        "9":"119",
        "10":"120",
        "11":"121",
        "12":"122",
        "13":"123",
        "14":"124",
        "15":"125",
        "16":"126",
        "17":"127",
        "18":"128",
        "19":"129",
        "20":"130",
        "21":"131",
        "22":"132",
        "X":"133",
        "Y":"134",
        "M":"135",
        "MT": "135"
    }

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

class FastaVCF:
    def __init__(self, fasta):
        self.fasta = fasta

        self.outputLocations = "Mert/DecipherVCF/GenomicVariants.tsv"
        self.outputLocations = "Mert/DecipherVCF/GenomicVariants.tsv"
        self.readFasta()

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
        output = open(self.outputLocations, "w")
        # output.write(self.vcfHeader)
        with open(self.fasta) as T:
            chro = None
            location = 0
            pattern = r'chromosome ([\w\d]+)'
            for line in T:
                if line.startswith(">NC_"):
                    match = re.search(pattern, line)
                    if match:
                        chro = match.group(1)
                    elif "mitochondrion" in line:
                        chro = "M"   
                    location = 0
                else:
                    for base in line.strip():
                        location += 1
                        if base in NCtCommon:
                            for Alt in Alternatives[base]:
                                try:
                                    id = decipherVariant(chro, str(location), base, Alt)
                                    if id != 1:
                                        output.write(f"{id}\t{base}:{Alt}\n")
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
                                            output.write(f"{id}\t{base}:{Alt}\n")
                                        else:
                                            print("ERROR: One variant was missed!")
                                            exit()
                                    except KeyError:
                                        pass

        output.close()

if __name__ == "__main__":
    x = decipherVariant('chr1', 12397401234, "A", "C")
    fasta = "/home/genwork2/Mert/DecipherVCF/GRCh38_latest_genomic.fna"
    t = FastaVCF(fasta)