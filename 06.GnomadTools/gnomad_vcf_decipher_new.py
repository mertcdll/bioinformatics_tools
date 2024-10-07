import vcf
import json
import os
import gzip
import copy

class Decipher:
    
    """This class deciphers variants according to their variant type.
    If it is a SNV (Single Nucleotide Variant), it assigns numbers to each base 
    (A: 1, T: 2, C: 3, G: 4). If it is a indel or any other variant type comprising multiple 
    base changes, it first provides the number of bases in the alternate (alt) and reference (ref) sequences, 
    and then provides the sum of the values assigned to the bases. 
    Each resulting ID includes a chromosome code, variant position in chromosome, and a base change code.
    """

    def __init__(self, chr_code_file_indel, chr_code_file_snv):
        
        """This method initializes the Decipher class"""

        self.chr_code_file_snv = self.importchrcodessnv(chr_code_file_snv)
        self.chr_code_file_indel = self.importchrcodesindel(chr_code_file_indel)

    def importchrcodesindel(self, chr_code_file):

        """This method imports chromosome codes for indels and structural variations."""

        with open(chr_code_file, 'r') as json_file:
            ints = json.load(json_file)

        ints["MT"] = "3135"

        return ints
    
   
    def importchrcodessnv(self, chr_code_file):

        """This method imports chromosome codes SNVs."""
        
        with open(chr_code_file, 'r') as json_file:
            ints = json.load(json_file)
        
        ints["MT"] = "1135"
        
        return ints

    
    def decipherSNVVariant(self, chrom, pos, ref, alt):

        """This method deciphers SNVs. It first gets chromosome codes designated for each chromosome, then 
           assigns numbers to each base (A: 1, T: 2, C: 3, G: 4). Each final ID includes, a chromosome chromosome code,
           variant position in chromosome and base code.""" 
        
        def getChromINTSNV(chrom):
            try:
                return self.chr_code_file_snv[chrom]
            except Exception as e:
                pass

        def getBaseINT(bases):
            ints = {
                "A": "1",
                "T": "2",
                "C": "3",
                "G": "4"
            }
            res = ""
            for base in bases:
                try:
                    res += ints[base.upper()]
                except Exception as e:
                    pass  
            return res

        chrom = getChromINTSNV(chrom=chrom)
        refid = getBaseINT(ref)
        altid = getBaseINT(alt)
        try:
            return int(chrom + str(pos) + str(refid) + str(altid))
        except Exception as e:
            pass
    
     
    def decipherINDELVariant(self, chrom, pos, ref, alt):

        """This method deciphers indels or other structural variants if exists. It first gets chromosome code designated for each chromosome, then provides the number
        of bases in reference (ref) sequences and the sum of values assigned to the bases. It does the same operation for alternate (alt) sequences. Each final ID includes a chromosome code,
        variant position in chromosome and base code"""

        def getChromINTINDEL(chrom):
            try:
                return self.chr_code_file_indel[chrom]
            except Exception as e:
                pass

        def getBaseSum(bases):
            ints = {
                "A": 1,
                "T": 2,
                "C": 3,
                "G": 4
            }
            res = 0
            for base in bases:
                try:
                    res += ints[base.upper()]
                except Exception as e:
                    pass  
            return res


        chrom = getChromINTINDEL(chrom)
        refid = getBaseSum(ref)
        altid = getBaseSum(alt)
        try:
            return int(chrom + str(pos) + str(len(ref)) + str(refid) + str(len(alt)) + str(altid))
        except Exception as e:
            pass


class decipherGnomadVCF:
    
    """This class is designed to process GnomAD VCF files and add deciphered variant IDs to them. 
    It utilizes the Decipher class to generate unique IDs for different types of variants, including SNVs and indels.
    The processed VCF file is exported in both tabular and compressed VCF formats, with additional information included.
    Processed VCF file only includes indels and structural variants if they exist. Tabular format, on the other hand, includes
    all variant data but not all information found in the info column. Only allele counts/number/frequency values for each race and histogram values 
    can be extracted from info column"""
    
    def __init__(self, vcf_file, chr_code_file_indel, chr_code_file_snv, output_dir):
        """This method initializes the decipherGnomadVCF class."""
        self.vcf_file = self.importVCF(vcf_file)
        self.output_dir = output_dir
        self.decipher = Decipher(chr_code_file_indel, chr_code_file_snv)
        self.vcf_name = vcf_file.split("/")[-1]


    def importVCF(self, vcf_file):
        """This method imports the GnomAD VCF file using PyVCF and returns a VCF Reader instance."""     
        try:
            vcf_reader = vcf.Reader(filename=vcf_file)
            return vcf_reader
        except UnicodeDecodeError:
            vcf_reader = vcf.Reader(open(vcf_file, "r", encoding="utf-8"))
            return vcf_reader

    
    def add_ids_export_vcf_txt(self):
        
        """This method processes VCF records, deciphers variant IDs, and exports the modified VCF file.
        It iterates through VCF records, determines the variant type (SNV, indel, or mixed), and
        adds deciphered 'ilyome_id' INFO fields to each record. It then exports the processed VCF data
        in both tabular (TXT) and compressed VCF (VCF.gz) formats."""

        ilyome_cutoff = 9223372036854775807
        ilyome_line = '##INFO=<ID=ilyome_id,Number=1,Type=Integer,Description="deciphered variant code for indels">'
        new_header_lines = copy.deepcopy(self.vcf_file._header_lines)

        new_header_lines.append(ilyome_line)
        print(new_header_lines)
        self.vcf_file._header_lines = new_header_lines
        
        infos = self.vcf_file.infos
        allelekeys = []

        exclude_substrings = ["seu", "bgr", "swe", "nwe", "jpn", "kor", "oea", "est", "onf"]

        for item in infos.keys():
            if not any(substring in item for substring in exclude_substrings):
                if item.startswith('AC') or item.startswith('AF') or item.startswith('AN'):
                    allelekeys.append(item)
                elif 'hist' in item:
                    allelekeys.append(item)
                elif 'nhomalt' in item: 
                    allelekeys.append(item)
                elif 'CTT_p_value' in item:
                    allelekeys.append(item)
                elif 'CMH'in item:
                    allelekeys.append(item)
                elif 'stat_union' in item:
                    allelekeys.append(item)
                elif 'not_called' in item:
                    allelekeys.append(item)

        
        allelekeys  = allelekeys + ["ilyome_id"]

        output_file_path = os.path.join(self.output_dir, f'{os.path.splitext(self.vcf_name)[0].replace(".vcf", "")}.dcp.vcf.gz')
        print(self.output_dir)
        txt_output_file_path = os.path.join(self.output_dir, f'{self.vcf_name.replace(".vcf.gz", "")}.dcp.txt')
        txt_output_file = open(txt_output_file_path, "w")
        print(txt_output_file_path)

        compressed_file = gzip.open(output_file_path, "wt")
        vcf_out = vcf.Writer(compressed_file, self.vcf_file)

        header = ["chrom", "pos", "ref", "alt"] + allelekeys

        header = [item.lower() for item in header]
        
        txt_output_file.write("\t".join(header) + "\n")

        for record in self.vcf_file:

            chrom = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = str(record.ALT[0])
            
            if ref and alt and pos:
                if len(ref) + len(alt) == 2:
                    id = self.decipher.decipherSNVVariant(chrom, pos, ref, alt)
                    try: 
                        record.INFO["ilyome_id"] = id
                        vcf_out.write_record(record)
                    except Exception as e:
                        pass
                else:
                    id =  self.decipher.decipherINDELVariant(chrom, pos, ref, alt)
                    try:
                        if id < ilyome_cutoff:   
                            record.INFO["ilyome_id"] = id
                            vcf_out.write_record(record)
                        else:
                            record.INFO["ilyome_id"] = None
                            vcf_out.write_record(record)
                    except Exception as e:
                        pass
            else:
                record.INFO["ilyome_id"] = None
                vcf_out.write_record(record)
  
            ac_list = record.INFO.get("AC", "")

            max_index = None 
            
            if len(ac_list) >= 2:
                 max_index = ac_list.index(max(ac_list))     
            
            info_values = []
            
            for key in allelekeys:
                info = record.INFO.get(key, "")
                
                if isinstance(info, float):

                    info_values.append(round(info, 6))
                
                elif isinstance(info, int):
                    
                    info_values.append(info)
                
                elif isinstance(info, list):

                    if max_index:
                        if all(isinstance(x, float) for x in info):
                            max_value = max(info)
                            info_values.append(round(max_value,6))

                        elif all(isinstance(x, int) for x in info):
                            max_value = max(info)
                            info_values.append(max_value)
                        
                        elif all(isinstance(x, str) for x in info):
                            info_values.append(info[max_index])
                    else:
                        if all(isinstance(x, float) for x in info):
                            info_values.append(round(info[0],6))
                        elif all(isinstance(x, int) for x in info):
                            info_values.append(info[0])
                        elif all(isinstance(x, str) for x in info):
                            info_values.append(info[0])                        

                else:
                    info_values.append(info)


            txt_output_file.write("\t".join([chrom, str(pos), ref, alt] + [str(value) for value in info_values]) + "\n")

        txt_output_file.close()


vcf_path = "/home/genwork2/Mert/gnomad4/exomesftrm/gnomad.exomes.v4.0.sites.chr1.dcp.ftrm.bg.vcf.gz"
chr_code_indel_path = "/home/genwork2/Mert/DecipherUCSCFasta/ints_indel.json"
chr_code_snv_path = "/home/genwork2/Mert/DecipherUCSCFasta/ints_snv.json"
output_directory =  "/home/genwork2/Mert/gnomad4joint"


decipher = decipherGnomadVCF(vcf_file=vcf_path, chr_code_file_indel=chr_code_indel_path, chr_code_file_snv= chr_code_snv_path, output_dir=output_directory)

decipher.add_ids_export_vcf_txt()
