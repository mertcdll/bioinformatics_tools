import vcf
import pandas as pd
import json
import os
import gzip


class decipherGnomadVCF:
    
    def __init__(self, vcf_file, chr_code_file_cnv, chr_code_file_snv, output_dir):
        self.vcf_file = self.importVCF(vcf_file)
        self.chr_code_file_snv = self.importchrcodessnv(chr_code_file_snv)
        self.chr_code_file_cnv = self.importchrcodescnv(chr_code_file_cnv)
        self.output_dir = output_dir

    def importchrcodescnv(self, chr_code_file):
        with open(chr_code_file, 'r') as json_file:
            ints = json.load(json_file)

        ints["MT"] = "3135"

        return ints
    
    def importchrcodessnv(self, chr_code_file):
        with open(chr_code_file, 'r') as json_file:
            ints = json.load(json_file)
        
        ints["MT"] = "1135"
        
        return ints
    
    def importVCF(self, vcf_file):
        clinvar_vcf = vcf.Reader(open(vcf_file, "r"))
        return clinvar_vcf

    
    def add_ids_export_cnv_vcf(self):

            
        def getChromINTCNV(chrom:str) -> str:
            try:
                return self.chr_code_file_cnv[chrom]
            except Exception as e:
                pass
        
        def getChromINTSNV(chrom:str) -> str:
            try:
                return self.chr_code_file_snv[chrom]
            except Exception as e:
                pass

        def getBaseSum(bases:str) -> str:

            ints = {
                "A":1,
                "T":2,
                "C":3,
                "G":4
            }
            res = 0
            for base in bases:
                try:
                    res += ints[base.upper()]
                except Exception as e:
                    pass # handel the exception later
            return res
        
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

        def decipherCNVVariant(chrom:str, pos:int, ref:str, alt:str) -> int:
            if chrom.startswith("chr"):
                chrom =  chrom.strip("chr")
            chrom = getChromINTCNV(chrom)
            refid= getBaseSum(ref)
            altid = getBaseSum(alt)
            try:
                return int(chrom + str(pos) + str(len(ref)) + str(refid) + str(len(alt)) + str(altid))
            except Exception as e:
                pass 

        def decipherSNVVariant(chrom: str, pos:int, ref:str, alt:str) -> int:
            if chrom.startswith("chr"):
                chrom =  chrom.strip("chr")
            chrom = getChromINTSNV(chrom=chrom)
            refid= getBaseINT(ref)
            altid = getBaseINT(alt)
            try:
                return int(chrom + str(pos) + str(refid) + str(altid))
            except Exception as e:
                pass

        ilyome_cutoff = 9223372036854775807
        ilyome_line = '##INFO=<ID=ILYOME_ID,Number=1,Type=Integer,Description="deciphered variant code for indels">'
        self.vcf_file._header_lines.append(ilyome_line)
        

        output_file_path = os.path.join(self.output_dir, "gnomad_cnvvcf_wid.vcf.gz")
        compressed_file = gzip.open(output_file_path, "wt")
        vcf_out = vcf.Writer(compressed_file, self.vcf_file)

        for record in self.vcf_file:

            variant_type = record.INFO.get("variant_type", None)
            chrom = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = str(record.ALT[0])

            if variant_type in ["snv", "multi-snv"]:
                id = decipherSNVVariant(chrom, pos, ref, alt)
                try: 
                    record.INFO["ilyome_id"] = id
                except Exception as e:
                    pass
            
            elif variant_type in ["indel", "multi-indel"]:
                id =  decipherCNVVariant(chrom, pos, ref, alt)
                try:
                    if id < ilyome_cutoff:   
                        record.INFO["ilyome_id"] = id
                        vcf_out.write_record(record)
                    else:
                        record.INFO["ilyome_id"] = None
                except Exception as e:
                    pass

                vcf_out.write_record(record)


            elif variant_type == "mixed":
                if len(ref) + len(alt) == 2:
                    id = decipherSNVVariant(chrom, pos, ref, alt)
                    try:
                        record.INFO["ilyome_id"] = id
                    except Exception as e:
                        pass
                
                else:
                    id = decipherCNVVariant(chrom, pos, ref, alt)
                    try:
                        if id < ilyome_cutoff:
                            record.INFO["ilyome_id"] = id
                            vcf_out.write_record(record)
                        else:
                            record.INFO["ilyome_id"] = None
                    except Exception as e:
                        pass
        

    def vcf_tabular(self):


        info_df = pd.DataFrame([record.INFO for record in self.all_records])

        met_df = pd.DataFrame({
            "CHROM": [record.CHROM for record in self.all_records],
            "POS": [record.POS for record in self.all_records],
            "ID": [record.ID for record in self.all_records],
            "REF": [record.REF for record in self.all_records],
            "ALT": [record.ALT for record in self.all_records],
            "QUAL": ".",
            "FILTER": "."
        })



        tab_cnv = pd.concat([met_df, info_df], axis = 1)

        tab_cnv.columns = tab_cnv.columns.str.lower()
        
        
        tab_cnv = tab_cnv.astype({col: 'object' for col in tab_cnv.columns})
        

        tab_cnv['ilyome_id'] = pd.to_numeric(tab_cnv['ilyome_id'], errors='coerce').astype(pd.Int64Dtype())   

        def stripper(item):
            if isinstance(item, str):
                if item.startswith("["):
                    item = item.strip("[]").replace("'","").replace(" _", "")
                else:
                    return item
            return item
        

        tab_clinvar  = tab_cnv.applymap(stripper)

        tab_clinvar.to_csv(os.path.join(self.output_dir, "clinvar_w_ilids_tabularnew.txt"), sep="\t", index=False)



vcf_path = "/home/genwork2/Mert/DecipherVCFAD/gnomad.exomes.r2.1.1.sites.vcf"
chr_code_cnv_path = "/home/genwork2/Mert/DecipherVCF/ints_cnv.json"
chr_code_snv_path = "/home/genwork2/Mert/DecipherVCF/ints_snv.json"
output_directory =  "/home/genwork2/Mert/DecipherVCFAD/gnomadvcfs"


decipher = decipherGnomadVCF(vcf_file=vcf_path, chr_code_file_cnv=chr_code_cnv_path, chr_code_file_snv= chr_code_snv_path ,output_dir=output_directory)

decipher.add_ids_export_cnv_vcf()

#decipher.vcf_tabular()