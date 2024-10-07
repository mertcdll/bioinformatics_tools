import vcf
import pandas as pd
import json
import os

class decipherClinvarVCF:
    
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

        self.all_records = []

        for record in self.vcf_file:
            self.all_records.append(record) 
        
        for record in self.all_records:

            variant_type = record.INFO.get("CLNVC", None)
            chrom = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = str(record.ALT[0])

            if variant_type == "single_nucleotide_variant":
                id = decipherSNVVariant(chrom, pos, ref, alt)
                try: 
                    record.INFO["ILYOME_ID"] = id
                except Exception as e:
                    pass
            
            elif variant_type not in ["single_nucleotide_variant", "Variant"]:
                id =  decipherCNVVariant(chrom, pos, ref, alt)
                try:
                    if id < ilyome_cutoff:   
                        record.INFO["ILYOME_ID"] = id
                    else:
                        record.INFO["ILYOME_ID"] = None
                except Exception as e:
                    pass

            else:
                try:
                    record.INFO["ILYOME_ID"] = None
                except Exception as e:
                    pass

        
        cnv_records = []

        for record in self.all_records:

            variant_type = record.INFO.get("CLNVC", None)
            ilyome_id = record.INFO["ILYOME_ID"]

            if variant_type not in ["single_nucleotide_variant", "Variant"] and ilyome_id is not None:

                cnv_records.append(record)


        ilyome_line = '##INFO=<ID=ILYOME_ID,Number=1,Type=Integer,Description="deciphered variant code for indels">'

        self.vcf_file._header_lines.append(ilyome_line)

        
        vcf_out = vcf.Writer(open(os.path.join(self.output_dir, "clinvar_w_ilidsnew.vcf"), "w"), self.vcf_file)

        for record in cnv_records:
            vcf_out.write_record(record)

        vcf_out.close()


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

        tab_cnv['ilyome_id'] = pd.to_numeric(tab_cnv['ilyome_id'], errors='coerce').astype(pd.Int64Dtype())   
        
        def stripper(item):
            if isinstance(item, str):
                if item.startswith("["):
                    item = item.strip("[]").replace("'","").replace(" _", "")
                else:
                    pass
            return item
        
        
        tab_cnv.to_csv(os.path.join(self.output_dir, "clinvar_w_ilids_tabular_all.txt"), sep="\t", index=False)

        tab_clin = pd.read_csv(os.path.join(self.output_dir, "clinvar_w_ilids_tabular_all.txt"), sep="\t")

        tab_clinvar  = tab_clin.applymap(stripper)

        tab_clinvar['ilyome_id'] = pd.to_numeric(tab_clinvar['ilyome_id'], errors='coerce').astype(pd.Int64Dtype())

        tab_clinvar.to_csv(os.path.join(self.output_dir, "clinvar_w_ilids_tabular.txt"), sep="\t", index=False)

        os.remove(os.path.join(self.output_dir, "clinvar_w_ilids_tabular_all.txt"))


vcf_path = "/home/genwork2/Mert/DecipherVCF/clinvar.vcf"
chr_code_cnv_path = "/home/genwork2/Mert/DecipherVCF/ints_cnv.json"
chr_code_snv_path = "/home/genwork2/Mert/DecipherVCF/ints_snv.json"
output_directory =  "/home/genwork2/Mert/DecipherVCFAD"


decipher = decipherClinvarVCF(vcf_file=vcf_path, chr_code_file_cnv=chr_code_cnv_path, chr_code_file_snv= chr_code_snv_path ,output_dir=output_directory)

decipher.add_ids_export_cnv_vcf()

decipher.vcf_tabular()

