import vcf
import pandas as pd
import json
import os

class decipherClinvarVCF:
    
    def __init__(self, vcf_file, chr_code_file, output_dir):
        self.vcf_file = self.importVCF(vcf_file)
        self.chr_code_file = self.importchrcodes(chr_code_file)
        self.output_dir = output_dir
        self.filtered_cnv_records = None

    def importchrcodes(self, chr_code_file):
        with open(chr_code_file, 'r') as json_file:
            ints = json.load(json_file)

        ints["MT"] = "3135"

        return ints
    
    def importVCF(self, vcf_file):
        clinvar_vcf = vcf.Reader(open(vcf_file, "r"))
        return clinvar_vcf

    
    def add_ids_export_cnv_vcf(self):

            
        def getChromINT(chrom:str) -> str:
            try:
                return self.chr_code_file[chrom]
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

        def decipherVariant(chrom:str, pos:int, ref:str, alt:str) -> int:
            if chrom.startswith("chr"):
                chrom =  chrom.strip("chr")
            chrom = getChromINT(chrom)
            refid= getBaseSum(ref)
            altid = getBaseSum(alt)
            try:
                return int(chrom + str(pos) + str(len(ref)) + str(refid) + str(len(alt)) + str(altid))
            except Exception as e:
                pass 

        cnv_records = []

        for record in self.vcf_file:

            variant_type = record.INFO.get("CLNVC", None)

            if variant_type not in ["single_nucleotide_variant", "Variant"]:

                cnv_records.append(record)


        ids = []        


        for record in cnv_records:

            chrom  = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = str(record.ALT[0])

            id = decipherVariant(chrom, pos, ref, alt)

            ids.append(id)


        indices = []
        for i, item in enumerate(ids):
            if item < 9223372036854775807:
                indices.append(i)

        
        filtered_ids = [ids[i] for i in indices]

        self.filtered_cnv_records  = [cnv_records[i] for i in indices]


        for record, id in zip(self.filtered_cnv_records, filtered_ids):
            record.INFO["ILYOMEID"] = id


        ilyome_line = '##INFO=<ID=ILYOMEID,Number=1,Type=Integer,Description="deciphered variant code for indels">'

        self.vcf_file._header_lines.append(ilyome_line)

        
        vcf_out = vcf.Writer(open(os.path.join(self.output_dir, "clinvar_w_ilids.vcf"), "w"), self.vcf_file)

        for record in self.filtered_cnv_records:
            vcf_out.write_record(record)

        vcf_out.close()

    

    def vcf_tabular(self):


        info_df = pd.DataFrame([record.INFO for record in self.filtered_cnv_records])

        met_df = pd.DataFrame({
            "CHROM": [record.CHROM for record in self.filtered_cnv_records],
            "POS": [record.POS for record in self.filtered_cnv_records],
            "ID": [record.ID for record in self.filtered_cnv_records],
            "REF": [record.REF for record in self.filtered_cnv_records],
            "ALT": [record.ALT for record in self.filtered_cnv_records],
            "QUAL": ".",
            "FILTER": "."
        })

        ### we can insert the 


        tab_cnv = pd.concat([met_df, info_df], axis = 1)

        tab_cnv.columns = tab_cnv.columns.str.lower()

        tab_cnv.to_csv(os.path.join(self.output_dir, "clinvar_w_ilids_tabular.txt"), sep="\t", index=False)





vcf_path = "/home/genwork2/Mert/DecipherVCF/clinvar.vcf"
chr_code_path = "/home/genwork2/Mert/DecipherVCF/ints.json"
output_directory =  "/home/genwork2/Mert/DecipherVCFAD"


decipher = decipherClinvarVCF(vcf_file=vcf_path, chr_code_file=chr_code_path, output_dir=output_directory)

decipher.add_ids_export_cnv_vcf()
decipher.vcf_tabular()