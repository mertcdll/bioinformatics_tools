import vcf
import gzip
import os
import requests
import copy

class GnomadGenomeVCFTrimmer:
    """This class is designed to process GnomAD VCF files and trim off the extra fields. make a smaller vcf for speed up the process."""
    
    def __init__(self, outputDir=os.getcwd()):

        lists = ['Y']
    
        for chromosome in lists:
            vcf_url = f"https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr{chromosome}.vcf.bgz"
            vcf_path = os.path.join(outputDir, f"gnomad.joint.v4.1.sites.chr{chromosome}.vcf.bgz")
            vcf_output_converted = os.path.join(outputDir,  f"gnomad.joint.v4.1.sites.chr{chromosome}.conv.vcf.gz")
            vcf_output_processed = os.path.join(outputDir,  f"gnomad.joint.v4.1.sites.chr{chromosome}.vcf.gz") 
            vcf_output_processed_ftrm = os.path.join(outputDir,  f"gnomad.joint.v4.1.sites.ftrm.chr{chromosome}.vcf.gz") 
            #Download the file
            print(f"For chromosome {chromosome}, processing started")
            self.download_and_convert_vcf_to_txt(vcf_url, vcf_path)
            # convert bgz > gz 
            print(f"Converting BGZ to GZ {vcf_path} ... ")
            os.system(f"zcat {vcf_path} | gzip > {vcf_output_converted}")
            #Read vcf file
            self.vcf_file = vcf.Reader(filename=vcf_output_converted)
            # Perper the writer
            compressed_file = gzip.open(vcf_output_processed, "wt")
            ftrm_compressed_file = gzip.open(vcf_output_processed_ftrm, "wt")
            self.vcf_file_writer = vcf.Writer(compressed_file, self.vcf_file)
            self.vcf_trimmed_writer = vcf.Writer(ftrm_compressed_file, self.vcf_file)
            
            print(f"Working on {vcf_path} ...")
            self._process()
            
            #cleaning Area
            os.remove(vcf_path)
            os.remove(vcf_output_converted)
            #os.system(f"zcat {vcf_output_processed} | bgzip -c > {vcf_path}")
            
            #os.system(f"tabix -p vcf {vcf_path}")

            #os.remove(vcf_output_processed)
            
            print(f"Finished successfuly: {vcf_path}")

    def download_and_convert_vcf_to_txt(self, url, output_txt_file):
        try:
            print("Downloading the file ...")
            response = requests.get(url, stream=True)
            response.raise_for_status()

            with open(output_txt_file, 'wb') as output_file:
                for chunk in response.iter_content(chunk_size=81922):
                    output_file.write(chunk)
        except Exception as e:
            print(e)
    
    def _process(self):
        """Reads a vcf and removes unnecessary fields and write it again to new vcf"""
        
        for record in self.vcf_file:

            fields_to_keep_short =  {}
            fields_to_keep_long = {}

            for key in record.INFO:

                if any(sub in key for sub in ('AC', 'AF', 'AN', 'hist', 'nhomalt', 'stat_union', 'CMH', 'CTT_p_value', 'not_called')):
                    fields_to_keep_long[key] = record.INFO[key]
                if key in ("AC_joint", "AF_joint", "AN_joint", "nhomalt_joint", "AC_exomes", "AF_exomes", "AN_exomes", "nhomalt_exomes", "AC_genomes", "AF_genomes", "AN_genomes", "nhomalt_genomes", 'CTT_p_value', 'CMH_p_value', 'stat_union_p_value','not_called_in_exomes','not_called_in_genomes'):
                    fields_to_keep_short[key] = record.INFO[key]
        
            print(fields_to_keep_long)
            print(fields_to_keep_short)
        
            new_record_long = copy.deepcopy(record)
            new_record_long.INFO.clear()
            new_record_long.INFO.update(fields_to_keep_long)

            new_record_short = copy.deepcopy(record)
            new_record_short.INFO.clear()
            new_record_short.INFO.update(fields_to_keep_short)

            self.vcf_file_writer.write_record(new_record_long)

            self.vcf_trimmed_writer.write_record(new_record_short)

        
        
if __name__ == "__main__":
    GnomadGenomeVCFTrimmer()



