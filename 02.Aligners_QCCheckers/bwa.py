import os

class Mapper:
    def _init_(self, **kwargs):
        
        for kwarg in kwargs.values():
            for key, value in kwarg.items():
                setattr(self, key, value)
        
        if not os.path.isdir(self.SampleOutputDir):
            os.mkdir(self.SampleOutputDir)

    def run(self):
        print(f"Aligning {self.rgsm} FASTQ files to the genome...")
        getbwa = f'bwa mem -M -t {self.max_workers} -R "@RG\\tID:{self.rgsm}\\tSM:{self.rgsm}\\tLB:Mylib\\tPU:Illumina" {self.Reference} {self.R1} {self.R2} > {self.output_sam}'
        getdragen = f'dragen-os -r {self.HashTable} -1 {self.R1} -2 {self.R2} --RGSM {self.rgsm} --RGID {self.rgsm} --num-threads {self.max_workers} > {self.output_sam}'
        try:
            # os.system(getbwa)
            os.system(getdragen)
            print(f'Alignment for {self.R1} and {self.R2} finished successfully!')
        except Exception as e:
            print(e)
            print(f'ERROR: Couldnt mapp the {self.R1} and {self.R2} to the genome reference!')
            return 1
        else:
            cmd = f"gatk --java-options '-Xms{self.memory}g -Xmx{self.memory}g' MarkDuplicatesSpark -I {self.output_sam} -O {self.out_bam_dedup} --metrics-file {self.Makr_dup_metrics} --spark-master local[{self.max_workers}] "
            try:
                os.system(cmd)
            except Exception as e:
                print(e)
                return 1
            else:
                if os.path.isfile(self.Makr_dup_metrics):
                    os.system(f"rm {self.output_sam}")
                    return 0
                else:
                    return 1
                
                
