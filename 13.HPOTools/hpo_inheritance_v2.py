import pandas as pd
import json
import os


class GeneInheritanceProcessor:
    """This class takes inheritance mode/code file and gene_to_phenotype file as inputs and creates two json files
    containing genes/inheritance modes and genes/inheritance mode codes"""

    def __init__(self, inheritance_file, genes_phenotype_file, geneid_name_keys, output_dir):
        self.inheritance_file = inheritance_file
        self.genes_phenotype_file = genes_phenotype_file
        self.geneid_name_keys = pd.read_csv(geneid_name_keys, sep="\t")
        self.inh_modes = pd.read_csv(inheritance_file, sep="\t")
        self.gene_phenotype = pd.read_csv(genes_phenotype_file, sep="\t")
        self.output_dir = output_dir

    def process(self):
        """This method processes gene_to_phenotype file to create two json files as described in class"""

        genes_inh = self.gene_phenotype.loc[self.gene_phenotype["hpo_name"].isin(self.inh_modes["modes"]),]
        
        genes_inh_unq = genes_inh[~genes_inh.duplicated(subset=["gene_symbol", "hpo_id"], keep="first")]
        
        genes_inh_unq.reset_index(inplace=True, drop=True)

        gene_hpo_dict = {}
        
        for index, row in genes_inh_unq.iterrows():
            gene_symbol = row["gene_symbol"]
            hpo_name = row["hpo_name"]

            if gene_symbol in gene_hpo_dict:
                gene_hpo_dict[gene_symbol].append(hpo_name)
            else:
                gene_hpo_dict[gene_symbol] = [hpo_name]

        gene_code_dict = {}

        for gene, inheritance_modes in gene_hpo_dict.items():
            codes = []
            for mode in inheritance_modes:
                code = self.inh_modes[self.inh_modes['modes'] == mode]['codes'].values[0]
                codes.append(code)
            gene_code_dict[gene] = codes

        for gene, codes in gene_code_dict.items():
            gene_code_dict[gene] = [int(code) for code in codes]

        with open(os.path.join(self.output_dir,"gene_inheritance_code.json"), "w") as json_file:
            json.dump(gene_code_dict, json_file)

        with open(os.path.join(self.output_dir,"gene_inheritance_mode.json"), "w") as file:
            json.dump(gene_hpo_dict, file)



inheritance_mode_file = "/home/genwork2/Mert/hpo_inheritance/inheritance_modes.txt"
genes_to_phenotype_file = "/home/genwork2/Mert/hpo_inheritance/genes_to_phenotype.txt"
output_directory = "/home/genwork2/Mert/hpo_inheritance"


processor = GeneInheritanceProcessor(inheritance_mode_file, genes_to_phenotype_file, output_directory)
processor.process()
