import pandas as pd
import os
import json
import csv

class feature_counter:
    def __init__(self, file_name, feature_column, count_column):
        self.file_name = file_name
        self.feature_column = feature_column
        self.count_column = count_column
        self.df = self.read_file_explode()
        self.df2 = self.read_file()
    # Reads the dataset, if there is multiple consequences in one cell, splits them and adds new rows (other columns will be identical).
    def read_file(self):
        df = pd.read_csv(self.file_name, sep='\t', header=0)
        df.drop('Unnamed: 0', axis=1, inplace=True)

        return df

    def read_file_explode(self): 
        df = pd.read_csv(self.file_name, sep='\t', header=0)
        df.drop('Unnamed: 0', axis=1, inplace=True)
        df["Consequence"] = df["Consequence"].str.split("&")
        exploded_df = df.explode("Consequence")
        exploded_df.reset_index(drop=True)
        # Merge some expressions.
        exploded_df['Clinvar_CLNSIG'] = exploded_df['Clinvar_CLNSIG'].str.replace(r'^Benign/Likely_benign.*', 'Benign', regex=True)
        exploded_df['Clinvar_CLNSIG'] = exploded_df['Clinvar_CLNSIG'].str.replace(r'^Pathogenic/Likely_pathogenic.*', 'Pathogenic', regex=True)

        return exploded_df
    # This is a count calculator, counts values in any feature of grouped objects.
    def count_calculator(self):
        counts = self.df.groupby(self.feature_column)[self.count_column].value_counts().unstack(fill_value=0)
        return counts
    # Creates pivot tables of clinical significance versus consequence for each gene.
    def combination_counter(self):
        grouped_by_gene = self.df.groupby(self.feature_column)
        gene_dict = {}
        clinsig_list = self.df['Clinvar_CLNSIG'].value_counts().index.to_list()
        pathogenicity_list = [item for item in clinsig_list if 'Pathogenic' in item or 'Likely_pathogenic' in item]
        benign_list = [item for item in clinsig_list if 'Benign' in item or 'Likely_benign' in item]
        total_non_vus_list = pathogenicity_list + benign_list

        for gene, group in grouped_by_gene:
            sub_grouped = group.groupby(['Clinvar_CLNSIG','Consequence']).size().unstack(fill_value=0)
            gene_dict[gene] = sub_grouped
        
        # Adds several rows to pivot table.
        for df in gene_dict.values():
            df['Total'] = df.sum(axis=1)

            df.loc['Total'] = df.sum()
            # Takes matches betwwen the lists and pivot table's index.    
            valid_total_non_vus_list = df.index.intersection(total_non_vus_list)
            if valid_total_non_vus_list.any():
                df.loc['Total_non_vus'] = df.loc[valid_total_non_vus_list].sum()
            else:
                df.loc['Total_non_vus'] = 0
            # Takes matches betwwen the lists and pivot table's index. 
            valid_pathogenicity_list = df.index.intersection(pathogenicity_list)
            if valid_pathogenicity_list.any():
                df.loc['Total_pathogenic'] = df.loc[valid_pathogenicity_list].sum()
            else:
                df.loc['Total_pathogenic'] = 0

            df.loc['Pathogenicity_ratio'] = df.loc['Total_pathogenic'] / df.loc['Total_non_vus']
            df.fillna(0, inplace=True)

        return (gene_dict)
    # Creates bed file
    def create_clsig_conseq_bed_file(self):
        gene_dict = self.combination_counter()
        gdict_df = pd.DataFrame({'gene_name': list(gene_dict.keys()), 'info': [json.dumps(val.to_dict()) for val in gene_dict.values()]})
        
        bed_df = pd.DataFrame({
            'gene': gdict_df['gene_name'],
            'start': pd.Series('1', index=range(len(gdict_df))),
            'end': pd.Series('4', index=range(len(gdict_df))),
            'info': gdict_df['info']                               
        })

        
        bed_file_path = self.file_name.rsplit('.', 1)[0] + '.conseq_clinsig.bed'
        bed_df.to_csv(bed_file_path, sep="\t", index=False, header=False,  quoting=csv.QUOTE_NONE, quotechar="", escapechar="\\")

        sorted_bed_file_path = bed_file_path.rsplit('.', 1)[0] + '.sorted.bed'
        gzipped_bed_file_path = sorted_bed_file_path + '.gz'

        os.system(f"sort -k1,1 {bed_file_path} > {sorted_bed_file_path}")
        os.system(f"bgzip -c {sorted_bed_file_path} > {gzipped_bed_file_path}")
        os.system(f"tabix -p bed {gzipped_bed_file_path}")
    
    def create_gene_pat_bed_file(self):
        pathogenic_df = self.df2[self.df2['Clinvar_CLNSIG'].str.contains('Pathogenic|Likely_pathogenic', case=True, na=False)]
        pathogenic_df['Clinvar_CLNSIG'] = pathogenic_df['Clinvar_CLNSIG'].str.replace('/', '-')
        pathogenic_df.reset_index(drop=True, inplace=True)
        
        
        pathogenic_df['Pos:Ref:Alt'] = pathogenic_df['POS'].astype(str) + ':' + pathogenic_df['REF'] + ':' + pathogenic_df['ALT']
        
        sorted_pathogenic_df = pathogenic_df.sort_values(by='MAX_AF', ascending=False)
        
        sorted_pathogenic_df.dropna(subset='SYMBOL', inplace=True)

        sorted_pathogenic_df.reset_index(drop=True, inplace=True)

        symbol_dict_df = sorted_pathogenic_df.groupby('SYMBOL').agg({
            'Clinvar_CLNHGVS': 'first',
            'Clinvar_CLNSIG': 'first',
            'MAX_AF': 'first'
        }).reset_index()
        


        symbol_dict_df2 = sorted_pathogenic_df.groupby('SYMBOL').agg({
            'Pos:Ref:Alt': 'first',
            'Clinvar_CLNSIG': 'first',
            'MAX_AF': 'first'
        }).reset_index()

        
        dict_df = symbol_dict_df.drop('SYMBOL', axis=1)
        dict_df2 = symbol_dict_df2.drop('SYMBOL', axis=1)
        
        dict_df_json = dict_df.apply(lambda row: row.to_json(), axis=1)
        dict_df_json_df = pd.DataFrame({'json_string': dict_df_json})

        dict_df2_json = dict_df2.apply(lambda row: row.to_json(), axis=1)
        dict_df2_json_df = pd.DataFrame({'json_string': dict_df2_json})


        bed1 = pd.DataFrame({
            'SYMBOL': symbol_dict_df['SYMBOL'],
            'START': pd.Series('1', index=range(len(symbol_dict_df))),
            'END': pd.Series('4', index=range(len(symbol_dict_df))),
            'INFO': dict_df_json_df['json_string']
        })
        
        bed1_file_path = self.file_name.rsplit('.', 1)[0] + '.pat_af_v1.bed'
        bed1.to_csv(bed1_file_path, sep="\t", index=False, header=False,  quoting=csv.QUOTE_NONE, quotechar="", escapechar="\\")

        sorted_bed1_file_path = bed1_file_path.rsplit('.', 1)[0] + '.sorted.bed'
        gzipped_bed1_file_path = sorted_bed1_file_path + '.gz'

        os.system(f"sort -k1,1 {bed1_file_path} > {sorted_bed1_file_path}")
        os.system(f"bgzip -c {sorted_bed1_file_path} > {gzipped_bed1_file_path}")
        os.system(f"tabix -p bed {gzipped_bed1_file_path}")
        
        bed2 = pd.DataFrame({
            'SYMBOL': symbol_dict_df2['SYMBOL'],
            'START': pd.Series('1', index=range(len(symbol_dict_df2))),
            'END': pd.Series('4', index=range(len(symbol_dict_df2))),
            'INFO': dict_df2_json_df['json_string']
        })

        bed2_file_path = self.file_name.rsplit('.', 1)[0] + '.pat_af_v2.bed'
        bed2.to_csv(bed2_file_path, sep="\t", index=False, header=False,  quoting=csv.QUOTE_NONE, quotechar="", escapechar="\\")

        sorted_bed2_file_path = bed2_file_path.rsplit('.', 1)[0] + '.sorted.bed'
        gzipped_bed2_file_path = sorted_bed2_file_path + '.gz'

        os.system(f"sort -k1,1 {bed2_file_path} > {sorted_bed2_file_path}")
        os.system(f"bgzip -c {sorted_bed2_file_path} > {gzipped_bed2_file_path}")
        os.system(f"tabix -p bed {gzipped_bed2_file_path}")
    # Saves gene dictionaries as json file.
    
    def save_as_json(self, file_name):
        gene_dict = self.combination_counter()
        json_dict = {gene: df.to_dict() for gene, df in gene_dict.items()}

        with open(file_name, "w") as json_file:
            json.dump(json_dict, json_file, indent=1)


counter = feature_counter('/home/genwork2/clinvar_result.txt', 'SYMBOL', 'IMPACT')

exploded_df = counter.read_file_explode()
result = counter.count_calculator()
gdict = counter.combination_counter()
counter.create_clsig_conseq_bed_file()



df = counter.read_file()
counter.create_gene_pat_bed_file()