import pandas as pd
import json
import numpy as np

twist_16udi = pd.read_csv("/home/genwork2/Mert/samplesheets/16-UDI.csv", sep=",")

twist_16udi.columns = ["index", "product_type", "product_part_number", "product_name", "plate", "well","i5_index_f","i5_index_r", "i7_index"]


twist_16udi


twist_96udi = pd.read_csv("/home/genwork2/Mert/samplesheets/96-UDI.csv", sep=",")

twist_96udi.columns = ["index", "product_type", "product_part_number", "product_name", "plate", "well","i5_index_f","i5_index_r", "i7_index"]


twist_96udi


#######################

pi_list = []
for index, row in twist_16udi.iterrows():
    product_info = {
        'index': row['index'],
        'plate': row['plate'],
        'well': row['well'],
        'i5_index_f': row['i5_index_f'],
        'i5_index_r': row['i5_index_r'],
        'i7_index': row['i7_index']
    }
    
    pi_list.append(product_info)

pi_list

###############################


groups = twist_96udi.groupby("product_name")

listofplates = []

for group_name, group_df in groups:
    pi_list = []
    for index, row in group_df.iterrows():
        product_info = {
            'index': row['index'],
            'plate': row['plate'],
            'well': row['well'],
            'i5_index_f': row['i5_index_f'],
            'i5_index_r': row['i5_index_r'],
            'i7_index': row['i7_index']
        }
    
        pi_list.append(product_info)

    product_type = group_df["product_type"].iloc[0]
    product_part_number = str(group_df["product_part_number"].iloc[0])
    product_name = group_df["product_name"].iloc[0]

    new_dict = {
        'product_type' : product_type,
        'product_part_number' : product_part_number,
        'product_name' : product_name,
        'product_info' : pi_list
    }

    listofplates.append(new_dict)


for item in listofplates:
    print(item)


with open("/home/genwork2/Mert/samplesheets/96-UDI-TWIST.json", 'w') as json_file:
    json.dump(listofplates, json_file, indent=2)


################ 96 twist >> maybe you can add product part number inside the product info?

pi_list = []    

for index, row in twist_96udi.iterrows():
    product_info = {
    'index': row['index'],
    'plate': row['plate'],
    'well': row['well'],
    'i5_index_f': row['i5_index_f'],
    'i5_index_r': row['i5_index_r'],
    'i7_index': row['i7_index']
}

    pi_list.append(product_info)



from collections import OrderedDict
twist_96_dict = OrderedDict({
    'product_type': twist_96udi["product_type"].iloc[0],
    'product_name': twist_96udi["product_name"][0].split(" Plate")[0],
    'num_plates': 4,
    'rows_per_plate': 8,
    'columns_per_plate': 12, 
    'total_indexes': 384,
    'product_info': pi_list
})



with open("/home/genwork2/Mert/samplesheets/96-UDI-TWIST_merged.json", 'w') as json_file:
    json.dump(twist_96_dict, json_file, indent=2)


####### 16 twist >> maybe you can add product part number inside the product info?


pi_list = []    

for index, row in twist_16udi.iterrows():
    product_info = {
    'index': row['index'],
    'product_part_number': row["product_part_number"],
    'plate': row['plate'],
    'well': row['well'],
    'i5_index_f': row['i5_index_f'],
    'i5_index_r': row['i5_index_r'],
    'i7_index': row['i7_index']
}

    pi_list.append(product_info)


pi_list

from collections import OrderedDict
twist_16_dict = OrderedDict({
    'product_type': twist_16udi["product_type"].iloc[0],
    'product_name': twist_16udi["product_name"].iloc[0],
    'num_plates': 1,
    'rows_per_plate': 8,
    'columns_per_plate': 2, 
    'total_indexes': 16,
    'product_info': pi_list
})




with open("/home/genwork2/Mert/samplesheets/16-UDI-TWIST_merged.json", 'w') as json_file:
    json.dump(twist_16_dict, json_file, indent=2)


twist_96_dict
twist_16_dict

###### Agilent htsli


def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement_sequence = ''.join(complement_dict[base] for base in reversed(dna_sequence))
    return reverse_complement_sequence




agilent_htsli_p515 = pd.read_csv("/home/genwork2/Mert/samplesheets/agilenthtsli_p5_1-5.txt", sep="\s+")

agilent_htsli_p51 = pd.read_csv("/home/genwork2/Mert/samplesheets/agilenthtsli_p5_1.txt", sep="\s+")

agilent_htsli_p7 = pd.read_csv("/home/genwork2/Mert/samplesheets/agilenthtsli_p7.txt", sep="\s+")


agilent_htsli_v1 = pd.merge(agilent_htsli_p51, agilent_htsli_p7, on=["Index", "Well"])

agilent_htsli_v1.columns = ["Index","Well","p5","p7"]

agilent_htsli_v1_5 = pd.merge(agilent_htsli_p515, agilent_htsli_p7, on=["Index", "Well"])

agilent_htsli_v1_5.columns = ["Index","Well","p5","p7"]


agilent_htsli_v1
agilent_htsli_v1_5


#### v1

pi_list = []    

for index, row in agilent_htsli_v1.iterrows():
    product_info = {
    'index': row["Index"],
    'product_part_number': None,
    'plate': "Green",
    'well': row["Well"],
    'i5_index_f': row["p5"],
    'i5_index_r': reverse_complement(row['p5']),
    'i7_index': row["p7"]
}

    pi_list.append(product_info)




from collections import OrderedDict
agilent_htsli_v1_dict = OrderedDict({
    'product_type': "XT HS LIDI v1",
    'product_name': "Agilent SureSelect XT HS Low Input Dual Indexes for NovaSeq (v 1.0 chemistry)",
    'num_plates': 1,
    'rows_per_plate': 8,
    'columns_per_plate': 12, 
    'total_indexes': 96,
    'product_info': pi_list
})





with open("/home/genwork2/Mert/samplesheets/AgilentXTHSLIDI_v1.json", 'w') as json_file:
    json.dump(agilent_htsli_v1_dict, json_file, indent=2)



#### v1.5


pi_list = []    

for index, row in agilent_htsli_v1_5.iterrows():
    product_info = {
    'index': row["Index"],
    'product_part_number': None,
    'plate': "Green",
    'well': row["Well"],
    'i5_index_f': row["p5"],
    'i5_index_r': reverse_complement(row['p5']),
    'i7_index': row["p7"]
}

    pi_list.append(product_info)


from collections import OrderedDict
agilent_htsli_v15_dict = OrderedDict({
    'product_type': "XT HS LIDI v1.5",
    'product_name': "Agilent SureSelect XT HS Low Input Dual Indexes for NovaSeq (v 1.5 chemistry)",
    'num_plates': 1,
    'rows_per_plate': 8,
    'columns_per_plate': 12, 
    'total_indexes': 96,
    'product_info': pi_list
})



with open("/home/genwork2/Mert/samplesheets/AgilentXTHSLIDI_v1_5.json", 'w') as json_file:
    json.dump(agilent_htsli_v15_dict, json_file, indent=2)


#### Agilent hts2


agilent_hts2 = pd.read_csv("/home/genwork2/Mert/samplesheets/agilent_hts2.txt", sep="\s+")


agilent_hts2
plate_list = ['Orange'] * 96 + ['Blue'] * 96 + ['Green'] * 96 + ['Red'] * 96

plate_list
agilent_hts2["Plate"] = plate_list


pi_list = []    


for index, row in agilent_hts2.iterrows():
    product_info = {
    'index': row["Index"],
    'product_part_number': None,
    'plate': row["Plate"],
    'well': row["Well"],
    'i5_index_f': row["P5_Index_Forward"],
    'i5_index_r': row["P5_Index_Reverse_Complement"],
    'i7_index': row["P7_Index_Forward"]
}

    pi_list.append(product_info)


from collections import OrderedDict
agilent_hts2_dict = OrderedDict({
    'product_type': "XT HS2",
    'product_name': "Agilent SureSelect XT HS2 Dual Indexes",
    'num_plates': 4,
    'rows_per_plate': 8,
    'columns_per_plate': 12, 
    'total_indexes': 384,
    'product_info': pi_list
})


with open("/home/genwork2/Mert/samplesheets/AgilentXTHS2.json", 'w') as json_file:
    json.dump(agilent_hts2_dict, json_file, indent=2)



 #### nextera

nextera = pd.read_csv("/home/genwork2/Mert/samplesheets/nextera.csv", sep=",")

nextera

pi_list = []

for index, row in nextera.iterrows():
    product_info = {
    'index': row["Sample_ID"],
    'product_part_number': None,
    'plate': row["Index_Plate"],
    'well': row["Index_Plate_Well"],
    'i5_index_f': row["index2"],
    'i5_index_r': reverse_complement(row['index2']),
    'i7_index': row["index"]
}

    pi_list.append(product_info)


from collections import OrderedDict
nextera_dict = OrderedDict({
    'product_type': "Nextera UDI-384",
    'product_name': "IDT-ILMN Nextera DNA UD Indexes (384)",
    'num_plates': 4,
    'rows_per_plate': 8,
    'columns_per_plate': 12, 
    'total_indexes': 384,
    'product_info': pi_list
})


with open("/home/genwork2/Mert/samplesheets/nextera.json", 'w') as json_file:
    json.dump(nextera_dict, json_file, indent=2)


###### Compile


import os

file_list = os.listdir("/home/genwork2/Mert/samplesheets")

all_kits = []

for file in file_list:
    if file.endswith(".json"):
        with open(os.path.join("/home/genwork2/Mert/samplesheets", file), 'r') as f:
            indexes = json.load(f)

        all_kits.append(indexes)



all_kits



with open("/home/genwork2/Mert/samplesheets/all_kits.json", 'w') as json_file:
    json.dump(all_kits, json_file, indent=2)