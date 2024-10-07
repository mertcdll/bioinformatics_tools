import pandas as pd 
import numpy as np
import json
from natsort import natsorted

fastaheaders = pd.read_csv("/home/genwork2/Mert/DecipherUCSCFasta/ucscfastaheaders.txt", header=None)

fastaheaders[0] = fastaheaders[0].str.replace(">", "")

fastaheaders = fastaheaders.rename(columns={fastaheaders.columns[0]: "chromosomes"})

chromosomes = fastaheaders["chromosomes"].to_list()

regular_chr = [item for item in chromosomes if len(item) <= 5]

contigs = [item for item in chromosomes if len(item) > 5]

sorted_rc = natsorted(regular_chr)

sorted_contigs = natsorted(contigs)

sorted_rc[22:] = ["chrX", "chrY", "chrM"]

sorted_rc

all_chromosomes = sorted_rc + sorted_contigs

all_chromosomes

len(all_chromosomes)

chr_ids = pd.DataFrame(
  {'chr_name': all_chromosomes,
   'code': np.arange(1111, 1111+len(chromosomes))}
)

chr_ids.to_csv("/home/genwork2/Mert/DecipherUCSCFasta/chr_ids.txt", sep="\t", index=False)
chr_ids

chr_ids['code'] = chr_ids['code'].astype(str)

ints = {key: value for key, value in zip(chr_ids["chr_name"], chr_ids["code"])}

json_ints = json.dumps(ints)

with open("/home/genwork2/Mert/DecipherUCSCFasta/ints_indel.json", "w") as json_file:
    json_file.write(json_ints)



for key in ints:
    ints[key] += 2000

ints

for key,value in ints.items():
    ints[key] = str(value)