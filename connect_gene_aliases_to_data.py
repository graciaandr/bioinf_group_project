from enum import unique
import pandas as pd
import allel as sc
import numpy as np
import matplotlib.pyplot as plt 

SNP_data = pd.read_csv("gene_names_filtered_for_PK.csv", sep = ",")
gene_alias_data = pd.read_csv("gene_name_alias.csv", sep = ",")

### Example for lncRNA ==> ask about this in tutorial
df_test = SNP_data[SNP_data['GENE']=='RP1-302D9.3']
print(df_test.head(2))

unique_genes = np.unique(SNP_data["GENE"])
keys = unique_genes.tolist()
dicts = {}

for k in keys:
    SNP_for_key = gene_alias_data[gene_alias_data['gene_symbol']==k]
    if (len(SNP_for_key) == 0):
        values = [k]
    else:
        values = list(SNP_for_key["alias"])
    dicts[k] = values

# print(dicts)
print({k: dicts[k] for k in list(dicts)[:5]})


