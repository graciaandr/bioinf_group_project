import allel as sc
import pandas as pd
import numpy as np
from enum import unique


SNP_data = pd.read_csv("filtered_SNPs_info.csv", sep = ",")
gene_alias_data = pd.read_csv("gene_name_alias.csv", sep = ",")
gene_data = pd.read_csv("gene_names_filtered_for_PK.csv", sep = ",")
gene_data.rename(columns={"ID": "rsID"})

mex_pop_data = pd.read_csv("freqs_mexican.csv", sep = ";")
mex_pop_data = mex_pop_data.add_prefix('MXL_')

luhya_pop_data = pd.read_csv("freqs_luhya.csv", sep = ";")
luhya_pop_data = luhya_pop_data.add_prefix('LWK_')


guj_pop_data = pd.read_csv("freqs_gujarati.csv", sep = ";")
guj_pop_data = guj_pop_data.add_prefix('GIH_')


gbr_pop_data = pd.read_csv("freqs_gbr.csv", sep = ";")
gbr_pop_data = gbr_pop_data.add_prefix('GBR_')


dai_chin_pop_data = pd.read_csv("freqs_dai_chinese.csv", sep = ";")  
dai_chin_pop_data = dai_chin_pop_data.add_prefix('CDX_')
# print(dai_chin_pop_data.head(5))

SNP_data = SNP_data.reset_index(drop=True).join(mex_pop_data)
SNP_data = SNP_data.reset_index(drop=True).join(dai_chin_pop_data)
SNP_data = SNP_data.reset_index(drop=True).join(luhya_pop_data)
SNP_data = SNP_data.reset_index(drop=True).join(gbr_pop_data)
SNP_data = SNP_data.reset_index(drop=True).join(guj_pop_data)
# print(SNP_data.head(2))

unique_genes = np.unique(gene_data["GENE"])
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
# print({k: dicts[k] for k in list(dicts)[:2]})

SNP_data['ALIAS'] = gene_data['GENE'].map(dicts)

SNP_data['ALIAS'] = pd.DataFrame([str(line).strip('[').strip(']').replace("'","") for line in SNP_data['ALIAS']])

print(SNP_data.head(5))

SNP_data.to_csv("complete_SNP_table.csv",index = False, header = True, sep = ',')



