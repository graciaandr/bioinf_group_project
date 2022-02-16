from cmath import nan
import allel as sc
import pandas as pd
import numpy as np
from enum import unique

# load all relevant data sets
SNP_data = pd.read_csv("filtered_SNPs_info.csv", sep = ",")
gene_alias_data = pd.read_csv("gene_name_alias.csv", sep = ",")
gene_data = pd.read_csv("gene_names_filtered_for_PK.csv", sep = ",")
gene_data = gene_data.rename(columns={"ID": "rsID"})
SNP_data = SNP_data.rename(columns={"PK_ID": "unique_SNP_ID"})

# drop duplicate columns (POS, REF, ALT are already in gene data, FK_ID is identical to PK_ID)
SNP_data.drop(['POS','REF', 'ALT'], inplace = True, axis = 1)
gene_data.drop(['FK_ID'], inplace = True, axis = 1)


# load all population data sets with allele frequencies
# and rename the columns by adding population code prefix
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

# merge / join all data sets
SNP_data = SNP_data.reset_index(drop=True).join(gene_data)
SNP_data = SNP_data.reset_index(drop=True).join(mex_pop_data)
SNP_data = SNP_data.reset_index(drop=True).join(dai_chin_pop_data)
SNP_data = SNP_data.reset_index(drop=True).join(luhya_pop_data)
SNP_data = SNP_data.reset_index(drop=True).join(gbr_pop_data)
SNP_data = SNP_data.reset_index(drop=True).join(guj_pop_data)

# take unique gene names and add gene aliases as extra column based on dictionary
# that reveals the gene aliases
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

# map gene aliases to genes in gene data frame and add to SNP data frame
SNP_data['GENE_ALIAS'] = gene_data['GENE'].map(dicts)
SNP_data['GENE_ALIAS'] = pd.DataFrame(
    [str(line).strip('[').strip(']').replace("'","") for line in SNP_data['GENE_ALIAS']])

# set all missing rsIDs to NAN instead of '.'
mask = SNP_data.rsID.eq('.')                                                                                     
SNP_data.loc[mask, 'rsID'] = np.nan   

# save this merged and finally edited data frame as CSV table
SNP_data.to_csv("final_complete_SNP_table.csv",index = False, header = True, sep = ',')



