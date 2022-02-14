import pandas as pd
import numpy as np
import allel

df_gbr= pd.read_csv("/data/scratch/bt211065/2022_02_01_Groupproject/test/gbr_freq.csv")

df_mexican= pd.read_csv("/data/scratch/bt211065/2022_02_01_Groupproject/test/mexican_freq.csv")

# extract ref and alt 
allele_count_gbr = df_gbr.loc[1:20,['AC_ref', 'AC_alt']]
# no nan values

#total variants gbr 
total_variantsgbr = allele_count_gbr.shape[0]

#exctract ref and alt
allele_count_mexican =df_mexican.loc[1:20,['AC_ref', 'AC_alt']] 
#no nan values 

#total varaints mxl
total_variantsmxl = allele_count_mexican.shape[0]

# gbr to array
gbr_array = allele_count_gbr.to_numpy()
# reshape to n_varaints, n_alleles
reshape_gbr = gbr_array.reshape(total_variantsgbr,2)

# mxl to array
mxl_array= allele_count_gbr.to_numpy()
#reshape mxl 
reshape_mxl = mxl_array.reshape(total_variantsmxl,2)


# pass allele count to hudson 
pop1,pop2= allel.hudson_fst(reshape_gbr,reshape_mxl)
# fst = pop1/pop2 diviude allele coumts
fst = np.sum(pop1) / np.sum(pop2)
print("FST VALUE",fst)
#print(pop1, pop2)

