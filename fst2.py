import pandas as pd
import numpy as np
import allel

df_gbr= pd.read_csv("/data/scratch/bt211065/2022_02_01_Groupproject/test/gbr_freq.csv")

df_mexican= pd.read_csv("/data/scratch/bt211065/2022_02_01_Groupproject/test/mexican_freq.csv")

# extract ref and alt 
#allele_count_gbr = df_gbr.loc[:,['AC_ref', 'AC_alt']]
# no nan values

#total variants gbr 
#total_variantsgbr = allele_count_gbr.shape[0]

#exctract ref and alt
#allele_count_mexican =df_mexican.loc[:,['AC_ref', 'AC_alt']] 
#no nan values 

#total varaints mxl
#total_variantsmxl = allele_count_mexican.shape[0]

# gbr to array
#gbr_array = allele_count_gbr.to_numpy()
# reshape to n_varaints, n_alleles
#reshape_gbr = gbr_array.reshape(total_variantsgbr,2)

# mxl to array
#mxl_array= allele_count_gbr.to_numpy()
#reshape mxl 
#reshape_mxl = mxl_array.reshape(total_variantsmxl,2)


# pass allele count to hudson 
#pop1,pop2= allel.hudson_fst(reshape_gbr,reshape_mxl)
# fst = pop1/pop2 diviude allele coumts
#fst = np.sum(pop1) / np.sum(pop2)
#print("FST VALUE",fst)
#print(pop1, pop2)


df_chinese= pd.read_csv("/data/scratch/bt211065/2022_02_01_Groupproject/test/chinese_freq.csv")

df_luhya= pd.read_csv("/data/scratch/bt211065/2022_02_01_Groupproject/test/luhya_freq.csv")

df_gujarati= pd.read_csv("/data/scratch/bt211065/2022_02_01_Groupproject/test/gujarati_freq.csv")


#### make function to extract ref and alt and convert to numpy####

def extract_and_makearray(accept_dataframe):
  # get ref and alt
  allele_count = accept_dataframe.loc[:,['AC_ref', 'AC_alt']]
  # get number of variants 
  total_variants = allele_count.shape[0]
  allelecount_array = allele_count.to_numpy()
  # reshape to n_varaints, n_alleles
  reshape_alelle = allelecount_array.reshape(total_variants,2)
  
  return reshape_alelle
  
array_chinese = extract_and_makearray(df_chinese)
array_luhya = extract_and_makearray(df_luhya)
array_gujarati = extract_and_makearray(df_gujarati)
array_gbr=extract_and_makearray(df_gbr)
array_mexican=extract_and_makearray(df_mexican)

######### fst calculation function ################



def calculate_fst(population1, population2):
  # pass array for 2 pops
  pop1,pop2= allel.hudson_fst(population1, population2)
  # fst = pop1/pop2 diviide by sum of  allele coumts for each population
  fst = np.sum(pop1) / np.sum(pop2)
  #print(pop1,pop1)
  print("The FST value is",fst)
  
  
#fst_chinese_luhya=calculate_fst(array_chinese,array_luhya)
#fst_chinese_gujarati=calculate_fst(array_chinese,array_gujarati)
#fst_chinese_gbr=calculate_fst(array_chinese,array_gbr)
#fst_chinese_mexican=calculate_fst(array_chinese,array_gujarati)


for one in [array_chinese, array_luhya, array_gujarati,array_gbr,array_mexican]:
  for two in [array_chinese, array_luhya, array_gujarati,array_gbr,array_mexican]:
      print(calculate_fst(one,two))
    
  
  
  



  