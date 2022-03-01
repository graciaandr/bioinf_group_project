import pandas as pd
import numpy as np
import allel
from matplotlib import pyplot as plt
from math import ceil

df_gbr= pd.read_csv("freqs_gbr.csv", sep = ";")
df_mexican= pd.read_csv("freqs_mexican.csv", sep = ";")
df_chinese= pd.read_csv("freqs_dai_chinese.csv", sep = ";")
df_luhya= pd.read_csv("freqs_luhya.csv", sep = ";")
df_gujarati= pd.read_csv("freqs_gujarati.csv", sep = ";")


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

# ######### fst calculation function ################
def calculate_fst(population1, population2):
  # pass array for 2 pops
  pop1,pop2= allel.hudson_fst(population1, population2)
  # fst = pop1/pop2 diviide by sum of  allele coumts for each population
  fst = np.sum(pop1) / np.sum(pop2)
  return (fst)
  
# a = calculate_fst(array_chinese, array_luhya)    
# print(a)  
  
########## compute sliding window ##########
def slide_for_fst(pop1,pop2, window_size = None):
    p1 = extract_and_makearray(pop1)
    p2 = extract_and_makearray(pop2)
    if (window_size == None):
      n = (p1.shape[0] + p2.shape[0])/2
      window_size = ceil(n/100)   
    print(window_size)
    # call function in moving hudson allel function 
    slide = allel.moving_hudson_fst(p1, p2,
                                    size= window_size)
    return slide

windows_average = slide_for_fst(df_gbr, df_mexican)
plt.plot(windows_average, 'purple')
plt.title("Sliding Window - FST stats")
plt.xlabel("Windows")
plt.ylabel("Average FST index")
plt.savefig('FST.png')
plt.show()  