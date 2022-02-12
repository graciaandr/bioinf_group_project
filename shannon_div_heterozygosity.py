from operator import index
import pandas as pd
import csv
import numpy as np
from numpy.lib.scimath import logn
from math import e
from statistics import mean

df = pd.read_csv("corrected_frequencies.csv", sep = ";")
n = df.shape[0]
print(n)

# Shannon Diversity
def calcShannonDiv(df, n):
    df['norm_ALT'] = df['AF_alt']/n
    df['ln_ALT'] = np.log(df['norm_ALT'])
    df['shannon'] = df['norm_ALT'] * df['ln_ALT']
    df['shannon'] = df['shannon'].fillna(0)
    shannondiv = df['shannon'].sum() * -1
    return(shannondiv)

shannondiv = calcShannonDiv(df, n)
print("Shannon Diversity Index:",shannondiv)

# Expected heterozygosity index
### Sources:
### https://popgen.nescent.org/StartSNP.html#genetic-diversity-observed-and-expected-heterozygosity
### according to formula in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5295611/pdf/671.pdf

df["heterozygosity_per_snp"] = 1-((df['AF_alt'])**2 + (df['AF_ref'])**2) 

print('Heterozygosity Index:' ,mean(df["heterozygosity_per_snp"]))

print(df.head(18))
print(df.shape[0])
df.to_csv('heterozygosity_df_GBR.csv', index = False, header = True)