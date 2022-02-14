from operator import index
import pandas as pd
import allel as sc
import numpy as np
from numpy.lib.scimath import logn
from math import e
from statistics import mean
import matplotlib.pyplot as plt 

# df = pd.read_csv("freqs_gbr.csv", sep = ";")
# df = pd.read_csv("freqs_dai_chinese.csv", sep = ";")
# df = pd.read_csv("freqs_gujarati.csv", sep = ";")
# df = pd.read_csv("freqs_luhya.csv", sep = ";")
df = pd.read_csv("freqs_mexican.csv", sep = ";")

# calculate Shannon Diversity for all SNPs in a data frame
def calcShannonDiv(df):
    n = df.shape[0] # nrows of data frame
    df['norm_ALT'] = df['AF_alt']/n # normalize alt. counts
    df['ln_ALT'] = np.log(df['norm_ALT']) # take natural log (ln) of alt. counts
    df['shannon'] = df['norm_ALT'] * df['ln_ALT'] # multiply
    df['shannon'] = df['shannon'].fillna(0) # for cases ln(0) = -inf, set the multiplication as 0
    shannondiv = df['shannon'].sum() * -1 # calculate sum and multiply with -1 to get shannon diversity
    df.drop(['norm_ALT', 'ln_ALT'], axis=1, inplace=True)
    return(shannondiv)

shannondiv = calcShannonDiv(df)
print("Shannon Diversity Index:", shannondiv)

# Expected Heterozygosity per locus / SNP
### Formula taken from:
### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5295611/pdf/671.pdf

def calcHeterozygosity(df):
    df["heterozygosity_per_snp"] = 1-((df['AF_alt'])**2 + (df['AF_ref'])**2) # calculate exp. het. according to formula (always only 2 alleles)
    df.to_csv('heterozygosity_df_GBR.csv', index = False, header = True)
    return(df)

# Tajima's D
arr = (df[['AC_ref', 'AC_alt']]).to_numpy() # take allele counts and save as numpy array
taj = sc.allel.tajima_d(arr) # calculate tajima's D using scikit-allel function
print("Tajima's D:", taj)


### Sliding Window Summary Stats ###

# Windowed Shannon Diversity
def calcWindowedShannonDiv(df, w = 10000):
    shannon_divs_per_w = []
    n = df.shape[0] # nrows of data frame
    df['norm_ALT'] = df['AF_alt']/n
    df['ln_ALT'] = np.log(df['norm_ALT'])
    df['shannon'] = df['norm_ALT'] * df['ln_ALT']
    df['shannon'] = df['shannon'].fillna(0)
    for i in range(0,n-w+1): 
        shannon_i = (df['shannon'][i:(i+w)]).sum() * -1 # calculate shannon diverseity per window
        shannon_divs_per_w.append(shannon_i) # store those values per window in a list/vector
    return(shannon_divs_per_w)

shannondiv = calcWindowedShannonDiv(df, w = 3)

# Plot windowed shannon diversity
# plt.plot(shannondiv,  'ro', markersize = 1)
# plt.title('Sliding Window - Shannon Diversity Index')
# plt.ylabel('Shannon Diversity')
# plt.xlabel('Windows')
# plt.show()


# Windowed Het. Div Index
def calcWindowedHetDiv(df, w = 1000):
    het_divs_per_w = []
    n = df.shape[0]
    df = calcHeterozygosity(df)
    for i in range(0,n-w+1):
        ### mean or median ? 
        het_i = mean(df['heterozygosity_per_snp'][i:(i+w)]) # take mean of the exp. het. for chosen window
        het_divs_per_w.append(het_i)
    return het_divs_per_w

hetDivs = calcWindowedHetDiv(df, w = 3)

# plt.plot(hetDivs,  'bo', markersize = 1)
# plt.title('Sliding Window - Heteozygosity Diversity Index')
# plt.ylabel('Heteozygosity Diversity')
# plt.xlabel('Windows')
# plt.ylim(0,0.03)
# plt.show()



## windowed Tajima's D
def windowedTajimasD(df, w):
    n = df.shape[0]
    wind_tajd = []
    arr = (df[['AC_ref', 'AC_alt']]).to_numpy() 
    for i in range(0, (n-w+1)):
        ac_window = arr[i:i+w] # set window of allele counts array
        D_i = sc.allel.tajima_d(ac_window)
        wind_tajd.append(D_i)
    return(wind_tajd)

windowed_tajimasD = windowedTajimasD(df, w = 3)

# plt.plot(windowed_tajimasD,  'go', markersize = 1)
# plt.title("Sliding Window - Tajima's D")
# plt.ylabel("Tajima's D")
# plt.xlabel('Windows')
# plt.show() 