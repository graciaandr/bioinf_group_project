from operator import index
import pandas as pd
import allel as sc
import numpy as np
from numpy.lib.scimath import logn
from math import e
from statistics import mean
import matplotlib.pyplot as plt 
from math import floor, ceil

# turning off warnings for now:
import sys
import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")
    
# df = pd.read_csv("freqs_gbr.csv", sep = ";")
# df = pd.read_csv("freqs_dai_chinese.csv", sep = ";")
# df = pd.read_csv("freqs_gujarati.csv", sep = ";")
# df = pd.read_csv("freqs_luhya.csv", sep = ";")
df = pd.read_csv("freqs_mexican.csv", sep = ";")
df = df.head(500)

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
    # df.to_csv('heterozygosity_df_GBR.csv', index = False, header = True)
    return(df)

# Tajima's D
arr = (df[['AC_ref', 'AC_alt']]).to_numpy() # take allele counts and save as numpy array
taj = sc.allel.tajima_d(arr) # calculate tajima's D using scikit-allel function
print("Tajima's D:", taj)


### Sliding Window Summary Stats ###

# Windowed Shannon Diversity
def calcWindowedShannonDiv(df, w = None):
    n = df.shape[0] # nrows of data frame
    if (w == None):
        w = ceil(n/10)
    shannon_divs_per_w = []
    df['norm_ALT'] = df['AF_alt']/n
    df['ln_ALT'] = np.log(df['norm_ALT'])
    df['shannon'] = df['norm_ALT'] * df['ln_ALT']
    df['shannon'] = df['shannon'].fillna(0)
    for i in range(0,n-w+1): 
        shannon_i = (df['shannon'][i:(i+w)]).sum() * -1 # calculate shannon diverseity per window
        shannon_divs_per_w.append(shannon_i) # store those values per window in a list/vector
    return(shannon_divs_per_w)

shannondiv = calcWindowedShannonDiv(df)

## Histogram
# hist, bin_edges = np.histogram(shannondiv, 
#                     range=(np.nanmin(shannondiv),
#                            np.nanmax(shannondiv)),
#                     bins='fd')                                        
# plt.hist(hist, bins='auto')  
# plt.xlabel('Shannon Div')
# plt.ylabel('Windows')
# plt.title('Histogram of Shannon Div')
# plt.grid(True)
# plt.show()

## Boxplot
# fig1, ax1 = plt.subplots()
# ax1.set_title('Basic Plot')
# ax1.boxplot(shannondiv)
# plt.show()

## Linegraph plot of windowed shannon diversity
plt.plot(shannondiv,  'r', markersize = 1)
plt.title('Sliding Window - Shannon Diversity Index')
plt.ylabel('Shannon Diversity')
plt.xlabel('Windows')
plt.show()


# Windowed Het. Div Index
def calcWindowedHetDiv(df, w = None):
    het_divs_per_w = []
    n = df.shape[0]
    if (w == None):
        w = ceil(n/10)
    df = calcHeterozygosity(df)
    for i in range(0,n-w+1):
        ### mean or median ? 
        het_i = mean(df['heterozygosity_per_snp'][i:(i+w)]) # take mean of the exp. het. for chosen window
        het_divs_per_w.append(het_i)
    return het_divs_per_w

hetDivs = calcWindowedHetDiv(df)

plt.plot(hetDivs,  'b', markersize = 1)
plt.title('Sliding Window - Heteozygosity Diversity Index')
plt.ylabel('Heteozygosity Diversity')
plt.xlabel('Windows')
plt.show()


## windowed Tajima's D
def windowedTajimasD(df, w = None):
    n = df.shape[0]
    if (w == None):
        w = ceil(n/10)
    wind_tajd = []
    arr = (df[['AC_ref', 'AC_alt']]).to_numpy() 
    for i in range(0, (n-w+1)):
        ac_window = arr[i:i+w] # set window of allele counts array
        D_i = sc.allel.tajima_d(ac_window)
        wind_tajd.append(D_i)
    return(wind_tajd)

windowed_tajimasD = windowedTajimasD(df)
plt.plot(windowed_tajimasD,  'g', markersize = 1)
plt.title("Sliding Window - Tajima's D")
plt.ylabel("Tajima's D")
plt.xlabel('Windows')
plt.show() 