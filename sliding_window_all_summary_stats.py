from operator import index
import pandas as pd
import allel as sc
import numpy as np
from numpy.lib.scimath import logn
from math import e
from statistics import mean
import matplotlib.pyplot as plt 

df = pd.read_csv("corrected_frequencies.csv", sep = ";")

# tajima d
# df = pd.read_csv('mock_freqs.csv', sep=',')
n = df.shape[0]
# print(df)
print(n)

arr = (df[['AC_ref', 'AC_alt']]).to_numpy()
taj = sc.allel.tajima_d(arr)
print("Tajima's D:",taj)

# windowed Tajima's D
def windowedTajimasD(df, n, w):
    wind_tajd = []
    for i in range(0, (n-w+1)):
        ac_window = arr[i:i+w]
        D_i = sc.allel.tajima_d(ac_window)
        wind_tajd.append(D_i)
    return(wind_tajd)

windowed_tajimasD = windowedTajimasD(df, n, w = 1000)
# print("Windowed Tajima's D:", windowed_tajimasD)

plt.plot(windowed_tajimasD,  'go', markersize = 1)
plt.title("Sliding Window - Tajima's D")
plt.ylabel("Tajima's D")
plt.xlabel('Windows')
# plt.ylim(0,0.03)
plt.show()

# df = df.iloc[17:25]
# print(df)



# Windowed Shannon Diversity
def calcWindowedShannonDiv(df, n, w = 1000):
    shannon_divs_per_w = []
    df['norm_ALT'] = df['AF_alt']/n
    df['ln_ALT'] = np.log(df['norm_ALT'])
    df['shannon'] = df['norm_ALT'] * df['ln_ALT']
    df['shannon'] = df['shannon'].fillna(0)
    for i in range(0,n-w+1): # + 1 or not?
        shannon_i = (df['shannon'][i:(i+w)]).sum() * -1
        shannon_divs_per_w.append(shannon_i)
    return(shannon_divs_per_w)

# shannondiv = calcWindowedShannonDiv(df, n, w = 1000)
# print("Shannon Diversity Vector:",shannondiv)

# plt.plot(shannondiv,  'ro', markersize = 1)
# plt.title('Sliding Window - Shannon Diversity Index')
# plt.ylabel('Shannon Diversity')
# plt.xlabel('Windows')
# plt.ylim(0,0.03)
# plt.show()


df["heterozygosity_per_snp"] = 1-((df['AF_alt'])**2 + (df['AF_ref'])**2) 

# Windowed Het. Div Index
def calcWindowedHetDiv(df, n, w = 1000):
    het_divs_per_w = []
    n = df.shape[0]
    for i in range(0,n-w+1):
        het_i = (df['heterozygosity_per_snp'][i:(i+w)])
        # print(het_i)
        het_divs_per_w.append(mean(het_i))
    return het_divs_per_w

hetDivs = calcWindowedHetDiv(df, n, w = 1000)
# print(hetDivs)

# plt.plot(hetDivs,  'bo', markersize = 1)
# plt.title('Sliding Window - Heteozygosity Diversity Index')
# plt.ylabel('Heteozygosity Diversity')
# plt.xlabel('Windows')
# plt.ylim(0,0.03)
# plt.show()
# print(df.head(18))
# print(df.shape[0])
# df.to_csv('heterozygosity_df_GBR.csv', index = False, header = True)