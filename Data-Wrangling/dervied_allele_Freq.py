from operator import index
import pandas as pd
import numpy as np
from math import floor, ceil

# load SNP data and derived frequency CSV files
complete_df = pd.read_csv("rounded_final_SNP_df.csv", sep = ",")
derivedAf_df= pd.read_csv("Derived_Freq.csv", sep = ";")

# remove multi-allelic SNPs
derivedAf_df = derivedAf_df[derivedAf_df['N_ALLELES']<=2]

# separate Allele Frequency columns into the respective allele and allele frequency
# A:0.999 => A | 0.999 (as two different column entries)
derivedAf_df[['Ancestral_Allele', 'AA_AF']] = derivedAf_df['Ancestral_Allele'].str.split(':', 1, expand=True)
derivedAf_df[['Derived_Allele', 'DA_AF']] = derivedAf_df['Derived_Allele1'].str.split(':', 1, expand=True)

# remove obsolete columns
derivedAf_df.drop(['CHROM', 'N_ALLELES', 'N_CHR', 'Derived_Allele1', 'Derived_Allele2'], axis=1, inplace=True)

# merge as left-join the modified derived allele frequency data frame to the SNP data 
# => left join as not all SNPs in the SNP dataframe

merged_df = pd.merge(complete_df, derivedAf_df, on='POS', how='left')


merged_df.to_csv("SNP_data_incl_derived_AF.csv",
                 index = False, 
                 header = True, 
                 sep = ',')
