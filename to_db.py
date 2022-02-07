import allel as sc
import pandas as pd
import csv, sqlite3

#get CHROM, ID, POS, REF, ALT from vcf file and store in csv file
sc.vcf_to_csv('chr22_pop_gt.vcf', '22_basic_info.csv', fields=['CHROM', 'POS', 'ID', 'REF', 'ALT'])

# extract the 11 relevant columns incl. gene_names from annotated txt file
gene_names = pd.read_csv('gene_names.txt', sep="\t", header=None)
first_n_columns = gene_names.iloc[:, :11]
first_n_columns.to_csv('22_gene_names.csv',sep = ',', index = None, header = False)

# store population info in csv file 
pop_df = pd.read_csv("population_samples.tsv", sep = "\t")
pop_df.to_csv("population_samples.csv", index = None)

# get genotype data of GBR samples
gt_df = pd.read_csv('GBR_gt.csv', sep="\t")
gt_df2 = gt_df.iloc[: ,9:]
# combine csv files
df1 = pd.read_csv('22_basic_info.csv')
(pd.concat([df1, gt_df2], axis=1).to_csv('22_GBR_genotype.csv', index=False, na_rep='N/A'))

############

con = sqlite3.connect("chr22.db")
cur = con.cursor()

# population_sample table
all_pop = pd.read_csv("population_samples.csv")
all_pop.to_sql("pop_samples", con, index = False)

# gene_names table
gene_names = pd.read_csv("22_gene_names.csv")
gene_names.to_sql("gene_names", con, index=False)

# genotype table
GBR_gt = pd.read_csv("22_GBR_genotype.csv")
GBR_gt.to_sql('GBR_genotype',con, index=False)

con.commit()
con.close()