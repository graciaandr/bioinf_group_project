import allel as sc
import pandas as pd
import database_query as dbq
import csv, sqlite3

#get CHROM, ID, POS, REF, ALT from vcf file and store in csv file
sc.vcf_to_csv('chr22_pop_gt.vcf', '22_basic_info.csv', fields=['CHROM', 'POS', 'ID', 'REF', 'ALT'])

# extract the 11 relevant columns incl. gene_names from annotated txt file
gene_names = pd.read_csv('raw_gene_names.txt', sep="\t", low_memory=False)
gene_df = gene_names.iloc[:, :11]
# insert foreign key column
cols = ['CHROM', 'POS', 'REF', 'ALT']
gene_df['foreign_key'] = gene_df[cols].apply(lambda row: ':'.join(row.values.astype(str)), axis=1)
# rearrange foreign key column
fkey = gene_df['foreign_key']
gene_df.drop(labels=['foreign_key'], axis=1,inplace = True)
gene_df.insert(0, 'foreign_key', fkey)
# remove unnecessary columns
gene_df.drop(['CHROM', 'POS', 'FILTER', 'ALLELE', 'EFFECT'], axis=1, inplace=True)
gene_df.to_csv('22_gene_names.csv',sep = ',', index = None)

# store population info in csv file 
pop_df = pd.read_csv("raw_population_samples.tsv", sep = "\t")
pop_df.rename(columns = {"Sample name" : "Sample",
                            "Biosample ID" : "Biosample_ID",
                            "Population code" : "Pop_code",
                            "Population name" : "Pop_name",
                            "Superpopulation code" : "Superpopulation_code",
                            "Superpopulation name" : "Superpopulation_name",
                            "Population elastic ID" : "Pop_elastic_ID",
                            "Data collections" : "Data"}, inplace=True)
pop_df.to_csv("population_samples.csv", index = None)

# get genotype data of GBR samples
gt_df1 = pd.read_csv('22_basic_info.csv')
gt_df1.drop(['ALT_2', 'ALT_3'], axis=1, inplace=True)
gt_df1.rename(columns={"ALT_1" : "ALT"}, inplace=True)
gt_df2 = pd.read_csv('GBR_gt.csv', sep="\t")
gt_df2 = gt_df2.iloc[: ,9:]
# rearrange column
id = gt_df1['ID']
gt_df1.drop(labels=['ID'], axis=1,inplace = True)
gt_df1.insert(0, 'PK_ID', id)
# combine csv
(pd.concat([gt_df1, gt_df2], axis=1).to_csv('22_GBR_genotype.csv', index=False, na_rep='N/A'))

################

con = sqlite3.connect("chr22.db")
cursor = con.cursor()

# population_sample table
cursor.execute("""CREATE TABLE pop_samples (Sample,Sex, Biosample_ID, Pop_code, Pop_name, 
                Superpopulation_code, Superpopulation_name, Pop_elastic_ID, Data, 
                PRIMARY KEY (Sample));""")
dbq.insert_table('population_samples.csv', 'pop')

# gene_names table
cursor.execute("""CREATE TABLE gene_names (foreign_key, ID, REF, ALT, IMPACT, GENE, GENEID,
                FOREIGN KEY (foreign_key));""")
dbq.insert_table('22_gene_names.csv', 'gene_names')

# genotype table
GBR_gt = pd.read_csv("22_GBR_genotype.csv")
GBR_gt.to_sql('GBR_genotype',con, index=False)

con.commit()
con.close()