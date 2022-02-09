import allel as sc
import pandas as pd
import database_query as dbq
import csv, sqlite3

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

# extract the 11 relevant columns incl. GENE ID from annotated txt file
genedf1 = pd.read_csv('raw_gene_names.txt', sep="\t", low_memory=False)
genedf1 = genedf1.iloc[:, :11]
# insert foreign key column
cols = ['CHROM', 'POS', 'REF', 'ALT']
genedf1['foreign_key'] = genedf1[cols].apply(lambda row: ':'.join(row.values.astype(str)), axis=1)
# rearrange foreign key column
fkey = genedf1['foreign_key']
genedf1.drop(labels=['foreign_key'], axis=1,inplace = True)
genedf1.insert(0, 'foreign_key', fkey)
print(genedf1.shape)
# ##### refine data
# keep single ALT
genedf2 = genedf1['ALT'].str.split(expand=True).apply(lambda x: x.str.len()) <= 1
genedf3 = genedf1[genedf2.any(axis=1)]
print(genedf3.shape)
# # keep SNPs only (single REF)
genedf4 = genedf3['REF'].str.split(expand=True).apply(lambda x: x.str.len()) <= 1
genedf5 = genedf3[genedf4.any(axis=1)]
print(genedf5.shape)
# # remove unnecessary columns
genedf5.drop(['CHROM', 'FILTER', 'ALLELE', 'EFFECT'], axis=1, inplace=True)
genedf5.to_csv('22_gene_names.csv',sep = ',', index = None)

#get CHROM, ID, POS, REF, ALT from vcf file and store in csv file
sc.vcf_to_csv('chr22_pop_gt.vcf', '22_basic_info.csv', fields=['CHROM', 'POS', 'ID', 'REF', 'ALT'])

# get genotype data of GBR samples
gt1 = pd.read_csv('22_basic_info.csv')
gt2 = pd.read_csv('GBR_gt.csv', sep="\t")
gt2 = gt2.iloc[: ,9:]
gt3 = pd.concat([gt1, gt2], axis=1)
### refine GT table
# filter out multi-allelic
gt3 = gt3[gt3['ALT_2'].isna()] 
gt3.drop(['ALT_2', 'ALT_3'], axis=1, inplace=True)
gt3.rename(columns={"ALT_1" : "ALT"}, inplace=True)
# keep single ALT
gt4 = gt3['ALT'].str.split(expand=True).apply(lambda x: x.str.len()) <= 1
gt5 = gt3[gt4.any(axis=1)]
# keep SNPs only (single REF)
gt6 = gt5['REF'].str.split(expand=True).apply(lambda x: x.str.len()) <= 1
gt7 = gt5[gt6.any(axis=1)]
# rearrange column
id = gt7['ID']
gt7.drop(labels=['ID'], axis=1,inplace = True)
gt7.insert(0, 'PK_ID', id)
# to csv
gt7.to_csv('22_GBR_genotype.csv', index=False)


################

con = sqlite3.connect("chr22.db")
cursor = con.cursor()

# population_sample table
cursor.execute("""CREATE TABLE pop_samples (Sample,Sex, Biosample_ID, Pop_code, Pop_name, 
                Superpopulation_code, Superpopulation_name, Pop_elastic_ID, Data, 
                PRIMARY KEY (Sample));""")
dbq.insert_table('population_samples.csv', 'pop')

# gene_names table
cursor.execute("""CREATE TABLE gene_names (foreign_key, POS, ID, REF, ALT, IMPACT, GENE, GENEID);""")
dbq.insert_table('22_gene_names.csv', 'gene_names')

# genotype table
GBR_gt = pd.read_csv("22_GBR_genotype.csv")
GBR_gt.to_sql('GBR_genotype',con, index=False)

con.commit()
con.close()