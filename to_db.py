import allel as sc
import pandas as pd
import numpy as np
from pyparsing import dbl_quoted_string
import functions_db_query as dbq
import csv, sqlite3

## population sample data
# store population info in csv file 
pop_gt7 = pd.read_csv("raw_population_samples.tsv", sep = "\t")
pop_gt7.rename(columns = {"Sample name" : "Sample",
                            "Biosample ID" : "Biosample_ID",
                            "Population code" : "Pop_code",
                            "Population name" : "Pop_name",
                            "Superpopulation code" : "Superpopulation_code",
                            "Superpopulation name" : "Superpopulation_name",
                            "Population elastic ID" : "Pop_elastic_ID",
                            "Data collections" : "Data"}, inplace=True)
pop_gt7.to_csv("population_samples.csv", index = False)


## refine gene names data
# extract the 11 relevant columns incl. GENE ID from annotated txt file
genegt1 = pd.read_csv('raw_gene_names.txt', sep="\t", low_memory=False)
genegt1 = genegt1.iloc[:, :11]
# insert foreign key column
cols = ['CHROM', 'POS', 'REF', 'ALT']
genegt1['FK_ID'] = genegt1[cols].apply(lambda row: ':'.join(row.values.astype(str)), axis=1)
# rearrange foreign key column
fkey = genegt1['FK_ID']
genegt1.drop(labels=['FK_ID'], axis=1,inplace = True)
genegt1.insert(0, 'FK_ID', fkey)
# SNPs only
genegt2 = dbq.snps_only(genegt1)
# # remove unnecessary columns
genegt3 = genegt2.drop(['CHROM', 'FILTER', 'ALLELE', 'EFFECT'], axis=1, inplace=True)
genegt3.to_csv('22_gene_data.csv',sep = ',', index = False)
# remove irrelevant vant Gene ID rows
df1 = dbq.gene_id_refinement(genegt3)
df1.to_csv("filtered_22_gene_data.csv", index = False, header = True) # save as extra df
# df1 = pd.read_csv('filtered_22_gene_data.csv', sep=',')
df1['FK_ID']= df1['FK_ID'].str.replace(r'^chr22', '22') # replace chr22:XXX:XXX with 22:XXX:XXX
df2 = pd.read_csv('filtered_SNPs_info.csv', sep=',') # load data set that has chrom PK_ID POS REF ALT (was in vcf)
df3 = df1[df1['FK_ID'].isin(df2['PK_ID'])] # filter gene name table for IDs that are in primary key table
df3.to_csv("gene_names_filtered_for_PK.csv", index = False, header = True) # save this new gene names table


#get CHROM, ID, POS, REF, ALT from vcf file and store in csv file
sc.vcf_to_csv('chr22_pop_gt.vcf', '22_basic_info.csv', fields=['CHROM', 'POS', 'ID', 'REF', 'ALT'])


# get genotype data of samples per population
gt1 = pd.read_csv('22_basic_info.csv')
gt2 = pd.read_csv('gt_mexican.csv', sep="\t")
gt2 = gt2.iloc[: ,9:]
gt3 = pd.concat([gt1, gt2], axis=1)
### refine GT table
# filter out multi-allelic
gt3 = dbq.single_alt(gt3)
# SNPs only
gt4 = dbq.snps_only(gt3)
# rearrange column
id = gt4['ID']
gt4.drop(labels=['ID'], axis=1,inplace = True)
gt4.insert(0, 'PK_ID', id)
# change genotype symbols to alphabet
gt5 = dbq.gt_to_alp(gt4, 5, 102)
gt5 = gt5.drop(columns =['CHROM', 'POS', 'REF', 'ALT'])
# to csv
gt5.to_csv('converted_gt_mexican.csv', index=False)


# gt al freq data
freq = pd.read_csv('22_GBR_gt1.csv')
freq = freq.iloc[:,:5]
freq2 = pd.read_csv('22_GBR_freq.csv', sep = ';')
freq2.rename(columns={"AF_ref":"GBR_AF_REF", "AF_alt":"GBR_AF_ALT",
                        "GT00":"GBR_GT00", "GT01":"GBR_GT01", "GT10":"GBR_GT10"},
                        inplace = True)
(pd.concat([freq, freq2], axis=1)).to_csv('22_gt_al_freq.csv', index=False)


################

con = sqlite3.connect("chr22.db")
cursor = con.cursor()

# population_sample table
cursor.execute("""CREATE TABLE pop_samples (Sample,Sex, Biosample_ID, Pop_code, Pop_name, 
                Superpopulation_code, Superpopulation_name, Pop_elastic_ID, Data, 
                PRIMARY KEY (Sample));""")
dbq.insert_table('population_samples.csv', 'chr22.db', 'pop_samples')

# genotype table
cursor.execute("""CREATE TABLE gt_al_freq (PK_ID,CHROM,POS,REF,ALT,GBR_AF_REF,GBR_AF_ALT,GBR_GT00,GBR_GT01,GBR_GT10,
                PRIMARY KEY (PK_ID));""")
dbq.insert_table('22_gt_al_freq.csv', 'chr22.db', 'gt_al_freq')

# gene_names table
cursor.execute("""CREATE TABLE gene_data (FK_ID, POS, ID, REF, ALT, IMPACT, GENE, GENEID,
                FOREIGN KEY (FK_ID) REFERENCES gt_al_freq(PK_ID));""")
dbq.insert_table('22_gene_data.csv', 'chr22.db', 'gene_data')




con.commit()
con.close()