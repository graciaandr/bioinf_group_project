# DATA MINING AND WRANGLING
# FROM VCFM,CSV, AND TXT FILES

import allel as sc
import pandas as pd
import numpy as np
import csv, sqlite3


################### functions ###################
#################################################

# filter out multi-allelic in vcf data
# vcf dataframe has 3 alt-allele columns: ALT_1, ALT_2, ALT_3
# extract rows with 1 alt allele => ALT_2 and ALT_3 is NA
def single_alt(df):
    df = df[df['ALT_2'].isna()] # if ALT_2 is NA then ALT_3 will be NA as well
    df.drop(['ALT_2', 'ALT_3'], axis=1, inplace=True) # remove ALT_2/3 columns
    df.rename(columns={"ALT_1" : "ALT"}, inplace=True) # rename column to ALT
    return df

# filter out irrelevant gene IDs
# take out those Gene IDs of the sorts "ENS\d+-ENS\d+""
def gene_id_refinement(df):
    pat = r'(\w+)\-(\w+)'
    df2 = df[df.GENEID.str.contains(pat) == False]
    return df2
    
# extract SNPs and biallelic
def snps_only(data):
    # filter single nucleotide alt and ref alleles
    data2 = data['ALT'].str.split(expand=True).apply(lambda x: x.str.len()) <= 1 # single letter column values: A/C/T/G
    data3 = data[data2.any(axis=1)]
    data4 = data3['REF'].str.split(expand=True).apply(lambda x: x.str.len()) <= 1
    data5 = data3[data4.any(axis=1)]
    return data5

# change genotype symbols to alphabet
# genotype 0 will be changed to REF allele while 1 will be changed to ALT allele
def gt_to_alp(data, start_col, end_col): # start/end col are column numbers where the genotypes are
    for pos in range(start_col, end_col):
        data.iloc[:,pos] = np.where(data.iloc[:,pos] == '0|1', data['REF'] + data['ALT'], data.iloc[:,pos])
        data.iloc[:,pos] = np.where(data.iloc[:,pos] == '1|0', data['ALT'] + data['REF'], data.iloc[:,pos])
        data.iloc[:,pos] = np.where(data.iloc[:,pos] == '1|1', data['ALT'] + data['ALT'], data.iloc[:,pos])
        data.iloc[:,pos] = np.where(data.iloc[:,pos] == '0|0', data['REF'] + data['REF'], data.iloc[:,pos])
    return data

# insert values from csv table into SQL table in database
def insert_table(csv_file, database_name, table_name):
    with sqlite3.connect(database_name) as con: # database connection using sqlite3 and cursor object
        cursor = con.cursor()
        with open (csv_file, 'r') as i: # open the csv file
            reader = csv.reader(i)
            columns = next(reader)
            # INSERT SQL statement to populate db table 
            query = "INSERT INTO " + table_name + " ({0}) VALUES ({1})" # assign query statement to a variable
            query = query.format(','.join(columns), ','.join('?' * len(columns)))
            # loop it through the whole rows in table
            for data in reader:
                cursor.execute(query, data) # call query and execute it
            con.commit()

############################################################################
####################### DATA MINING AND WRANGLING ##########################

## TABLE 1: population and samples data
pop_gt7 = pd.read_csv("raw_population_samples.tsv", sep = "\t")
# rename columns to not have space
pop_gt7.rename(columns = {"Sample name" : "Sample",
                            "Biosample ID" : "Biosample_ID",
                            "Population code" : "Pop_code",
                            "Population name" : "Pop_name",
                            "Superpopulation code" : "Superpopulation_code",
                            "Superpopulation name" : "Superpopulation_name",
                            "Population elastic ID" : "Pop_elastic_ID",
                            "Data collections" : "Data"}, inplace=True)
pop_gt7.to_csv("population_samples.csv", index = False) # save populations and samples data as csv


## TABLE 2: refine gene names data
# extract the 11 relevant columns incl. GENE ID from annotated txt file
genegt1 = pd.read_csv('raw_gene_names.txt', sep="\t", low_memory=False)
genegt1 = genegt1.iloc[:, :11]
# add foreign keys in a new column which consist of chrom, pos, ref, alt separated by colon (:)
cols = ['CHROM', 'POS', 'REF', 'ALT']
genegt1['FK_ID'] = genegt1[cols].apply(lambda row: ':'.join(row.values.astype(str)), axis=1)
# rearrange foreign key column to be the first column in table
fkey = genegt1['FK_ID']
genegt1.drop(labels=['FK_ID'], axis=1,inplace = True)
genegt1.insert(0, 'FK_ID', fkey)
# call function to filter SNPs only samples
genegt2 = snps_only(genegt1)
# remove unnecessary columns
genegt3 = genegt2.drop(['CHROM', 'FILTER', 'ALLELE', 'EFFECT'], axis=1, inplace=True)
genegt3.to_csv('22_gene_data.csv',sep = ',', index = False)
# remove irrelevant Gene ID rows
df1 = gene_id_refinement(genegt3)
df1.to_csv("filtered_22_gene_data.csv", index = False, header = True) # save as extra df
# more refinement 
df1['FK_ID']= df1['FK_ID'].str.replace(r'^chr22', '22') # replace chr22:XXX:XXX with 22:XXX:XXX
df2 = pd.read_csv('filtered_SNPs_info.csv', sep=',') # load data set that has chrom PK_ID POS REF ALT (was in vcf)
df3 = df1[df1['FK_ID'].isin(df2['PK_ID'])] # filter gene name table for IDs that are in primary key table
df3.to_csv("gene_names_filtered_for_PK.csv", index = False, header = True) # save this new gene names table


#get CHROM, ID, POS, REF, ALT from vcf file and store in csv file
sc.vcf_to_csv('chr22_pop_gt.vcf', '22_basic_info.csv', fields=['CHROM', 'POS', 'ID', 'REF', 'ALT'])


## TABLE 3: genotype data
# get genotype data of samples per population
gt1 = pd.read_csv('22_basic_info.csv')
gt2 = pd.read_csv('gt_mexican.csv', sep="\t") # gt_population.csv was generated using bcftools in linux
gt2 = gt2.iloc[: ,9:]
gt3 = pd.concat([gt1, gt2], axis=1) # combine basic info data with the genotype data
### refine genotype table
# filter out multi-allelic
gt3 = single_alt(gt3)
# SNPs only
gt4 = snps_only(gt3)
# rearrange column
id = gt4['ID']
gt4.drop(labels=['ID'], axis=1,inplace = True)
gt4.insert(0, 'PK_ID', id)
# change genotype symbols to alphabet
gt5 = gt_to_alp(gt4, 5, 102)
# remove unecessary columns because there is unique SNP ID
gt5 = gt5.drop(columns =['CHROM', 'POS', 'REF', 'ALT'])
# to csv
gt5.to_csv('converted_gt_mexican.csv', index=False)


## TABLE 4: genotype and alele frequenc
# genotype, allele freq data
freq = pd.read_csv('22_GBR_gt1.csv')
freq = freq.iloc[:,:5]
freq2 = pd.read_csv('22_GBR_freq.csv', sep = ';')
# rename columns to include which populations
freq2.rename(columns={"AF_ref":"GBR_AF_REF", "AF_alt":"GBR_AF_ALT",
                        "GT00":"GBR_GT00", "GT01":"GBR_GT01", "GT10":"GBR_GT10"},
                        inplace = True)
(pd.concat([freq, freq2], axis=1)).to_csv('22_gt_al_freq.csv', index=False)


##############################################################
################ Insert csv into SQL database ################

# connecting to sql database. if database with the name below doesn't exist, it will automatically generate a new db with that name
connection = sqlite3.connect("chromosome22_snps.db")
cursor = connection.cursor() # create cursor object
# SQL CREATE TABLE statement 
cursor.execute("""CREATE TABLE snp_table (unique_SNP_ID,CHROM,POS,rsID,REF,ALT,IMPACT,GENE,GENEID,
                MXL_AF_ref,MXL_AF_alt,MXL_AC_ref,MXL_AC_alt,MXL_GT00,MXL_GT0110,MXL_GT11,
                CDX_AF_ref,CDX_AF_alt,CDX_AC_ref,CDX_AC_alt,CDX_GT00,CDX_GT0110,CDX_GT11,
                LWK_AF_ref,LWK_AF_alt,LWK_AC_ref,LWK_AC_alt,LWK_GT00,LWK_GT0110,LWK_GT11,
                GBR_AF_ref,GBR_AF_alt,GBR_AC_ref,GBR_AC_alt,GBR_GT00,GBR_GT0110,GBR_GT11,
                GIH_AF_ref,GIH_AF_alt,GIH_AC_ref,GIH_AC_alt,GIH_GT00,GIH_GT0110,GIH_GT11,GENE_ALIAS,
                Ancestral_Allele,AA_AF,Derived_Allele,DA_AF,
                PRIMARY KEY (unique_SNP_ID));""")
# insert table values from csv file
insert_table('SNP_data_derived_AF.csv', 'chromosome22_snps.db', 'snp_table') 

connection.commit()
connection.close() #close connection