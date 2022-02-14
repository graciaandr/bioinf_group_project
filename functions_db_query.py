import sqlite3, csv
import pandas as pd
import numpy as np
import re

# change genotype symbols to alphabet
def gt_to_alp(data, start_col, end_col):
    for pos in range(start_col, end_col):
        data.iloc[:,pos] = np.where(data.iloc[:,pos] == '0|1', data['REF'] + data['ALT'], data.iloc[:,pos])
        data.iloc[:,pos] = np.where(data.iloc[:,pos] == '1|0', data['ALT'] + data['REF'], data.iloc[:,pos])
        data.iloc[:,pos] = np.where(data.iloc[:,pos] == '1|1', data['ALT'] + data['ALT'], data.iloc[:,pos])
        data.iloc[:,pos] = np.where(data.iloc[:,pos] == '0|0', data['REF'] + data['REF'], data.iloc[:,pos])
    return data

# filter out multi-allelic
def single_alt(df):
    df = df[df['ALT_2'].isna()] 
    df.drop(['ALT_2', 'ALT_3'], axis=1, inplace=True)
    df.rename(columns={"ALT_1" : "ALT"}, inplace=True)
    return df

# filter out irrelevant gene IDs
def gene_id_refinement(df):
    # take out those Gene IDs of the sorts "ENS\d+-ENS\d+""
    pat = r'(\w+)\-(\w+)'
    df2 = df[df.GENEID.str.contains(pat) == False]
    return df2

# extract SNPs and biallelic
def snps_only(data):
    data2 = data['ALT'].str.split(expand=True).apply(lambda x: x.str.len()) <= 1
    data3 = data[data2.any(axis=1)]
    data4 = data3['REF'].str.split(expand=True).apply(lambda x: x.str.len()) <= 1
    data5 = data3[data4.any(axis=1)]
    return data5

# insert values into table in db
def insert_table(csv_file, table_name):
    with sqlite3.connect("chr22.db") as con:
        cursor = con.cursor()
        with open (csv_file, 'r') as i:
            reader = csv.reader(i)
            columns = next(reader) 
            query = "INSERT INTO " + table_name + " ({0}) VALUES ({1})"
            query = query.format(','.join(columns), ','.join('?' * len(columns)))
            for data in reader:
                cursor.execute(query, data)
            con.commit()

# load data with gene names for chromosome 22
df = pd.read_csv('22_gene_data.csv', sep=',')
df1 = gene_id_refinement(df) # remove irrelevant vant Gene ID rows
df1.to_csv("filtered_22_gene_data.csv", index = False, header = True) # save as extra df

# df1 = pd.read_csv('filtered_22_gene_data.csv', sep=',')
df1['FK_ID']= df1['FK_ID'].str.replace(r'^chr22', '22') # replace chr22:XXX:XXX with 22:XXX:XXX

df2 = pd.read_csv('filtered_SNPs_info.csv', sep=',') # load data set that has chrom PK_ID POS REF ALT (was in vcf)
df3 = df1[df1['FK_ID'].isin(df2['PK_ID'])] # filter gene name table for IDs that are in primary key table
# print(df3.shape)
df3.to_csv("gene_names_filtered_for_PK.csv", index = False, header = True) # save this new gene names table


############### query functions ################

def use_id(search_value):
    with sqlite3.connect("chr22.db") as connection:
        search_value = search_value.lower()
        cursor = connection.cursor()
        result = cursor.execute(f"SELECT * FROM gene_data WHERE ID = '{search_value}'").fetchall()
    return result

def use_gene(search_value):
    with sqlite3.connect("chr22.db") as connection:
        search_value = search_value.upper()
        cursor = connection.cursor()
        result = cursor.execute(f"SELECT * FROM gene_data WHERE GENE = '{search_value}'").fetchall()
    return result

def use_pos(search_value):
    with sqlite3.connect("chr22.db") as connection:
        cursor = connection.cursor()
        result = cursor.execute(f"SELECT * FROM gt_al_freq WHERE POS = '{search_value}'").fetchall()
    return result

def search_db(search_type, search_value):
    if search_type == "snp_id":
        return use_id(search_value)
    elif search_type == "gene_name":
        return use_gene(search_value)
    elif search_type == "genomic_coordinate":
        return use_pos(search_value)
    else:
        return "404: Not Found"