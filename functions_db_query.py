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
    df = df[df['ALT_2'].isna()] #extract rows with 1 ALT
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
def insert_table(csv_file, database_name, table_name):
    with sqlite3.connect(database_name) as con:
        cursor = con.cursor()
        with open (csv_file, 'r') as i:
            reader = csv.reader(i)
            columns = next(reader) 
            query = "INSERT INTO " + table_name + " ({0}) VALUES ({1})"
            query = query.format(','.join(columns), ','.join('?' * len(columns)))
            for data in reader:
                cursor.execute(query, data)
            con.commit()



############### query functions ################

def use_id(search_value):
    with sqlite3.connect("test.db") as connection:
        search_value = search_value.lower()
        cursor = connection.cursor()
        result = cursor.execute(f"SELECT * FROM gene_data WHERE ID = '{search_value}'").fetchall()
    return result

def use_gene(search_value):
    with sqlite3.connect("test.db") as connection:
        search_value = search_value.upper()
        cursor = connection.cursor()
        result = cursor.execute(f"SELECT * FROM gene_data WHERE GENE = '{search_value}'").fetchall()
    return result

def use_pos(search_value):
    # 2 types of input: ranged number separated by hyphen or a position number
    range = r"(^\d+)\-(\d+$)"
    position = r"(^\d+$)"
    with sqlite3.connect("test.db") as connection:
        cursor = connection.cursor()
        if re.search(range, f"{search_value}") is not None:
            input_nums = (str(search_value)).split('-')
            result = cursor.execute(f"SELECT * FROM gene_data WHERE POS BETWEEN '{input_nums[0]}' AND '{input_nums[1]}'").fetchall()
            return result
        elif re.search(position, f"{search_value}") is not None:
            result = cursor.execute(f"SELECT * FROM gene_data WHERE POS = '{search_value}'").fetchall()
            return result
        else:
            return "Enter a valid position range" #invalid input

def search_db(search_type, search_value):
    if search_type == "snp_id":
        return use_id(search_value)
    elif search_type == "gene_name":
        return use_gene(search_value)
    elif search_type == "pos_numbers":
        return use_pos(search_value)
    else:
        return "404: Not Found"



################ population selection ################

# col_name = "FK_ID,POS,ID,REF,ALT,IMPACT,GENE,GENEID,AF_ALT_GBR,AF_ALT_CDX,AF_ALF_MXL"
col_name = ""
with sqlite3.connect('test.db') as connection:
    cursor = connection.cursor()
    col_data = cursor.execute("PRAGMA table_info(gene_data);").fetchall()
for value in col_data:
    col_name += value[1] + ','
col_name = col_name[:-1]
col_list = col_name.split(',')

# print(col_list)

def select_pop(pop_list, data):
    data = pd.DataFrame(data, columns=col_list)
    filtered_data = data.iloc[:,:8]
    for pop in pop_list:
        if pop == 'gbr':
            pop_cols = [col for col in data.columns if 'GBR' in col]
            filtered_data = filtered_data.join(data[pop_cols])
        elif pop == 'lwk':
            pop_cols = [col for col in data.columns if 'LWK' in col]
            filtered_data = filtered_data.join(data[pop_cols])
        elif pop == 'mxl':
            pop_cols = [col for col in data.columns if 'MXL' in col]
            filtered_data = filtered_data.join(data[pop_cols])
        elif pop == 'cdx':
            pop_cols = [col for col in data.columns if 'CDX' in col]
            filtered_data = filtered_data.join(data[pop_cols])
        elif pop == 'gih':
            pop_cols = [col for col in data.columns if 'GIH' in col]
            filtered_data = filtered_data.join(data[pop_cols])
    return filtered_data