from os import stat
import sqlite3
import csv
import pandas as pd
import numpy as np
import re
import io
import base64
import allel as sc
from math import ceil
from statistics import mean, median
import matplotlib.pyplot as plt


# change genotype symbols to alphabet
def gt_to_alp(data, start_col, end_col):
    for pos in range(start_col, end_col):
        data.iloc[:, pos] = np.where(
            data.iloc[:, pos] == '0|1', data['REF'] + data['ALT'], data.iloc[:, pos])
        data.iloc[:, pos] = np.where(
            data.iloc[:, pos] == '1|0', data['ALT'] + data['REF'], data.iloc[:, pos])
        data.iloc[:, pos] = np.where(
            data.iloc[:, pos] == '1|1', data['ALT'] + data['ALT'], data.iloc[:, pos])
        data.iloc[:, pos] = np.where(
            data.iloc[:, pos] == '0|0', data['REF'] + data['REF'], data.iloc[:, pos])
    return data

# filter out multi-allelic


def single_alt(df):
    df = df[df['ALT_2'].isna()]  # extract rows with 1 ALT
    df.drop(['ALT_2', 'ALT_3'], axis=1, inplace=True)
    df.rename(columns={"ALT_1": "ALT"}, inplace=True)
    return df

# filter out irrelevant gene IDs


def gene_id_refinement(df):
    # take out those Gene IDs of the sorts "ENS\d+-ENS\d+""
    pat = r'(\w+)\-(\w+)'
    df2 = df[df.GENEID.str.contains(pat) == False]
    return df2

# extract SNPs and biallelic


def snps_only(data):
    data2 = data['ALT'].str.split(expand=True).apply(
        lambda x: x.str.len()) <= 1
    data3 = data[data2.any(axis=1)]
    data4 = data3['REF'].str.split(
        expand=True).apply(lambda x: x.str.len()) <= 1
    data5 = data3[data4.any(axis=1)]
    return data5

# insert values into table in db


def insert_table(csv_file, database_name, table_name):
    with sqlite3.connect(database_name) as con:
        cursor = con.cursor()
        with open(csv_file, 'r') as i:
            reader = csv.reader(i)
            columns = next(reader)
            query = "INSERT INTO " + table_name + " ({0}) VALUES ({1})"
            query = query.format(
                ','.join(columns), ','.join('?' * len(columns)))
            for data in reader:
                cursor.execute(query, data)
            con.commit()


############### query functions ################

def use_id(search_value):
    with sqlite3.connect("chromosome22_snps.db") as connection:
        search_value = search_value.lower()
        cursor = connection.cursor()
        result = cursor.execute(
            f"SELECT * FROM snp_table WHERE rsID = '{search_value}'").fetchall()
    return result


def use_gene(search_value):
    with sqlite3.connect("chromosome22_snps.db") as connection:
        search_value = search_value.upper()
        cursor = connection.cursor()
        result = cursor.execute(f"""SELECT * FROM snp_table WHERE GENE = '{search_value}'
                                    OR GENE_ALIAS LIKE '%{search_value}%'
                                    """).fetchall()
    return result


def use_pos(search_value):
    # 2 types of input: ranged number separated by hyphen or a position number
    range = r"(^\d+)\-(\d+$)"
    position = r"(^\d+$)"
    with sqlite3.connect("chromosome22_snps.db") as connection:
        cursor = connection.cursor()
        if re.search(range, f"{search_value}") is not None:
            input_nums = (str(search_value)).split('-')
            result = cursor.execute(
                f"SELECT * FROM snp_table WHERE POS BETWEEN '{input_nums[0]}' AND '{input_nums[1]}'").fetchall()
            return result
        elif re.search(position, f"{search_value}") is not None:
            result = cursor.execute(
                f"SELECT * FROM snp_table WHERE POS = '{search_value}'").fetchall()
            return result
        else:
            return "Enter a valid position range"  # invalid input


def search_db(search_type, search_value):
    if search_type == "snp_id":
        return use_id(search_value)
    elif search_type == "gene_name":
        return use_gene(search_value)
    elif search_type == "position":
        return use_pos(search_value)
    else:
        return "404: Not Found"


################ statistics ################

# all statistics need dataframe as input
def to_df(data):  # data is list of tuples
    # get column names from database
    col_name = ""
    with sqlite3.connect('chromosome22_snps.db') as connection:
        cursor = connection.cursor()
        col_table = cursor.execute("PRAGMA table_info(snp_table);").fetchall()
    for value in col_table:
        col_name += value[1] + ','
    col_name = col_name[:-1]  # remove the last comma
    col_list = col_name.split(',')  # put col names in a list
    data = pd.DataFrame(data, columns=col_list)
    return data

# make a new dataframe for stats


def calc_stats(df, stats_list, pop_list):
    stats_df = pd.DataFrame()
    for pop in pop_list:
        for stats in stats_list:
            if stats == 'shannon':
                shannon = windowedShannonDiv(df, pop)
                stats_df[pop+'_shannon'] = shannon
            elif stats == 'tajima':
                tajima = windowedTajimasD(df, pop)
                stats_df[pop+'_tajima'] = tajima
            elif stats == 'hetero':
                heterow = windowedHetDiv(df, pop)
                stats_df[pop+"_hetero"] = heterow
    # GRACIA
    # add positions for x-axis of distribution plot 
    pos_per_w = []
    n = df.shape[0]
    w = ceil(n/10)
    for i in range(0, n-w+1):
        # take median position of the window
        pos_i = median(df['POS'][i:(i+w)]) # maybe need to round up? 
        pos_per_w.append(pos_i)
    stats_df['positions'] = pos_per_w
    return stats_df

# shannon
def calcShannonDiv(df, pop):
    AF_pop = [col for col in df.columns if pop in col]
    AF_col = [col for col in AF_pop if 'AF_alt' in col]
    n = df.shape[0]  # nrows of data frame
    df[AF_col] = df[AF_col].astype(str).astype(float)  # column type as float
    df['norm_alt'] = df[AF_col]/n  # normalize alt. counts
    # take natural log (ln) of alt. counts
    df['ln_alt'] = np.log(df['norm_alt'])
    df['shannon'] = df['norm_alt'] * df['ln_alt']  # multiply
    # for cases ln(0) = -inf, set the multiplication as 0
    df['shannon'] = df['shannon'].fillna(0)
    # calculate sum and multiply with -1 to get shannon diversity
    shannondiv = df['shannon'].sum() * -1
    df.drop(['norm_alt', 'ln_alt'], axis=1, inplace=True)
    return(shannondiv)
# Tajima's D


def calcTajimaD(df, pop):
    AC_pop = [col for col in df.columns if pop in col]  # filter population
    AC_ref = [col for col in AC_pop if 'AC_ref' in col]  # ref allele count
    AC_alt = [col for col in AC_pop if 'AC_alt' in col]  # alt allele count
    AC_ref = ','.join(AC_ref)
    AC_alt = ','.join(AC_alt)
    df[AC_ref] = df[AC_ref].astype(str).astype(int)
    df[AC_alt] = df[AC_alt].astype(str).astype(int)
    # take allele counts and save as numpy array
    arr = (df[[AC_ref, AC_alt]]).to_numpy()
    # calculate tajima's D using scikit-allel function
    taj = sc.allel.tajima_d(arr)
    return (taj)
# heterozygosity


def calcHeterozygosity(df, pop):
    AF_pop = [col for col in df.columns if pop in col]  # filter population
    AF_ref = [col for col in AF_pop if 'AF_ref' in col]  # ref allele freq
    AF_alt = [col for col in AF_pop if 'AF_alt' in col]
    AF_ref = ','.join(AF_ref)
    AF_alt = ','.join(AF_alt)
    df[AF_ref] = df[AF_ref].astype(str).astype(float)
    df[AF_alt] = df[AF_alt].astype(str).astype(float)
    # calculate exp. het. according to formula (always only 2 alleles)
    het = 1-((df[AF_alt])**2 + (df[AF_ref])**2)
    return(het)


# sliding window

def windowedShannonDiv(df, pop, w=None):
    AF_pop = [col for col in df.columns if pop in col]
    AF_col = [col for col in AF_pop if 'AF_alt' in col]
    df[AF_col] = df[AF_col].astype(str).astype(float)  # column type as float
    n = df.shape[0]  # nrows of data frame
    if (w == None):
        w = ceil(n/10)
    shannon_divs_per_w = []
    df['norm_alt'] = df[AF_col]/n  # normalize alt. counts
    # take natural log (ln) of alt. counts
    df['ln_alt'] = np.log(df['norm_alt'])
    df['shannon'] = df['norm_alt'] * df['ln_alt']  # multiply
    # for cases ln(0) = -inf, set the multiplication as 0
    df['shannon'] = df['shannon'].fillna(0)
    for i in range(0, n-w+1):
        # calculate shannon diverseity per window
        shannon_i = (df['shannon'][i:(i+w)]).sum() * -1
        # store those values per window in a list/vector
        shannon_divs_per_w.append(shannon_i)
    df.drop(['norm_alt', 'ln_alt'], axis=1, inplace=True)
    return(shannon_divs_per_w)


def windowedTajimasD(df, pop, w=None):
    AC_pop = [col for col in df.columns if pop in col]  # filter population
    AC_ref = [col for col in AC_pop if 'AC_ref' in col]  # ref allele count
    AC_alt = [col for col in AC_pop if 'AC_alt' in col]  # alt allele count
    AC_ref = ','.join(AC_ref)
    AC_alt = ','.join(AC_alt)
    df[AC_ref] = df[AC_ref].astype(str).astype(int)
    df[AC_alt] = df[AC_alt].astype(str).astype(int)
    n = df.shape[0]
    if (w == None):
        w = ceil(n/10)
    wind_tajd = []
    # take allele counts and save as numpy array
    arr = (df[[AC_ref, AC_alt]]).to_numpy()
    for i in range(0, (n-w+1)):
        ac_window = arr[i:i+w]  # set window of allele counts array
        D_i = sc.allel.tajima_d(ac_window)
        wind_tajd.append(D_i)
    wind_tajd = np.nan_to_num(wind_tajd) # GRACIA
    return(wind_tajd)


def windowedHetDiv(df, pop, w=None):
    het_divs_per_w = []
    n = df.shape[0]
    if (w == None):
        w = ceil(n/10)
    het_vec = calcHeterozygosity(df, pop)
    for i in range(0, n-w+1):
        # take mean of the exp. het. for chosen window
        het_i = mean(het_vec[i:(i+w)])
        het_divs_per_w.append(het_i)
    return het_divs_per_w


# plot distributions
# def get_plot(stats_df, pop, stats):
#     if stats == 'shannon':
#         return get_shannon_plot(stats_df, pop)
#     elif stats == 'tajima':
#         return get_tajima_plot(stats_df, pop)
#     elif stats == 'hetero':
#         return get_hetero_plot(stats_df, pop)

def get_plot(stats_df, pop_list, stats):
    if stats == 'shannon':
        return get_shannon_plot(stats_df, pop_list)
    elif stats == 'tajima':
        return get_tajima_plot(stats_df, pop_list)
    elif stats == 'hetero':
        return get_hetero_plot(stats_df, pop_list)


def get_shannon_plot(stats_df, pop_list):
    for pop in pop_list:
        shan_list = stats_df[pop+'_shannon'].tolist()
        plt.plot(shan_list,  'r', markersize=1)
        plt.title('Sliding Window - Shannon Diversity Index')
        plt.ylabel('Shannon Diversity')
        plt.xlabel('Windows')
        # encode
        img = io.BytesIO()
        plt.savefig(img, format='png')
        img.seek(0)
        plot_url = base64.b64encode(img.getvalue()).decode()
    plt.close('all')  # closing the plot window
    return plot_url

# def get_shannon_plot(stats_df, pop):
#     shan_list = stats_df[pop+'_shannon'].tolist()
#     plt.plot(shan_list,  'r', markersize = 1)
#     plt.title('Sliding Window - Shannon Diversity Index')
#     plt.ylabel('Shannon Diversity')
#     plt.xlabel('Windows')
#     # encode
#     img = io.BytesIO()
#     plt.savefig(img, format='png')
#     img.seek(0)
#     plot_url = base64.b64encode(img.getvalue()).decode()
#     plt.close('all') ### closing the plot window?
#     return plot_url


def get_tajima_plot(stats_df, pop):
    taj_list = stats_df[pop+'_tajima'].tolist()
    plt.plot(taj_list,  'g', markersize=1)
    plt.title("Sliding Window - Tajima's D")
    plt.ylabel("Tajima's D")
    plt.xlabel('Windows')
    # encode
    img = io.BytesIO()
    plt.savefig(img, format='png')
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    plt.close('all')  # closing the plot window?
    return plot_url


def get_hetero_plot(stats_df, pop):
    het_list = stats_df[pop+'_hetero'].tolist()
    plt.plot(het_list,  'b', markersize=1)
    plt.title('Sliding Window - Heterozygosity Diversity Index')
    plt.ylabel('Heteozygosity Diversity')
    plt.xlabel('Windows')
    # encode
    img = io.BytesIO()
    plt.savefig(img, format='png')
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    plt.close('all')  # closing the plot window?
    return plot_url


# fst
def extract_and_makearray(df, pop):
    AC_pop = [col for col in df.columns if pop in col]
    AC_ref = [col for col in AC_pop if 'AC_ref' in col]
    AC_alt = [col for col in AC_pop if 'AC_alt' in col]
    AC_ref = ','.join(AC_ref)
    AC_alt = ','.join(AC_alt)
    df[AC_ref] = df[AC_ref].astype(str).astype(int)
    df[AC_alt] = df[AC_alt].astype(str).astype(int)
    # get ref and alt
    allele_count = df[[AC_ref, AC_alt]]
    # get number of variants
    total_variants = allele_count.shape[0]
    allelecount_array = allele_count.to_numpy()
    # reshape to n_varaints, n_alleles
    reshape_alelle = allelecount_array.reshape(total_variants, 2)
    return reshape_alelle


def calcFst(df, pop_list):
    fst_dict = {}
  # pass array for 2 pops
    for i in range(0, len(pop_list)-1):
        for j in range(i+1, len(pop_list)):
            pop1_array = extract_and_makearray(df, pop_list[i])
            pop2_array = extract_and_makearray(df, pop_list[j])
            pop1, pop2 = sc.hudson_fst(pop1_array, pop2_array)
            fst = np.sum(pop1) / np.sum(pop2)
            key = pop_list[i] + "-" + pop_list[j]
            value = fst
            fst_dict[key] = value
    fst_df = pd.DataFrame.from_dict(fst_dict, orient='index', columns=['fst'])
    return (fst_df)


# def select_pop(pop_list, data):
#     data = pd.DataFrame(data, columns=col_list)
#     filtered_data = data.iloc[:,:8]
#     for pop in pop_list:
#         if pop == 'GBR':
#             pop_cols = [col for col in data.columns if 'GBR' in col]
#             filtered_data = filtered_data.join(data[pop_cols])
#         elif pop == 'LWK':
#             pop_cols = [col for col in data.columns if 'LWK' in col]
#             filtered_data = filtered_data.join(data[pop_cols])
#         elif pop == 'MXL':
#             pop_cols = [col for col in data.columns if 'MXL' in col]
#             filtered_data = filtered_data.join(data[pop_cols])
#         elif pop == 'CDX':
#             pop_cols = [col for col in data.columns if 'CDX' in col]
#             filtered_data = filtered_data.join(data[pop_cols])
#         elif pop == 'GIH':
#             pop_cols = [col for col in data.columns if 'GIH' in col]
#             filtered_data = filtered_data.join(data[pop_cols])
#     return filtered_data


#####
