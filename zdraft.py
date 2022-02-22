from cgi import test
# from curses import ACS_RTEE
from weakref import ref
import pandas as pd
import numpy as np
import allel as sc
import sqlite3, csv, re
from scipy.stats import entropy
from sqlalchemy import false
import functions_db_query as dbq
from statistics import mean
from math import floor, ceil, e
import matplotlib.pyplot as plt 
# import sliding_window_all_summary_stats as stats


df = pd.read_csv('final_complete_SNP_table.csv', sep=',')
df = df.head(100)
stats_list = ['shannon', 'tajima','hetero']
pop_list = ['GBR', 'CDX']


# if 'shannon' or 'tajima' in stats_list:
#     shan_taj_df = pd.DataFrame()
#     print(shan_taj_df)
# elif 'hetero' in stats_list:
#     het_df = pd.DataFrame()
#     print(het_df)


stats_df = pd.DataFrame()
for pop in pop_list:
    for stats in stats_list:
        if stats == 'shannon':
            shannon = dbq.windowedShannonDiv(df, pop)
            stats_df[pop+'_shannon']=shannon
        elif stats == 'tajima':
            tajima = dbq.windowedTajimasD(df, pop)
            stats_df[pop+'_tajima']=tajima
        elif stats == 'hetero':
            hetero = dbq.calcHeterozygosity(df, pop)
            heterow = dbq.windowedHetDiv(df, pop)
            stats_df[pop+"_hetero"] = heterow
print(stats_df)
            
# if 'shannon' in stats_list:
#     shannon_df = pd.DataFrame({'pop' : pop_list}) # make df for shannon value per population
#     shannon_list = []
#     for pop in pop_list: # loop to get shannon per pop
#         shannon = dbq.calcShannonDiv(data_df, pop)
#         shannon_list.append(shannon)
#     shannon_df["shannon"] = shannon_list # input shannon values into df
#     shannon = shannon_df.to_numpy() # df to list of tuples for html
# elif 'tajima' in stats_list:
#     tajima_df = pd.DataFrame({'pop' : pop_list})
#     tajima_list = []
#     for pop in pop_list:
#         tajima = dbq.calcTajimaD(data_df, pop)
#         tajima_list.append(tajima)
#     tajima_df["tajima"] = tajima_list
#     tajima = tajima_df.to_numpy()
# elif 'hetero' in stats_list:
#     hetero_df = pd.DataFrame({'pop' : pop_list})
#     hetero_list = []
#     for pop in pop_list:
#         hetero = dbq.calcHeterozygosity(data_df, pop)
#         hetero_list.append(hetero)
#     hetero_df["heterozygosity"] = hetero_list
#     hetero = hetero_df.to_numpy()

# else:
#     shannon="Choose statistics to calculate"


# shannondiv = dbq.windowedShannonDiv(df, pop)

# plt.plot(shannondiv,  'r', markersize = 1)
# plt.title('Sliding Window - Shannon Diversity Index')
# plt.ylabel('Shannon Diversity')
# plt.xlabel('Windows')
# plt.show()






# tajima needs ref and alt
# if 'tajima' in stats_list:
#     AC_pop = [col for col in df.columns if pop in col]
#     AC_ref = [col for col in AC_pop if 'AC_ref' in col]
#     AC_alt = [col for col in AC_pop if 'AC_alt' in col]
#     AC_ref = ','.join(AC_ref)
#     AC_alt = ','.join(AC_alt)
#     df[AC_ref] = df[AC_ref].astype(str).astype(float)
#     df[AC_alt] = df[AC_alt].astype(str).astype(float)
#     print(df[[AC_ref, AC_alt]])
#     arr = (df[[AC_ref, AC_alt]]).to_numpy() # take allele counts and save as numpy array
#     taj = sc.allel.tajima_d(arr) # calculate tajima's D using scikit-allel function
# print(taj)


# pop_list = ['gbr', 'chx']
# col_name = ""
# with sqlite3.connect('test.db') as connection:
#     cursor = connection.cursor()
#     col_data = cursor.execute("PRAGMA table_info(gene_data);").fetchall()
# for value in col_data:
#     col_name += value[1] + ','
# col_name = col_name[:-1]
# col_list = col_name.split(',')

# # print(col_list)

# data = pd.DataFrame(data, columns=col_list)
# filtered_data = data.iloc[:,:8]
# for pop in pop_list:
#     if pop == 'gbr':
#         gbr_cols = [col for col in data.columns if 'GBR' in col]
#         filtered_data = filtered_data.join(data[gbr_cols])
#     elif pop == 'lwk':
#         filtered_data = filtered_data.join(data.iloc[:,[8]])
#     elif pop == 'mxl':
#         filtered_data = filtered_data.join(data.iloc[:,[9]])
#     elif pop == 'cdx':
#         filtered_data = filtered_data.join(data.iloc[:,[10]])
#     elif pop == 'gih':
#         filtered_data = filtered_data.join(data.iloc[:,[10]])

# print(filtered_data)



# base_data = data.iloc[:,:8]
# for pop in population:
#     if pop == 'gbr':
#         base_data = base_data.join(data.iloc[:,[8]])
#     elif pop == 'chx':
#         base_data = base_data.join(data.iloc[:,[9]])
#     elif pop == 'mxn':
#         base_data = base_data.join(data.iloc[:,[10]])
# print(base_data)





# def use_pos(pos_numbers):
#     # 2 types of input: ranged number separated by hyphen or a position number
#     range = r"(^\d+)\-(\d+$)"
#     position = r"(^\d+$)"
#     with sqlite3.connect("chr22.db") as connection:
#         cursor = connection.cursor()
#         if re.search(range, f"{pos_numbers}") is not None:
#             input_nums = (str(pos_numbers)).split('-')
#             result = cursor.execute(f"SELECT * FROM data WHERE POS BETWEEN '{input_nums[0]}' AND '{input_nums[1]}'").fetchall()
#             return result
#         elif re.search(position, f"{pos_numbers}") is not None:
#             result = cursor.execute(f"SELECT * FROM data WHERE POS = '{pos_numbers}'").fetchall()
#             return result
#         else:
#             return('Enter valid position range') #invalid input
# print(use_pos("10519276-10519413"))

# ^[0-9]+\-[0-9]+$
# ^(\d+-?)+\d+$
# r"([A-Za-z0-9\.]+)@[A-Za-z0-9\.]*qmul\.ac\.uk"

# # SELECT STATEMENT FOR POSITION
# with sqlite3.connect("chr22.db") as connection:
#         cursor = connection.cursor()
#         result = cursor.execute(f"SELECT * FROM table WHERE POS BETWEEN '{start_pos}' AND '{end_pos}'").fetchall()

# print (result)

# def pos_range(pos_numbers):
#     pattern = r"(^\d+)\-(\d+$)"
#     if re.search(pattern, f"{pos_numbers}"):
#         return 'Found a match!'
#     else:
#         return('Not matched!')


# shannon using genotype
# df = pd.read_csv('mockdiversity.csv', sep=',')
# df = df.iloc[:,5:]
# df2 = list(df.loc[1])
# pd_series = pd.Series(df2)
# counts = pd_series.value_counts()
# entropy = entropy(counts)
# print(entropy)

# data.values.tolist() - for all rows

# pd_series = pd.Series(data)
# counts = pd_series.value_counts()
# entropy = entropy(counts)
# print(entropy)

# # rearrange column
# id = gt5['ID']
# gt5.drop(labels=['ID'], axis=1,inplace = True)
# gt5.insert(0, 'PK_ID', id)
# # to csv
# gt5.to_csv('22_GBR_genotype.csv', index=False)





############
# con = sqlite3.connect("chromosome22_snps.db")
# cursor = con.cursor()

# result = cursor.execute("""SELECT * FROM snp_table WHERE rsID = 'NCAPH2'
#                             UNION
#                             SELECT * FROM snp_table WHERE GENE_ALIAS = 'U2'
#                         """).fetchall()
# print(result)

# con.commit()
# con.close()
#############




# with open ('population_samples.csv', 'r') as i:
#     reader = csv.reader(i)
#     columns = next(reader) 
#     query = 'INSERT INTO pop_samples ({0}) VALUES ({1})'
#     query = query.format(','.join(columns), ','.join('?' * len(columns)))
#     for data in reader:
#         cursor.execute(query, data)

# check duplicate rows by column value
# df = pd.read_csv("22_GBR_genotype.csv")
# boo = df["ID"].duplicated().any() #false means no dup
# uni = df['ID'].unique()
# print(boo)

# df = pd.read_csv("22_gene_names.csv")
# print(df.shape)

# def db_connection():
#     conn = sqlite3.connect('chr22.db')
#     conn.row_factory = sqlite3.Row
#     return conn
# #get CHROM, ID, POS, REF, ALT from vcf file and store in csv file
# sc.vcf_to_csv('genotype_population.vcf', 'basic_vcf_info.csv', fields=['CHROM', 'POS', 'ID', 'REF', 'ALT'])

# # get genotype data of GBR samples
# gt_df = pd.read_csv('GBR_data.csv', sep="\t")
# gt_df2 = gt_df.iloc[: ,9:]
# # combine csv files
# df1 = pd.read_csv('basic_vcf_info.csv')
# (pd.concat([df1, gt_df2], axis=1).to_csv('GBR_genotype.csv', index=False, na_rep='N/A'))


# # extract the 11 relevant columns incl. gene_names from annotated txt file
# gene_names = pd.read_csv('gene_names.txt', sep="\t", header=None)
# first_n_columns = gene_names.iloc[:, :11]
# # print(first_n_columns.head(8))
# first_n_columns.to_csv('gene_names.csv',sep = ',', index = None, header = False)


# # store population info in csv file 
# pop_df = pd.read_csv("population_samples.tsv", sep = "\t")
# pop_df.to_csv("population_samples.csv", index = None)


# # df_GBR = pop_df.loc[pop_df['Population code'] == "GBR"]


# ######

# con = sqlite3.connect("chr22.db")
# cur = con.cursor()

# cur.execute("""SELECT * FROM GBR_genotype WHERE ALT_2 IS NULL """)
# print(cur.fetchall())

# con.commit()
# con.close()


# population_sample table
# pop = pd.read_csv("population_samples.csv")
# pop.to_sql("pop_samples", con, index = False)

# cur.execute("CREATE TABLE pop_sample (Sample_name,Sex,BiosampleID,Population_code,Population_name,Superpopulation_code,Superpopulation_name,Population_elastic_ID,Data_collections);") # use your column names here

# with open('population_samples.csv','r') as fin:
#     dr = csv.DictReader(fin) # comma is default delimiter
#     to_db = [(i['Sample name'], i['Sex'], i['Biosample ID'], i['Population code'], i['Population name'], i['Superpopulation code'], i['Superpopulation name'], i['Population elastic ID'], i['Data collections']) for i in dr]
# cur.executemany("INSERT INTO pop_sample (Sample_name,Sex,BiosampleID,Population_code,Population_name,Superpopulation_code,Superpopulation_name,Population_elastic_ID,Data_collections) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);", to_db)


# gene_names table
# cur.execute("CREATE TABLE gene_names (CHROM, POS, ID, REF, ALT, FILTER, ALLELE, EFFECT, IMPACT, GENE, GENEID);") # use your column names here

# with open('gene_names.csv','r') as fin:
#     dr = csv.DictReader(fin) # comma is default delimiter
#     to_db = [(i['CHROM'], i['POS'], i['ID'], i['REF'], i['ALT'], i['FILTER'], i['ALLELE'], i['EFFECT'], i['IMPACT'], i['GENE'], i['GENEID']) for i in dr]
# cur.executemany("INSERT INTO gene_names (CHROM, POS, ID, REF, ALT, FILTER, ALLELE, EFFECT, IMPACT, GENE, GENEID) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", to_db)



# # genotype table
# GBR_gt = pd.read_csv("GBR_genotype.csv")
# GBR_gt.to_sql('GBR_genotype',con, index=False)


# cur.execute("""SELECT * FROM gene_names INNER JOIN basic_info ON gene_names.POS = basic_info.POS
#             WHERE gene_names.POS = '10519410'  """)
# print(cur.fetchall())

# cur.execute("""SELECT * FROM GBR_genotype WHERE POS  = '50807929'""")
# print(cur.fetchall())

# cur.execute("DROP TABLE GBR_genotype")

