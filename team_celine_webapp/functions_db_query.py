import sqlite3, re, io , base64
import pandas as pd
import numpy as np
import allel as sc
import matplotlib.pyplot as plt
import seaborn as sb
from markupsafe import Markup
from math import ceil
from statistics import mean, median


#################################################
################ QUERY FUNCTIONS ################
# make database connection to flask using sqlite3.connect and connection cursor

# user enter snp id (rs)
def use_id(search_value):
    with sqlite3.connect("chromosome22_snps.db") as connection:
        search_value = search_value.lower()
        cursor = connection.cursor()
        result = cursor.execute(f"SELECT * FROM snp_table WHERE rsID = '{search_value}'").fetchall()
    return result

# user enter gene name or gene alias
def use_gene(search_value):
    with sqlite3.connect("chromosome22_snps.db") as connection:
        search_value = search_value.upper()
        cursor = connection.cursor()
        result = cursor.execute(f"""SELECT * FROM snp_table WHERE GENE = '{search_value}'
                                    OR GENE_ALIAS LIKE '%{search_value},%'
                                    """).fetchall()
    return result

# user enter genomic position
# 2 types of input: ranged number separated by hyphen or a position number
def use_pos(search_value):
    range = r"(^\d+)\-(\d+$)"
    position = r"(^\d+$)"
    with sqlite3.connect("chromosome22_snps.db") as connection:
        cursor = connection.cursor()
        if re.search(range, f"{search_value}") is not None:
            input_nums = (str(search_value)).split('-')
            result = cursor.execute(f"SELECT * FROM snp_table WHERE POS BETWEEN '{input_nums[0]}' AND '{input_nums[1]}'").fetchall()
            return result
        elif re.search(position, f"{search_value}") is not None:
            result = cursor.execute(f"SELECT * FROM snp_table WHERE POS = '{search_value}'").fetchall()
            return result

# search type: radio buttons
# search value: text box value from html form
def search_db(search_type, search_value):
    if search_type == "snp_id":
        return use_id(search_value)
    elif search_type == "gene_name":
        return use_gene(search_value)
    elif search_type == "position":
        return use_pos(search_value)


############################################
################ STATISTICS DF ################

# all statistics need dataframe as input
def to_df(data): # data is list of tuples
    # get column names from database 
    col_name = ""
    with sqlite3.connect ('chromosome22_snps.db') as connection:
        cursor = connection.cursor()
        col_table = cursor.execute("PRAGMA table_info(snp_table);").fetchall()
    for value in col_table:
        col_name += value[1] + ','
    col_name = col_name[:-1] # remove the last comma
    col_list = col_name.split(',') # put col names in a list
    # build the dataframe
    data_df = pd.DataFrame(data, columns=col_list)
    return data_df

# make a new dataframe for stats
def calc_stats(df, stats_list, pop_list):
    stats_df = pd.DataFrame()
    for pop in pop_list:
        for stats in stats_list:
            if stats == 'shannon':
                shannon = windowedShannonDiv(df, pop)
                stats_df[pop+'_shannon']=shannon
            elif stats == 'tajima':
                tajima = windowedTajimasD(df, pop)
                stats_df[pop+'_tajima']=tajima
            elif stats == 'hetero':
                hetero = windowedHetDiv(df, pop)
                stats_df[pop+"_hetero"] = hetero
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


# heterozygosity per SNP
def calcHeterozygosity(df, pop):
    AF_pop = [col for col in df.columns if pop in col] # filter population
    AF_ref = [col for col in AF_pop if 'AF_ref' in col] # ref allele frequency
    AF_alt = [col for col in AF_pop if 'AF_alt' in col] # alt allele frequency
    AF_ref = ','.join(AF_ref)
    AF_alt = ','.join(AF_alt) 
    df[AF_ref] = df[AF_ref].astype(str).astype(float) # column values type as float
    df[AF_alt] = df[AF_alt].astype(str).astype(float)
    het = 1-((df[AF_alt])**2 + (df[AF_ref])**2) # calculate exp. het. according to formula (always only 2 alleles)
    return(het)


################################################
################ SLIDING WINDOW ################

def windowedShannonDiv(df, pop, w = None):
    AF_pop = [col for col in df.columns if pop in col] #filter population
    AF_col = [col for col in AF_pop if 'AF_alt' in col] #alt allele frequency
    df[AF_col] = df[AF_col].astype(str).astype(float) # column values type as float
    n = df.shape[0] # nrows of data frame
    if (w == None):
        w = ceil(n/10)
    shannon_divs_per_w = []
    df['norm_alt'] = df[AF_col]/n # normalize alt. counts
    df['ln_alt'] = np.log(df['norm_alt']) # take natural log (ln) of alt. counts
    df['shannon'] = df['norm_alt'] * df['ln_alt'] # multiply
    df['shannon'] = df['shannon'].fillna(0) # for cases ln(0) = -inf, set the multiplication as 0
    for i in range(0,n-w+1): 
        shannon_i = (df['shannon'][i:(i+w)]).sum() * -1 # calculate shannon diverseity per window
        shannon_divs_per_w.append(shannon_i) # store those values per window in a list/vector
    df.drop(['norm_alt', 'ln_alt'], axis=1, inplace=True)
    return(shannon_divs_per_w)

def windowedTajimasD(df, pop, w = None):
    AC_pop = [col for col in df.columns if pop in col] # filter population
    AC_ref = [col for col in AC_pop if 'AC_ref' in col] # ref allele count
    AC_alt = [col for col in AC_pop if 'AC_alt' in col] # alt allele count
    AC_ref = ','.join(AC_ref)
    AC_alt = ','.join(AC_alt)
    df[AC_ref] = df[AC_ref].astype(str).astype(int) #tajima needs integer
    df[AC_alt] = df[AC_alt].astype(str).astype(int)
    n = df.shape[0]
    if (w == None):
        w = ceil(n/10)
    wind_tajd = []
    arr = (df[[AC_ref, AC_alt]]).to_numpy() # take allele counts and save as numpy array
    for i in range(0, (n-w+1)):
        ac_window = arr[i:i+w] # set window of allele counts array
        D_i = sc.allel.tajima_d(ac_window) #calculate tajima's D using scikit-allel function
        wind_tajd.append(D_i)
    # change na to 0
    wind_tajd = np.nan_to_num(wind_tajd)
    return(wind_tajd)

def windowedHetDiv(df, pop, w = None):
    het_divs_per_w = []
    n = df.shape[0]
    if (w == None):
        w = ceil(n/10)
    het_vec = calcHeterozygosity(df, pop)
    for i in range(0,n-w+1):
        ### mean or median ? 
        het_i = mean(het_vec[i:(i+w)]) # take mean of the exp. het. for chosen window
        het_divs_per_w.append(het_i)
    return het_divs_per_w


####################################################
################ PLOT DISTRIBUTIONS ################

def shannon_plot(stats_df, pop_list):
    plt.figure(figsize=(16, 5))
    # loop through all selected populations to plot as one figure
    for pop in pop_list:
        shan_list = stats_df[pop+'_shannon'].tolist()
        x = stats_df['positions'].astype(str).astype(int) # x axis is genomic coordinates
        plt.plot(x, shan_list, markersize = 1, label=pop)
        plt.legend()
    plt.title('Sliding Window - Shannon Diversity Index')
    plt.ylabel('Shannon Diversity')
    plt.xlabel('Genomic Coordinates')
    # x-axis will not be shown in exponential form
    ax=plt.gca()
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    # format x-axis ticks and labels
    n = ceil(stats_df.shape[0])
    plt.xticks(np.arange(min(x), max(x), n/0.5), rotation='vertical')
    # encode using base64 to be passed on to flask html
    img = io.BytesIO()
    plt.savefig(img, format='png', bbox_inches='tight')
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    plt.close('all')
    return plot_url

def tajima_plot(stats_df, pop_list):
    plt.figure(figsize=(16, 5))
    # loop through all selected populations to plot as one figure
    for pop in pop_list:
        taj_list = stats_df[pop+'_tajima'].tolist()
        x = stats_df['positions'].astype(str).astype(int) # x axis is genomic coordinates
        plt.plot(x,taj_list, markersize = 1, label=pop)
        plt.legend()
    plt.title("Sliding Window - Tajima's D")
    plt.ylabel("Tajima's D")
    plt.xlabel('Genomic Coordinates')
    # x-axis will not be shown in exponential form
    ax=plt.gca()
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    # format x-axis ticks and labels
    n = ceil(stats_df.shape[0])
    plt.xticks(np.arange(min(x), max(x), n/0.5), rotation='vertical')
    # encode using base64 to be passed on to flask html
    img = io.BytesIO()
    plt.savefig(img, format='png', bbox_inches='tight')
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    plt.close('all')
    return plot_url

def hetero_plot(stats_df, pop_list):
    plt.figure(figsize=(16, 5))
    # loop through all selected populations to plot as one figure
    for pop in pop_list:
        het_list = stats_df[pop+'_hetero'].tolist()
        x = stats_df['positions'].astype(str).astype(int) # x axis is genomic coordinates
        plt.plot(x, het_list, markersize = 1, label=pop)
        plt.legend()
    plt.title('Sliding Window - Heterozygosity Diversity Index')
    plt.ylabel('Heteozygosity Diversity')
    plt.xlabel('Genomic Coordinates')
    # x-axis will not be shown in exponential form
    ax=plt.gca()
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    # format x-axis ticks and labels
    n = ceil(stats_df.shape[0])
    plt.xticks(np.arange(min(x), max(x), n/0.5), rotation='vertical')
    # encode using base64 to be passed on to flask html
    img = io.BytesIO()
    plt.savefig(img, format='png', bbox_inches='tight')
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    plt.close('all')
    return plot_url

# get_plot function that will be called in flask route script
def get_plot(stats_df, pop_list, stats):
    if stats == 'shannon':
        return shannon_plot(stats_df, pop_list)
    elif stats == 'tajima':
        return tajima_plot(stats_df, pop_list)
    elif stats == 'hetero':
        return hetero_plot(stats_df, pop_list)

# save summary statistics plots into a variable
def summary_stats_plot (stats_df, stats_list, pop_list):
    all_plots = []
    for stats in stats_list:
        plot_url = get_plot(stats_df, pop_list, stats)
        plot = Markup('<img src="data:image/png;base64,{}" width: 800px; height: 300px>'.format(plot_url))
        all_plots.append(plot)
    return all_plots


#####################################
################ FST ################

# make allele counts of chose populations as array
def extract_and_makearray(df, pop):
    # filter allele count columns of chosen populations
    AC_pop = [col for col in df.columns if pop in col]
    AC_ref = [col for col in AC_pop if 'AC_ref' in col]
    AC_alt = [col for col in AC_pop if 'AC_alt' in col]
    AC_ref = ','.join(AC_ref)
    AC_alt = ','.join(AC_alt)
    # columns values need to be integers
    df[AC_ref] = df[AC_ref].astype(str).astype(int)
    df[AC_alt] = df[AC_alt].astype(str).astype(int)
    # get ref and alt columns
    allele_count = df[[AC_ref, AC_alt]]
    # get number of variants 
    total_variants = allele_count.shape[0]
    allelecount_array = allele_count.to_numpy() # make into array
    reshape_alelle = allelecount_array.reshape(total_variants,2)
    return reshape_alelle

# calculate fst and plot heatmap
def calcFst(df, pop_list):
    # make a dataframe matrix
    fst_df = pd.DataFrame(index=pop_list, columns=pop_list)
    # pass array for all possible combinations of 2 populations
    for i in range(0, len(pop_list)):
        for j in range(i, len(pop_list)):
            pop1_array = extract_and_makearray(df, pop_list[i])
            pop2_array = extract_and_makearray(df, pop_list[j])
            pop1,pop2= sc.hudson_fst(pop1_array, pop2_array)
            fst = np.sum(pop1) / np.sum(pop2)
            fst_df.iat[i,j] = fst
            fst_df.iat[j,i] = fst
    fst_df = fst_df.apply(pd.to_numeric)
    # plot heatmap
    fig = sb.heatmap(fst_df, cmap = 'crest', annot = True)
    heatmap = fig.get_figure()
    plt.close('all')
    # encode using base64 to be passed on to flask html
    img = io.BytesIO()
    heatmap.savefig(img, format='png')
    img.seek(0)
    heatmap_url = base64.b64encode(img.getvalue()).decode()
    heatmap = Markup('<img src="data:image/png;base64,{}" width: 400; height: 300>'.format(heatmap_url))
    return (fst_df, heatmap)