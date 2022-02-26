import pandas as pd
import numpy as np
import allel
from matplotlib import pyplot as plt
import seaborn as sb


##### Read db as df #######
complete_df= pd.read_csv("/data/scratch/bt211065/2022_02_01_Groupproject/test2/SNP_data_incl_derived_AF.csv")


########################################################################
###############extract ref/alt count from df ###########################
########################################################################

# function accepts df and any number of key value arguments i.e. pop1 = 'MXL'
def extract_and_makelist(accept_dataframe, *args):
  # get ref and alt count and make df for each pop 
  allele_count_MXL = accept_dataframe.loc[:,'MXL_AC_ref': 'MXL_AC_alt']
  allele_count_GBR = accept_dataframe.loc[:,'GBR_AC_ref': 'GBR_AC_alt']
  allele_count_CDX = accept_dataframe.loc[:,'CDX_AC_ref': 'CDX_AC_alt']
  allele_count_LWK = accept_dataframe.loc[:,'LWK_AC_ref': 'LWK_AC_alt']
  allele_count_GIH = accept_dataframe.loc[:,'GIH_AC_ref': 'GIH_AC_alt']
  
  list_counts = []
  dict_pop = {}
  
  # check if gbr passed in key value args and grab pop counts to convert to numpt then to list 
  # convert to numpy/to list to store pop df as numpy array then convert to list with format as list of lists. Each element in list is list of AC ref/alt for each   snp [[], [], [] ]
  
  if 'GBR' in args:
    array_GBR=allele_count_GBR.to_numpy()
    GBR_list=array_GBR.tolist()
    dict_pop['GBR'] = GBR_list
  # store as dict- key is pop name and value list of counts 
  
  if 'MXL' in args:
    array_MXL=allele_count_MXL.to_numpy()
    MXL_list=array_MXL.tolist()
    dict_pop['MXL'] = MXL_list
    
  
  if 'CDX' in args:
    array_CDX= allele_count_CDX.to_numpy()
    CDX_list=array_CDX.tolist()
    dict_pop['CDX'] = CDX_list
    
  
  if 'LWK' in args:
    array_LWK=allele_count_LWK.to_numpy()
    LWK_list=array_LWK.tolist()
    dict_pop['LWK'] = LWK_list

    
  if 'GIH' in  args:
    array_GIH=allele_count_GIH.to_numpy()
    GIH_list=array_GIH.tolist()
    dict_pop['GIH'] = GIH_list
  
  
#  for k,v in dict_pop.items():
 #   print(len(v))
  #print(len(dict_pop.items()))
  return dict_pop
  
#extract_and_makelist(complete_df,'GBR',  'MXL' , 'CDX', 'LWK', 'GIH')

  
####################################################
######### fst calculation function ################
####################################################
 

def calculate_fst(accept_dataframe, *args):

  ########  call  extract function######### 
  dict_allelecounts = extract_and_makelist(accept_dataframe, *args)
  
  
  ##### make list of pop names from dict.keys() #########
  names = []
  list_names =dict_allelecounts.keys()
  
  for k in list_names:
    names.append(k)
  
   
  
  
  ########## duplicate list for each pop and  make dict each key coresponding to population ########
  
  dict_arraydups = {}
  
  
  for key, value in dict_allelecounts.items():
    if key == key:
      single_list = dict_allelecounts[key] # grab value of key 
      dict_arraydups[key]=[np.asarray(single_list) for x in range(len(names))] #make list of duplicated arrays + each key has value of duplicated list of arrays
      # len of values of each key equal length of keys/pops
  
  
      
  ###### make list of values in dict  and convert list to arrays for nondup#############
  
  
  list_values = dict_allelecounts.values()
  
  array_pop = []
  for i in list_values:
    array_pop.append((np.asarray(i))) ## convert each list in list values to array and append to list 
    
 
    
  ############# Calculat fst for each array in each dict value and store in dict ######

  fst_store={}
  fst_list = []
  
  for key,value in dict_arraydups.items(): # loop over key and value in dict with dupllicated arrays as values ( each value correspond to correct pop-key)
    if key == key: # each key and value pased for each elmenet in 3 lists
      for v,arr,n in zip(value,array_pop,names): 
      #dict value is looped over again as each value has list of arrays and each arrray(1array) needs to be passed into fstcal
         pop1,pop2= allel.hudson_fst(v ,arr)
         fst = np.sum(pop1)/np.sum(pop2)
         fst_list.append([fst,n]) # store n in names in fst list as correspond to nonduplicated array in arraypop passed - n is used as df row 
      fst_store[key]=fst_list # list with 2 val append after each inner loop to update dictionary - each key is used for col name in df 
      fst_list= [] # list empty before next inner loop
      

  return fst_store

#calculate_fst(complete_df, 'GBR', 'MXL', 'CDX', 'GIH', 'LWK')  


################################################################
########### making dataframe from fst val dictionary############ 
################################################################

def make_df(accept_dataframe, *args):
  fst_stores = calculate_fst(accept_dataframe, *args) # get fst store dict
  
  list_keys= fst_stores.keys() # get list of keys 
  
  
  
  fst_df = pd.DataFrame(index= list_keys , columns = list_keys) # make empty df  with index and column names as list of keys 
 
  
  
  ############# insert value in dict according to column/row in dataframe ###########
  i = 0 # set count to 0 - list index 
  fst_val = []
  for key, value in fst_stores.items():
    if key == key:
      for n in list_keys:
        col = key # assign key to col
        row = n # assign name of pop in list keys to row
        val = value[i][0] # extract every fst value in index 0 every nth sublist in list value
        i+=1 # increment by 1 for every inner loop
        fst_df.at[row,col] = val # insert fst value in cell by row and col
      i= 0 # set to 0 before next outer loop/key 
  
  
  return fst_df  
  
  
#make_df(complete_df, 'MXL', 'GBR', 'LWK', 'GIH', 'CDX')

#################################################################
######make heatplot from fst dataframe ##########################
#################################################################

def make_heatplot(accept_dataframe, *args):
  
  fst_df = make_df(complete_df,  *args) # get fst_df
  
   #### size of figure ####
  figure, ax = plt.subplots(figsize = (10,5))
  fst_df = fst_df.apply(pd.to_numeric) # convert each col in df to float from object for seaborn plot
  fig=sb.heatmap(fst_df, cmap="Reds", annot = True) 
  heatmap = fig.get_figure()
  heatmap.savefig('test.png')
  
make_heatplot(complete_df, 'MXL', 'CDX' , 'GBR', 'LWK')
  
  
  