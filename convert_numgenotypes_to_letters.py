import csv
from numpy import var 


def lines(file):
    header = '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00106	HG00113	HG00118	HG00120	HG00125	HG00132	HG00137	HG00149	HG00151	HG00233	HG00238	HG00240	HG00245	HG00252	HG00257	HG00264	HG00097	HG00100	HG00101	HG00105	HG00112	HG00117	HG00129	HG00131	HG00136	HG00143	HG00148	HG00150	HG00155	HG00232	HG00237	HG00244	HG00251	HG00256	HG00263	HG02215	HG00109	HG00111	HG00116	HG00123	HG00128	HG00130	HG00142	HG00154	HG00159	HG00157	HG00231	HG00236	HG00243	HG00234	HG00239	HG00246	HG00253	HG00258	HG00260	HG00265	HG00099	HG00102	HG00107	HG00114	HG00119	HG00121	HG00250	HG00126	HG00255	HG00133	HG00262	HG00138	HG00140	HG00145	HG00096	HG01789	HG01791	HG00158	HG00160	HG00235	HG00242	HG00254	HG00259	HG00261	HG00103	HG00108	HG00110	HG00115	HG00122	HG00127	HG00139	HG00141	HG00146	HG01334	HG01790'
    list_csv = []
    csv_file= open(file, 'r') 
    csvfile_read= csv.reader(csv_file)
        # csvfile_read is reader object
    for record in csvfile_read:
        list_csv.append(','.join(record))
        # get index of header in list_csv 
    header_list = header.split()
    #print(header_list)
    index_header= list_csv.index(header)
    
    # extract all records from header 
    record= list_csv[index_header:]
    return record


###########extracting fields############

import pandas 

def extract_fields(function,file):
    # call function read 10 lines 
    recs= lines(file) 
    header= '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00106	HG00113	HG00118	HG00120	HG00125	HG00132	HG00137	HG00149	HG00151	HG00233	HG00238	HG00240	HG00245	HG00252	HG00257	HG00264	HG00097	HG00100	HG00101	HG00105	HG00112	HG00117	HG00129	HG00131	HG00136	HG00143	HG00148	HG00150	HG00155	HG00232	HG00237	HG00244	HG00251	HG00256	HG00263	HG02215	HG00109	HG00111	HG00116	HG00123	HG00128	HG00130	HG00142	HG00154	HG00159	HG00157	HG00231	HG00236	HG00243	HG00234	HG00239	HG00246	HG00253	HG00258	HG00260	HG00265	HG00099	HG00102	HG00107	HG00114	HG00119	HG00121	HG00250	HG00126	HG00255	HG00133	HG00262	HG00138	HG00140	HG00145	HG00096	HG01789	HG01791	HG00158	HG00160	HG00235	HG00242	HG00254	HG00259	HG00261	HG00103	HG00108	HG00110	HG00115	HG00122	HG00127	HG00139	HG00141	HG00146	HG01334	HG01790'
    list_keys=header.split()


    # convert each record to list 
    list_recs= []
    for r in recs[1:]: 
        list_recs.append(r.split())
    print(list_recs)
    

     
    #grab ref alleles
    ref = []
    for rec in list_recs:
        ref.append(rec[3])
    
    variant_dict= {'Ref':ref}
    #print(variant_dict)
    
    # grab alt alleles 
    alt = []
    for al in list_recs:
        alt.append(al[3])
    
    #print(alt)
    

    variant_dict['Alt']= alt
    #print(variant_dict)

    #grab genotpye section for each list/record in a list 
    genotype = []
    for sl in list_recs:
        genotype.append(sl[9:])
    
    print(list_recs)
    #iterate throught each record in list  
    assigned_gt=[]
    for gt in genotype:
        # for each record zip sample names
        assigned_gt.append(tuple(zip(gt,list_keys[9:])))
    samples= list_keys[9:]
    #print(assigned_gt)

    ##### convert tuple assignt gt to list ####
    record_sample_gt = [] # store and append as record
    list_gt_sample_assigned = [] # store append as pairs 
    # store list 
    # iterated over each record in assigned gt 
    for record in assigned_gt:  
        # itearate each pair of record   
        for pair in record:
            list_gt_sample_assigned.append(list(pair))
        #convert each pair to list and append
        # at end of first inner loop append to different nelits and clear storage for pair to be used in next outer loop         
        record_sample_gt.append(list_gt_sample_assigned)
        list_gt_sample_assigned = []    
            
    #print(record_sample_gt)
    ############assign key as sample names with list of values correpsonding to samples/indviduals##############
    # should have 5 (test) dict items/keys,values 
    i=0
    n=0
    
    genotypedict=dict.fromkeys(samples,[])    
    
     #iteratcte over each record in list record_sample gt
   # for rec in record_sample_gt:
        # iteartef oover each gentoype assigend sample  
              
    for rec in record_sample_gt:
        for items in genotypedict.items():
            if items[0] == rec[i][1]:
                genotypedict.update({items[0]:rec[i][0]+'\t'})
                i+=1
        i=0
        for key in genotypedict.keys():
            if key == rec[n][1]:
                genotypedict[key]+=rec[n][0]
                n +=1
        
        n=0
    
    #dict2={key:list(value) for key,value in genotypedict.items()}  
    print(genotypedict)
    #print(dict2)
    #### merge dictionaries####
    
    
    #completedict = variant_dict.copy()
    #completedict.update(dict2)
    
    

    #return variant_dict, dict2



    
        
def joind(function,file):
    variant_dict , dict2 = extract_fields(lines,file)
    completedict= variant_dict|dict2

    return completedict



def make_dataframe(function,file):
    allelsdf=pandas.DataFrame(joind(joind,file))
    print(allelsdf)
    print(type(allelsdf))
    
    
       



#lines('filtered_GBR_biallelic.csv')
extract_fields(lines,'filtered_GBR_biallelic.csv')
#joind(extract_fields,'filtered_GBR_biallelic.csv')
#make_dataframe(joind, 'filtered_GBR_biallelic.csv')


## converting numbers to letters###





