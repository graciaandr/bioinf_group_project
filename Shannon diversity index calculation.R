install.packages("vegan")
library('vegan')
install.packages('tidyr')
library('tidyr')


############load data############# 


df<- read.csv('~/Queen mary uni/Groupproject_test/frequencies.csv')

#split single column to multiple based on ;
df2<- data.frame(AF_ref.AF_alt.GT00.GT01.GT10=df[1])
df3<-separate(data=df2, col=AF_ref.AF_alt.GT00.GT01.GT10, into=c('REF', 'ALT', 'GT00', 'GT01', 'GT10'), sep = ";")


# grabs only genotype freq
genotype_freq<- df3[,3:5]
str(genotype_freq)

#make test data
#mock_datafreq<- genotype_freq[1:20,]


## make input data numeric 
#mock_datafreq$GT00<-as.numeric(mock_datafreq$GT00)
#mock_datafreq$GT01<-as.numeric(mock_datafreq$GT01)
#mock_datafreq$GT10<-as.numeric(mock_datafreq$GT10)
#str(mock_datafreq)
#str(mock_datafreq)

#### requires dataframe input with column as numeric to work
#shannon_stats<-diversity(mock_datafreq,index = 'shannon')

#print(shannon_stats)


### make function to intake any dataframe of multiple snps from one population to find shannon diveirty index ####

calculate_shannonindex<-function(df_withGTfreq){
  for (i in colnames(df_withGTfreq)){ 
    df_withGTfreq[i]<-as.numeric(unlist(df_withGTfreq[i]))
    
  }  
  shannon_stats<-diversity(df_withGTfreq,index = 'shannon')
  
}
#ex<-calculate_shannonindex(mock_datafreq)
# check function output correct
#all.equal(ex,mock_datafreq)
#True
#ex_stats<-calculate_shannonindex(mock_datafreq)

# check function output correct

#all.equal(ex_stats,shannon_stats)

#True

full_set_stats<-calculate_shannonindex(genotype_freq)