#install and load packages vegan, adegenet and pegas as 
##they have functions to calculate genetic and haplotype diversity,... 
###neutrality and Fst.

library(vegan)
library(pegas)
library(adegenet)
library(DBI)
library(dplyr)

#make a function that connects files from the sql database
portaldb_inpt <- dbconnect(RSQLite::SQLite(), '')

dbListTable(portaldb_inpt, '')
dbListFields(portaldb_inpt, '') 
surveys <- tbl(portaldb_inpt, '')

#Extract table from connected SQLite database
surveys.df <- collect(surveys) 

#calculate summary statistics
tajimas.D <- tajima.test(surveys.df)
haplotype_diversity <- hap.div(surveys.df)
genetic_diversity <- diversity(surveys.df, index  = 'shannon')

#make a summary statistics table/output file

#copy results back to SQLite database
copy_to(portaldb_inpt, #'name of table'#, temporary = FALSE) 
        
#check whether table is in SQLite database
dbListTable(portaldb_inpt, '') 









