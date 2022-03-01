library(dplyr)
library(magrittr)
library(tidyverse)
library(DBI)

file = 'gene_names_filtered_for_PK.csv'
df <- read.csv(file, header = T)
nrow(df)
head(df)
summary(df)

# install Bioconductor and the database 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")

# load the annotation database
library(org.Hs.eg.db)

# set up query genes
queryGeneNames <- unique(df$GENE)

# convert to string aand gene names separated by logical or (for query later)
queryGeneNames <- (paste(queryGeneNames, collapse="|^")) 

# use sql to get alias table and gene_info table (contains the symbols)
# first open the(paste(queryGeneNames, collapse="|^")) database connection
dbCon <- org.Hs.eg_dbconn()
# SQL query
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
aliasSymbol <- DBI::dbGetQuery(dbCon, sqlQuery)

# subset to get your results
result <-aliasSymbol[grepl(queryGeneNames,aliasSymbol$symbol),]
result <- result %>% dplyr::select(-"_id")
View(result)

# rename column names (old: "alias_symbol", "_id", "gene_name", "symbol")
colnames(result) <- c("alias", "index", "complete_gene_name", "gene_symbol")
# write.csv(result, file = "gene_name_alias.csv", row.names = F)

# script modified from https://www.biostars.org/p/14971/
