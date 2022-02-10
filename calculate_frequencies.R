file = 'C:/Users/andri/Documents/Uni London/QMUL/SemesterB/BIO727P/Coding/22_GBR_renamed_genotype.csv'
# df <- read.csv(file, header = T)
head(df)

library(dplyr)
library(magrittr)
library(tidyverse)
library(genetics)


df1 <- df %>% dplyr::select(-PK_ID, -CHROM, -POS, -REF, -ALT) 
df_test <- df1 #[1:50,]
afs_ref = c()
afs_alt = c()
gt_00 = c()
gt_01 = c()
gt_10 = c()

n = (nrow(df_test))
# n =20
for (i in (1:n)) {
  t = as.character(as.vector(df_test[i,]))
  # print(t)
  # print(i)
  g1 <- genotype(t, sep = "")
  sum_g1 <- summary(genotype(g1))
  af_g1 <- data.frame(t(sum_g1$allele.freq))
  final_af_g1_ref <- as.numeric(as.vector(af_g1[2,]))[1]
  final_af_g1_alt <- as.numeric(as.vector(af_g1[2,]))[2]
  
  gt_g1 <- data.frame(sum_g1$genotype.freq)
  # print(dim(gt_g1))
  if ((dim(gt_g1)[2]) >1) gt_g1 = t(gt_g1)
  # print(gt_g1)
  
  final_gt_g1_00 <- as.numeric(as.vector(gt_g1[2,]))[1]
  final_gt_g1_01 <- as.numeric(as.vector(gt_g1[2,]))[2]
  final_gt_g1_10 <- as.numeric(as.vector(gt_g1[2,]))[3]
  
  
  afs_alt <- c(afs_alt, final_af_g1_alt)
  afs_ref <- c(afs_ref, final_af_g1_ref)
  gt_00 <- c(gt_00, final_gt_g1_00)
  gt_01 <- c(gt_01, final_gt_g1_01)
  gt_10 <- c(gt_10, final_gt_g1_10)
  
}
# 
# print(afs_ref)
# print(afs_alt)
# 
# print(gt_00)
# print(gt_01)
# gt_10

d = data.frame(AF_ref = afs_ref, AF_alt = afs_alt, GT00 = gt_00, GT01 = gt_01, GT10 = gt_10)
d[is.na(d)] <- 0

d
