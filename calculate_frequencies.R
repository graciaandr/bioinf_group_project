file = 'C:/Users/andri/Documents/Uni London/QMUL/SemesterB/BIO727P/Coding/22_GBR_finalGT.csv'
# df <- read.csv(file, header = T)
head(df)

library(dplyr)
library(magrittr)
library(tidyverse)
library(genetics)


df1 <- df %>% dplyr::select(-PK_ID, -CHROM, -POS, -REF, -ALT) 
df_test <- df1 # [1:1000,]
afs_ref = c()
afs_alt = c()
gt_00 = c()
gt_01_10 = c()
gt_11 = c()

n = (nrow(df_test))
start_time <- Sys.time()

for (i in (1:n)) {
  t = as.character(as.vector(df_test[i,]))
  print(i)
  g1 <- genotype(t, sep = "")
  sum_g1 <- summary(genotype(g1))
  af_g1 <- data.frame(t(sum_g1$allele.freq))
  final_af_g1_ref <- as.numeric(as.vector(af_g1[2,]))[1]
  final_af_g1_alt <- as.numeric(as.vector(af_g1[2,]))[2]
  
  gt_g1 <- data.frame(sum_g1$genotype.freq)

  if ((dim(gt_g1)[2]) >1) gt_g1 = t(gt_g1)

  final_gt_g1_00 <- as.numeric(as.vector(gt_g1[2,]))[1]
  final_gt_g1_01_10 <- as.numeric(as.vector(gt_g1[2,]))[2]
  final_gt_g1_11 <- as.numeric(as.vector(gt_g1[2,]))[3]
  
  
  afs_alt <- c(afs_alt, final_af_g1_alt)
  afs_ref <- c(afs_ref, final_af_g1_ref)
  gt_00 <- c(gt_00, final_gt_g1_00)
  gt_01_10 <- c(gt_01_10, final_gt_g1_01_10)
  gt_11 <- c(gt_11, final_gt_g1_11)
  
}

end_time <- Sys.time()
print(end_time - start_time)

dd = data.frame(AF_ref = afs_ref, AF_alt = afs_alt, GT00 = gt_00, GT0110 = gt_01_10, GT11 = gt_11)
dd[is.na(dd)] <- 0

View(dd) 
write.table(dd, file = "corrected_frequencies.csv", sep = ";", col.names = TRUE, row.names = FALSE)
