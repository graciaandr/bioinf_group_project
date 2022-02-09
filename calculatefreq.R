gt<-read.csv('~/Queen mary uni/Group project coding/gt.csv',sep=",")
install.packages('genetics')

library('genetics')

vec_gt <-  c()


for (i in colnames(gt)){
  vec_gt<-c(vec_gt,as.vector(gt[i]))

  
}

vec_gt1.nosep <- c()
for (i in vec_gt)
  #print(i)
  { vec_gt1.nosep<-c(vec_gt1.nosep,i)
  #  #vec_gt3.nosep<-c(vec_gt2.nosep)
  
}


#print(vec_gt1.nosep)
#vec_gt2<-unlist(vec_gt)

#gt3<- c()

#for (i in vec_gt1.nosep)
 # { gt3<-c(gt3,genotype(i,sep = ""))}

genotype(vec_gt1.nosep)
gt_stats=c()

for (i in gt3)
  {gt_stats<- c(gt_stats, summary(i))}



gt_stats2<-as.matrix(gt_stats)


