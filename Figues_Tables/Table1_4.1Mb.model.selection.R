setwd("~/Dropbox/paper1_genomic_summary/data/1Mb")

total_1Mb=read.table("1Mb.summary.txt",header=T)
total_1Mb=total_1Mb[which(total_1Mb$Thetas_tremula_numsites>100000),]
total_1Mb$tremula_ldhat_new=total_1Mb$Ldhat_tremula/1000
total_1Mb$tremuloides_ldhat_new=total_1Mb$Ldhat_tremuloides/1000
total_1Mb$trichocarpa_ldhat_new=total_1Mb$Ldhat_trichocarpa/1000
#total_1Mb$tremula_alpha=(total_1Mb$tremula_trichocarpa_zero_fold_dxy/total_1Mb$Thetas_tremula_zero_fold_tP)/(total_1Mb$tremula_trichocarpa_four_fold_dxy/total_1Mb$Thetas_tremula_four_fold_tP)

total_1Mb[which(total_1Mb$Ldhat_tremula_num<100),]$tremula_ldhat_new="NA"
total_1Mb[which(total_1Mb$Ldhat_tremuloides_num<100),]$tremuloides_ldhat_new="NA"
total_1Mb[which(total_1Mb$Ldhat_trichocarpa_num<100),]$trichocarpa_ldhat_new="NA"
total_1Mb$tremula_ldhat_new=as.numeric(as.character(total_1Mb$tremula_ldhat_new))
total_1Mb$tremuloides_ldhat_new=as.numeric(as.character(total_1Mb$tremuloides_ldhat_new))
total_1Mb$trichocarpa_ldhat_new=as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new))

attach(total_1Mb)
detach(total_1Mb)
model1=lm(Thetas_tremula_four_fold_tP~GC+Gene_num+tremula_ldhat_new+tremula_trichocarpa_four_fold_fixed+tremula_trichocarpa_zero_fold_fixed)
summary(model1)
model2=aov(Thetas_tremula_four_fold_tP~GC+Gene_num+tremula_ldhat_new+tremula_trichocarpa_four_fold_dxy+tremula_trichocarpa_zero_fold_dxy)
summary(model2)


###partial correlations
####neutral diveristy~Gene number+divergence (mutation)+GC content+Number of neutral sites in the window+recombination rate
#install.packages("ppcor")
library(ppcor)
tremula_1Mb_four_fold=total_1Mb[,c("Thetas_tremula_four_fold_tP","Thetas_tremula_four_fold_numsites","tremula_ldhat_new","GC","Gene_num","tremula_trichocarpa_four_fold_fixed")]
tremula_1Mb_four_fold_new=tremula_1Mb_four_fold[complete.cases(tremula_1Mb_four_fold),]

cor.test(tremula_1Mb_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_four_fold_new$tremula_ldhat_new)
pcor.test(tremula_1Mb_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_four_fold_new$tremula_ldhat_new,tremula_1Mb_four_fold_new[,-c(1,3)])
cor.test(tremula_1Mb_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_four_fold_new$tremula_ldhat_new,method="spearman")
pcor.test(tremula_1Mb_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_four_fold_new$tremula_ldhat_new,tremula_1Mb_four_fold_new[,-c(1,3)],method="spearman")
cor.test(tremula_1Mb_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_four_fold_new$Gene_num,method="spearman")
pcor.test(tremula_1Mb_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_four_fold_new$Gene_num,tremula_1Mb_four_fold_new[,-c(1,5)],method="spearman")

##small and large gene numbers
tremula_1Mb_four_fold_new_small=tremula_1Mb_four_fold_new[which(tremula_1Mb_four_fold_new$Gene_num<85),]
cor.test(tremula_1Mb_four_fold_new_small$Thetas_tremula_four_fold_tP,tremula_1Mb_four_fold_new_small$Gene_num,method="spearman")
pcor.test(tremula_1Mb_four_fold_new_small$Thetas_tremula_four_fold_tP,tremula_1Mb_four_fold_new_small$Gene_num,tremula_1Mb_four_fold_new_small[,-c(1,5)],method="spearman")
pcor.test(tremula_1Mb_four_fold_new_small$Thetas_tremula_four_fold_tP,tremula_1Mb_four_fold_new_small$Gene_num,tremula_1Mb_four_fold_new_small$tremula_ldhat_new,method="spearman")



tremula_1Mb_four_fold_new_large=tremula_1Mb_four_fold_new[which(tremula_1Mb_four_fold_new$Gene_num>=85),]
cor.test(tremula_1Mb_four_fold_new_large$Thetas_tremula_four_fold_tP,tremula_1Mb_four_fold_new_large$Gene_num,method="spearman")
pcor.test(tremula_1Mb_four_fold_new_large$Thetas_tremula_four_fold_tP,tremula_1Mb_four_fold_new_large$Gene_num,tremula_1Mb_four_fold_new_large[,-c(1,5)],method="spearman")
pcor.test(tremula_1Mb_four_fold_new_large$Thetas_tremula_four_fold_tP,tremula_1Mb_four_fold_new_large$Gene_num,tremula_1Mb_four_fold_new_large$tremula_ldhat_new,method="spearman")



tremula_1Mb_intergenic=total_1Mb[,c("Thetas_tremula_intergenic_tP","Thetas_tremula_intergenic_numsites","tremula_ldhat_new","GC","Gene_num","tremula_trichocarpa_intergenic_fixed")]
tremula_1Mb_intergenic_new=tremula_1Mb_intergenic[complete.cases(tremula_1Mb_intergenic),]
cor.test(tremula_1Mb_intergenic_new$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_new$tremula_ldhat_new)
pcor.test(tremula_1Mb_intergenic_new$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_new$tremula_ldhat_new,tremula_1Mb_intergenic_new[,-c(1,3)],method="pearson")
cor.test(tremula_1Mb_intergenic_new$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_new$tremula_ldhat_new,method="spearman")
pcor.test(tremula_1Mb_intergenic_new$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_new$tremula_ldhat_new,tremula_1Mb_intergenic_new[,-c(1,3)],method="spearman")
cor.test(tremula_1Mb_intergenic_new$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_new$Gene_num,method="spearman")
pcor.test(tremula_1Mb_intergenic_new$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_new$Gene_num,tremula_1Mb_intergenic_new[,-c(1,5)],method="spearman")

##small and large gene numbers
tremula_1Mb_intergenic_new_small=tremula_1Mb_intergenic_new[which(tremula_1Mb_intergenic_new$Gene_num<85),]
cor.test(tremula_1Mb_intergenic_new_small$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_new_small$Gene_num,method="spearman")
pcor.test(tremula_1Mb_intergenic_new_small$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_new_small$Gene_num,tremula_1Mb_intergenic_new_small[,-c(1,5)],method="spearman")
pcor.test(tremula_1Mb_intergenic_new_small$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_new_small$Gene_num,tremula_1Mb_intergenic_new_small$tremula_ldhat_new,method="spearman")


tremula_1Mb_intergenic_new_large=tremula_1Mb_intergenic_new[which(tremula_1Mb_intergenic_new$Gene_num>=85),]
cor.test(tremula_1Mb_intergenic_new_large$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_new_large$Gene_num,method="spearman")
pcor.test(tremula_1Mb_intergenic_new_large$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_new_large$Gene_num,tremula_1Mb_intergenic_new_large[,-c(1,5)],method="spearman")
pcor.test(tremula_1Mb_intergenic_new_large$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_new_large$Gene_num,tremula_1Mb_intergenic_new_large$tremula_ldhat_new,method="spearman")


tremula_1Mb_four_zero_fold=total_1Mb[,c("Thetas_tremula_four_fold_tP","Thetas_tremula_four_fold_numsites","tremula_ldhat_new","GC","Gene_num","tremula_trichocarpa_four_fold_fixed","tremula_trichocarpa_zero_fold_fixed")]
tremula_1Mb_four_zero_fold_new=tremula_1Mb_four_zero_fold[complete.cases(tremula_1Mb_four_zero_fold),]
cor.test(tremula_1Mb_four_zero_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_four_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,method="spearman")
pcor.test(tremula_1Mb_four_zero_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_four_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,tremula_1Mb_four_zero_fold_new[,-c(1,7)],method="spearman")

tremula_1Mb_intergenic_zero_fold=total_1Mb[,c("Thetas_tremula_intergenic_tP","Thetas_tremula_intergenic_numsites","tremula_ldhat_new","GC","Gene_num","tremula_trichocarpa_intergenic_fixed","tremula_trichocarpa_zero_fold_fixed")]
tremula_1Mb_intergenic_zero_fold_new=tremula_1Mb_intergenic_zero_fold[complete.cases(tremula_1Mb_intergenic_zero_fold),]
cor.test(tremula_1Mb_intergenic_zero_fold_new$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,method="spearman")
pcor.test(tremula_1Mb_intergenic_zero_fold_new$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,tremula_1Mb_intergenic_zero_fold_new[,-c(1,7)],method="spearman")



#####tremuloides
tremuloides_1Mb_four_fold=total_1Mb[,c("Thetas_tremuloides_four_fold_tP","Thetas_tremuloides_four_fold_numsites","tremuloides_ldhat_new","GC","Gene_num","tremuloides_trichocarpa_four_fold_fixed")]
tremuloides_1Mb_four_fold_new=tremuloides_1Mb_four_fold[complete.cases(tremuloides_1Mb_four_fold),]
cor.test(tremuloides_1Mb_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_four_fold_new$tremuloides_ldhat_new)
pcor.test(tremuloides_1Mb_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_four_fold_new$tremuloides_ldhat_new,tremuloides_1Mb_four_fold_new[,-c(1,3)])
cor.test(tremuloides_1Mb_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_four_fold_new$tremuloides_ldhat_new,method="spearman")
pcor.test(tremuloides_1Mb_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_four_fold_new$tremuloides_ldhat_new,tremuloides_1Mb_four_fold_new[,-c(1,3)],method="spearman")
cor.test(tremuloides_1Mb_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_four_fold_new$Gene_num,method="spearman")
pcor.test(tremuloides_1Mb_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_four_fold_new$Gene_num,tremuloides_1Mb_four_fold_new[,-c(1,5)],method="spearman")

##small and large gene numbers
tremuloides_1Mb_four_fold_new_small=tremuloides_1Mb_four_fold_new[which(tremuloides_1Mb_four_fold_new$Gene_num<85),]
cor.test(tremuloides_1Mb_four_fold_new_small$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_four_fold_new_small$Gene_num,method="spearman")
pcor.test(tremuloides_1Mb_four_fold_new_small$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_four_fold_new_small$Gene_num,tremuloides_1Mb_four_fold_new_small[,-c(1,5)],method="spearman")
pcor.test(tremuloides_1Mb_four_fold_new_small$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_four_fold_new_small$Gene_num,tremuloides_1Mb_four_fold_new_small$tremuloides_ldhat_new,method="spearman")



tremuloides_1Mb_four_fold_new_large=tremuloides_1Mb_four_fold_new[which(tremuloides_1Mb_four_fold_new$Gene_num>=85),]
cor.test(tremuloides_1Mb_four_fold_new_large$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_four_fold_new_large$Gene_num,method="spearman")
pcor.test(tremuloides_1Mb_four_fold_new_large$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_four_fold_new_large$Gene_num,tremuloides_1Mb_four_fold_new_large[,-c(1,5)],method="spearman")
pcor.test(tremuloides_1Mb_four_fold_new_large$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_four_fold_new_large$Gene_num,tremuloides_1Mb_four_fold_new_large$tremuloides_ldhat_new,method="spearman")


##intergenic
tremuloides_1Mb_intergenic=total_1Mb[,c("Thetas_tremuloides_intergenic_tP","Thetas_tremuloides_intergenic_numsites","tremuloides_ldhat_new","GC","Gene_num","tremuloides_trichocarpa_intergenic_fixed")]
tremuloides_1Mb_intergenic_new=tremuloides_1Mb_intergenic[complete.cases(tremuloides_1Mb_intergenic),]
cor.test(tremuloides_1Mb_intergenic_new$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_new$tremuloides_ldhat_new)
pcor.test(tremuloides_1Mb_intergenic_new$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_new$tremuloides_ldhat_new,tremuloides_1Mb_intergenic_new[,-c(1,3)])
cor.test(tremuloides_1Mb_intergenic_new$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_new$tremuloides_ldhat_new,method="spearman")
pcor.test(tremuloides_1Mb_intergenic_new$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_new$tremuloides_ldhat_new,tremuloides_1Mb_intergenic_new[,-c(1,3)],method="spearman")
cor.test(tremuloides_1Mb_intergenic_new$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_new$Gene_num,method="spearman")
pcor.test(tremuloides_1Mb_intergenic_new$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_new$Gene_num,tremuloides_1Mb_intergenic_new[,-c(1,5)],method="spearman")

##small and large gene numbers
tremuloides_1Mb_intergenic_new_small=tremuloides_1Mb_intergenic_new[which(tremuloides_1Mb_intergenic_new$Gene_num<85),]
cor.test(tremuloides_1Mb_intergenic_new_small$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_new_small$Gene_num,method="spearman")
pcor.test(tremuloides_1Mb_intergenic_new_small$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_new_small$Gene_num,tremuloides_1Mb_intergenic_new_small[,-c(1,5)],method="spearman")
pcor.test(tremuloides_1Mb_intergenic_new_small$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_new_small$Gene_num,tremuloides_1Mb_intergenic_new_small$tremuloides_ldhat_new,method="spearman")


tremuloides_1Mb_intergenic_new_large=tremuloides_1Mb_intergenic_new[which(tremuloides_1Mb_intergenic_new$Gene_num>=85),]
cor.test(tremuloides_1Mb_intergenic_new_large$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_new_large$Gene_num,method="spearman")
pcor.test(tremuloides_1Mb_intergenic_new_large$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_new_large$Gene_num,tremuloides_1Mb_intergenic_new_large[,-c(1,5)],method="spearman")
pcor.test(tremuloides_1Mb_intergenic_new_small$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_new_small$Gene_num,tremuloides_1Mb_intergenic_new_small$tremuloides_ldhat_new,method="spearman")


#attach(tremuloides_1Mb_intergenic_new)
#detach(tremuloides_1Mb_intergenic_new)
#model2=lm(Thetas_tremuloides_intergenic_tP~GC+Gene_num+tremuloides_ldhat_new+Thetas_tremuloides_intergenic_numsites+tremuloides_trichocarpa_intergenic_fixed)

#####trichocarpa
trichocarpa_1Mb_four_fold=total_1Mb[,c("Thetas_trichocarpa_four_fold_tP","Thetas_trichocarpa_four_fold_numsites","trichocarpa_ldhat_new","GC","Gene_num","tremula_trichocarpa_four_fold_fixed")]
trichocarpa_1Mb_four_fold_new=trichocarpa_1Mb_four_fold[complete.cases(trichocarpa_1Mb_four_fold),]
cor.test(trichocarpa_1Mb_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_new$trichocarpa_ldhat_new)
pcor.test(trichocarpa_1Mb_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_new$trichocarpa_ldhat_new,trichocarpa_1Mb_four_fold_new[,-c(1,3)])
cor.test(trichocarpa_1Mb_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_new$trichocarpa_ldhat_new,method="spearman")
pcor.test(trichocarpa_1Mb_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_new$trichocarpa_ldhat_new,trichocarpa_1Mb_four_fold_new[,-c(1,3)],method="spearman")
cor.test(trichocarpa_1Mb_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_new$Gene_num,method="spearman")
pcor.test(trichocarpa_1Mb_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_new$Gene_num,trichocarpa_1Mb_four_fold_new[,-c(1,5)],method="spearman")

##small and large gene numbers
trichocarpa_1Mb_four_fold_new_small=trichocarpa_1Mb_four_fold_new[which(trichocarpa_1Mb_four_fold_new$Gene_num<85),]
cor.test(trichocarpa_1Mb_four_fold_new_small$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_new_small$Gene_num,method="spearman")
pcor.test(trichocarpa_1Mb_four_fold_new_small$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_new_small$Gene_num,trichocarpa_1Mb_four_fold_new_small[,-c(1,5)],method="spearman")
pcor.test(trichocarpa_1Mb_four_fold_new_small$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_new_small$Gene_num,trichocarpa_1Mb_four_fold_new_small$trichocarpa_ldhat_new,method="spearman")



trichocarpa_1Mb_four_fold_new_large=trichocarpa_1Mb_four_fold_new[which(trichocarpa_1Mb_four_fold_new$Gene_num>=85),]
cor.test(trichocarpa_1Mb_four_fold_new_large$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_new_large$Gene_num,method="spearman")
pcor.test(trichocarpa_1Mb_four_fold_new_large$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_new_large$Gene_num,trichocarpa_1Mb_four_fold_new_large[,-c(1,5)],method="spearman")
pcor.test(trichocarpa_1Mb_four_fold_new_large$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_new_large$Gene_num,trichocarpa_1Mb_four_fold_new_large$trichocarpa_ldhat_new,method="spearman")

###intergenic 
trichocarpa_1Mb_intergenic=total_1Mb[,c("Thetas_trichocarpa_intergenic_tP","Thetas_trichocarpa_intergenic_numsites","trichocarpa_ldhat_new","GC","Gene_num","tremula_trichocarpa_intergenic_fixed")]
trichocarpa_1Mb_intergenic_new=trichocarpa_1Mb_intergenic[complete.cases(trichocarpa_1Mb_intergenic),]
cor.test(trichocarpa_1Mb_intergenic_new$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_new$trichocarpa_ldhat_new)
pcor.test(trichocarpa_1Mb_intergenic_new$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_new$trichocarpa_ldhat_new,trichocarpa_1Mb_intergenic_new[,-c(1,3)])
cor.test(trichocarpa_1Mb_intergenic_new$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_new$trichocarpa_ldhat_new,method="spearman")
pcor.test(trichocarpa_1Mb_intergenic_new$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_new$trichocarpa_ldhat_new,trichocarpa_1Mb_intergenic_new[,-c(1,3)],method="spearman")
cor.test(trichocarpa_1Mb_intergenic_new$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_new$Gene_num,method="spearman")
pcor.test(trichocarpa_1Mb_intergenic_new$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_new$Gene_num,trichocarpa_1Mb_intergenic_new[,-c(1,5)],method="spearman")

##small and large gene numbers
trichocarpa_1Mb_intergenic_new_small=trichocarpa_1Mb_intergenic_new[which(trichocarpa_1Mb_intergenic_new$Gene_num<85),]
cor.test(trichocarpa_1Mb_intergenic_new_small$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_new_small$Gene_num,method="spearman")
pcor.test(trichocarpa_1Mb_intergenic_new_small$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_new_small$Gene_num,trichocarpa_1Mb_intergenic_new_small[,-c(1,5)],method="spearman")
pcor.test(trichocarpa_1Mb_intergenic_new_small$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_new_small$Gene_num,trichocarpa_1Mb_intergenic_new_small$trichocarpa_ldhat_new,method="spearman")


trichocarpa_1Mb_intergenic_new_large=trichocarpa_1Mb_intergenic_new[which(trichocarpa_1Mb_intergenic_new$Gene_num>=85),]
cor.test(trichocarpa_1Mb_intergenic_new_large$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_new_large$Gene_num,method="spearman")
pcor.test(trichocarpa_1Mb_intergenic_new_large$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_new_large$Gene_num,trichocarpa_1Mb_intergenic_new_large[,-c(1,5)],method="spearman")
pcor.test(trichocarpa_1Mb_intergenic_new_large$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_new_large$Gene_num,trichocarpa_1Mb_intergenic_new_large$trichocarpa_ldhat_new,method="spearman")


######tajD

tremula_1Mb_four_fold_tajD=total_1Mb[,c("Thetas_tremula_four_fold_tajD","Thetas_tremula_four_fold_numsites","tremula_ldhat_new","GC","Gene_num","tremula_trichocarpa_four_fold_fixed")]
tremula_1Mb_four_fold_tajD_new=tremula_1Mb_four_fold_tajD[complete.cases(tremula_1Mb_four_fold_tajD),]
cor.test(tremula_1Mb_four_fold_tajD_new$Thetas_tremula_four_fold_tajD,tremula_1Mb_four_fold_tajD_new$tremula_ldhat_new)
pcor.test(tremula_1Mb_four_fold_tajD_new$Thetas_tremula_four_fold_tajD,tremula_1Mb_four_fold_tajD_new$tremula_ldhat_new,tremula_1Mb_four_fold_tajD_new[,-c(1,3)])
cor.test(tremula_1Mb_four_fold_tajD_new$Thetas_tremula_four_fold_tajD,tremula_1Mb_four_fold_tajD_new$tremula_ldhat_new,method="spearman")
pcor.test(tremula_1Mb_four_fold_tajD_new$Thetas_tremula_four_fold_tajD,tremula_1Mb_four_fold_tajD_new$tremula_ldhat_new,tremula_1Mb_four_fold_tajD_new[,-c(1,3)],method="spearman")
cor.test(tremula_1Mb_four_fold_tajD_new$Thetas_tremula_four_fold_tajD,tremula_1Mb_four_fold_tajD_new$Gene_num,method="spearman")
pcor.test(tremula_1Mb_four_fold_tajD_new$Thetas_tremula_four_fold_tajD,tremula_1Mb_four_fold_tajD_new$Gene_num,tremula_1Mb_four_fold_tajD_new[,-c(1,5)],method="spearman")

tremula_1Mb_intergenic_tajD=total_1Mb[,c("Thetas_tremula_intergenic_tajD","Thetas_tremula_intergenic_numsites","tremula_ldhat_new","GC","Gene_num","tremula_trichocarpa_intergenic_fixed")]
tremula_1Mb_intergenic_tajD_new=tremula_1Mb_intergenic_tajD[complete.cases(tremula_1Mb_intergenic),]
cor.test(tremula_1Mb_intergenic_tajD_new$Thetas_tremula_intergenic_tajD,tremula_1Mb_intergenic_tajD_new$tremula_ldhat_new)
pcor.test(tremula_1Mb_intergenic_tajD_new$Thetas_tremula_intergenic_tajD,tremula_1Mb_intergenic_tajD_new$tremula_ldhat_new,tremula_1Mb_intergenic_tajD_new[,-c(1,3)],method="pearson")
cor.test(tremula_1Mb_intergenic_tajD_new$Thetas_tremula_intergenic_tajD,tremula_1Mb_intergenic_tajD_new$tremula_ldhat_new,method="spearman")
pcor.test(tremula_1Mb_intergenic_tajD_new$Thetas_tremula_intergenic_tajD,tremula_1Mb_intergenic_tajD_new$tremula_ldhat_new,tremula_1Mb_intergenic_tajD_new[,-c(1,3)],method="spearman")
cor.test(tremula_1Mb_intergenic_tajD_new$Thetas_tremula_intergenic_tajD,tremula_1Mb_intergenic_tajD_new$Gene_num,method="spearman")
pcor.test(tremula_1Mb_intergenic_tajD_new$Thetas_tremula_intergenic_tajD,tremula_1Mb_intergenic_tajD_new$Gene_num,tremula_1Mb_intergenic_tajD_new[,-c(1,5)],method="spearman")

tremula_1Mb_four_zero_fold=total_1Mb[,c("Thetas_tremula_four_fold_tP","Thetas_tremula_four_fold_numsites","tremula_ldhat_new","GC","Gene_num","tremula_trichocarpa_four_fold_fixed","tremula_trichocarpa_zero_fold_fixed")]
tremula_1Mb_four_zero_fold_new=tremula_1Mb_four_zero_fold[complete.cases(tremula_1Mb_four_zero_fold),]
cor.test(tremula_1Mb_four_zero_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_four_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,method="spearman")
pcor.test(tremula_1Mb_four_zero_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_four_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,tremula_1Mb_four_zero_fold_new[,-c(1,7)],method="spearman")

tremula_1Mb_intergenic_zero_fold=total_1Mb[,c("Thetas_tremula_intergenic_tP","Thetas_tremula_intergenic_numsites","tremula_ldhat_new","GC","Gene_num","tremula_trichocarpa_intergenic_fixed","tremula_trichocarpa_zero_fold_fixed")]
tremula_1Mb_intergenic_zero_fold_new=tremula_1Mb_intergenic_zero_fold[complete.cases(tremula_1Mb_intergenic_zero_fold),]
cor.test(tremula_1Mb_intergenic_zero_fold_new$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,method="spearman")
pcor.test(tremula_1Mb_intergenic_zero_fold_new$Thetas_tremula_intergenic_tP,tremula_1Mb_intergenic_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,tremula_1Mb_intergenic_zero_fold_new[,-c(1,7)],method="spearman")



#####tremuloides
tremuloides_1Mb_four_fold_tajD=total_1Mb[,c("Thetas_tremuloides_four_fold_tajD","Thetas_tremuloides_four_fold_numsites","tremuloides_ldhat_new","GC","Gene_num","tremuloides_trichocarpa_four_fold_fixed")]
tremuloides_1Mb_four_fold_tajD_new=tremuloides_1Mb_four_fold_tajD[complete.cases(tremuloides_1Mb_four_fold),]
cor.test(tremuloides_1Mb_four_fold_tajD_new$Thetas_tremuloides_four_fold_tajD,tremuloides_1Mb_four_fold_tajD_new$tremuloides_ldhat_new)
pcor.test(tremuloides_1Mb_four_fold_tajD_new$Thetas_tremuloides_four_fold_tajD,tremuloides_1Mb_four_fold_tajD_new$tremuloides_ldhat_new,tremuloides_1Mb_four_fold_tajD_new[,-c(1,3)])
cor.test(tremuloides_1Mb_four_fold_tajD_new$Thetas_tremuloides_four_fold_tajD,tremuloides_1Mb_four_fold_tajD_new$tremuloides_ldhat_new,method="spearman")
pcor.test(tremuloides_1Mb_four_fold_tajD_new$Thetas_tremuloides_four_fold_tajD,tremuloides_1Mb_four_fold_tajD_new$tremuloides_ldhat_new,tremuloides_1Mb_four_fold_tajD_new[,-c(1,3)],method="spearman")
cor.test(tremuloides_1Mb_four_fold_tajD_new$Thetas_tremuloides_four_fold_tajD,tremuloides_1Mb_four_fold_tajD_new$Gene_num,method="spearman")
pcor.test(tremuloides_1Mb_four_fold_tajD_new$Thetas_tremuloides_four_fold_tajD,tremuloides_1Mb_four_fold_tajD_new$Gene_num,tremuloides_1Mb_four_fold_tajD_new[,-c(1,5)],method="spearman")

tremuloides_1Mb_intergenic_tajD=total_1Mb[,c("Thetas_tremuloides_intergenic_tajD","Thetas_tremuloides_intergenic_numsites","tremuloides_ldhat_new","GC","Gene_num","tremuloides_trichocarpa_intergenic_fixed")]
tremuloides_1Mb_intergenic_tajD_new=tremuloides_1Mb_intergenic_tajD[complete.cases(tremuloides_1Mb_intergenic),]
cor.test(tremuloides_1Mb_intergenic_tajD_new$Thetas_tremuloides_intergenic_tajD,tremuloides_1Mb_intergenic_tajD_new$tremuloides_ldhat_new)
pcor.test(tremuloides_1Mb_intergenic_tajD_new$Thetas_tremuloides_intergenic_tajD,tremuloides_1Mb_intergenic_tajD_new$tremuloides_ldhat_new,tremuloides_1Mb_intergenic_tajD_new[,-c(1,3)])
cor.test(tremuloides_1Mb_intergenic_tajD_new$Thetas_tremuloides_intergenic_tajD,tremuloides_1Mb_intergenic_tajD_new$tremuloides_ldhat_new,method="spearman")
pcor.test(tremuloides_1Mb_intergenic_tajD_new$Thetas_tremuloides_intergenic_tajD,tremuloides_1Mb_intergenic_tajD_new$tremuloides_ldhat_new,tremuloides_1Mb_intergenic_tajD_new[,-c(1,3)],method="spearman")
cor.test(tremuloides_1Mb_intergenic_tajD_new$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_tajD_new$Gene_num,method="spearman")
pcor.test(tremuloides_1Mb_intergenic_tajD_new$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_intergenic_tajD_new$Gene_num,tremuloides_1Mb_intergenic_tajD_new[,-c(1,5)],method="spearman")

#attach(tremuloides_1Mb_intergenic_tajD_new)
#detach(tremuloides_1Mb_intergenic_tajD_new)
#model2=lm(Thetas_tremuloides_intergenic_tP~GC+Gene_num+tremuloides_ldhat_new+Thetas_tremuloides_intergenic_numsites+tremuloides_trichocarpa_intergenic_fixed)

#####trichocarpa
trichocarpa_1Mb_four_fold_tajD=total_1Mb[,c("Thetas_trichocarpa_four_fold_tajD","Thetas_trichocarpa_four_fold_numsites","trichocarpa_ldhat_new","GC","Gene_num","tremula_trichocarpa_four_fold_fixed")]
trichocarpa_1Mb_four_fold_tajD_new=trichocarpa_1Mb_four_fold_tajD[complete.cases(trichocarpa_1Mb_four_fold),]
cor.test(trichocarpa_1Mb_four_fold_tajD_new$Thetas_trichocarpa_four_fold_tajD,trichocarpa_1Mb_four_fold_tajD_new$trichocarpa_ldhat_new)
pcor.test(trichocarpa_1Mb_four_fold_tajD_new$Thetas_trichocarpa_four_fold_tajD,trichocarpa_1Mb_four_fold_tajD_new$trichocarpa_ldhat_new,trichocarpa_1Mb_four_fold_tajD_new[,-c(1,3)])
cor.test(trichocarpa_1Mb_four_fold_tajD_new$Thetas_trichocarpa_four_fold_tajD,trichocarpa_1Mb_four_fold_tajD_new$trichocarpa_ldhat_new,method="spearman")
pcor.test(trichocarpa_1Mb_four_fold_tajD_new$Thetas_trichocarpa_four_fold_tajD,trichocarpa_1Mb_four_fold_tajD_new$trichocarpa_ldhat_new,trichocarpa_1Mb_four_fold_tajD_new[,-c(1,3)],method="spearman")
cor.test(trichocarpa_1Mb_four_fold_tajD_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_tajD_new$Gene_num,method="spearman")
pcor.test(trichocarpa_1Mb_four_fold_tajD_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_four_fold_tajD_new$Gene_num,trichocarpa_1Mb_four_fold_tajD_new[,-c(1,5)],method="spearman")

trichocarpa_1Mb_intergenic_tajD=total_1Mb[,c("Thetas_trichocarpa_intergenic_tajD","Thetas_trichocarpa_intergenic_numsites","trichocarpa_ldhat_new","GC","Gene_num","tremula_trichocarpa_intergenic_fixed")]
trichocarpa_1Mb_intergenic_tajD_new=trichocarpa_1Mb_intergenic_tajD[complete.cases(trichocarpa_1Mb_intergenic),]
cor.test(trichocarpa_1Mb_intergenic_tajD_new$Thetas_trichocarpa_intergenic_tajD,trichocarpa_1Mb_intergenic_tajD_new$trichocarpa_ldhat_new)
pcor.test(trichocarpa_1Mb_intergenic_tajD_new$Thetas_trichocarpa_intergenic_tajD,trichocarpa_1Mb_intergenic_tajD_new$trichocarpa_ldhat_new,trichocarpa_1Mb_intergenic_tajD_new[,-c(1,3)])
cor.test(trichocarpa_1Mb_intergenic_tajD_new$Thetas_trichocarpa_intergenic_tajD,trichocarpa_1Mb_intergenic_tajD_new$trichocarpa_ldhat_new,method="spearman")
pcor.test(trichocarpa_1Mb_intergenic_tajD_new$Thetas_trichocarpa_intergenic_tajD,trichocarpa_1Mb_intergenic_tajD_new$trichocarpa_ldhat_new,trichocarpa_1Mb_intergenic_tajD_new[,-c(1,3)],method="spearman")
cor.test(trichocarpa_1Mb_intergenic_tajD_new$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_tajD_new$Gene_num,method="spearman")
pcor.test(trichocarpa_1Mb_intergenic_tajD_new$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_intergenic_tajD_new$Gene_num,trichocarpa_1Mb_intergenic_tajD_new[,-c(1,5)],method="spearman")


####Gene number and recombination rate

#tremula
tremula_1Mb_gene=total_1Mb[,c("Thetas_tremula_numsites","tremula_ldhat_new","GC","Gene_num")]
tremula_1Mb_gene_new=tremula_1Mb_gene[complete.cases(tremula_1Mb_gene),]
tremula_1Mb_gene_new_small=tremula_1Mb_gene_new[which(tremula_1Mb_gene_new$Gene_num<85),]
tremula_1Mb_gene_new_large=tremula_1Mb_gene_new[which(tremula_1Mb_gene_new$Gene_num>=85),]

cor.test(tremula_1Mb_gene_new$Gene_num,tremula_1Mb_gene_new$tremula_ldhat_new)
cor.test(tremula_1Mb_gene_new$Gene_num,tremula_1Mb_gene_new$tremula_ldhat_new,method="spearman")
pcor.test(tremula_1Mb_gene_new$Gene_num,tremula_1Mb_gene_new$tremula_ldhat_new,tremula_1Mb_gene_new[,-c(2,4)],method="spearman")

cor.test(tremula_1Mb_gene_new_small$Gene_num,tremula_1Mb_gene_new_small$tremula_ldhat_new,method="spearman")
pcor.test(tremula_1Mb_gene_new_small$Gene_num,tremula_1Mb_gene_new_small$tremula_ldhat_new,tremula_1Mb_gene_new_small[,-c(2,4)],method="spearman")

cor.test(tremula_1Mb_gene_new_large$Gene_num,tremula_1Mb_gene_new_large$tremula_ldhat_new,method="spearman")
pcor.test(tremula_1Mb_gene_new_large$Gene_num,tremula_1Mb_gene_new_large$tremula_ldhat_new,tremula_1Mb_gene_new_large[,-c(2,4)],method="spearman")


cor.test(tremula_1Mb_gene_new$Gene_num,tremula_1Mb_gene_new$tremula_ldhat_new,method="spearman")
cor.test(tremula_1Mb_gene_new$GC,tremula_1Mb_gene_new$tremula_ldhat_new,method="spearman")
cor.test(tremula_1Mb_gene_new$Thetas_tremula_numsites,tremula_1Mb_gene_new$tremula_ldhat_new,method="spearman")

cor.test(tremula_1Mb_gene_new$Gene_num,tremula_1Mb_gene_new$GC,method="spearman")
cor.test(tremula_1Mb_gene_new$Gene_num,tremula_1Mb_gene_new$Thetas_tremula_numsites,method="spearman")

cor.test(tremula_1Mb_gene_new$GC,tremula_1Mb_gene_new$Thetas_tremula_numsites,method="spearman")

##tremuloides
tremuloides_1Mb_gene=total_1Mb[,c("Thetas_tremuloides_numsites","tremuloides_ldhat_new","GC","Gene_num")]
tremuloides_1Mb_gene_new=tremuloides_1Mb_gene[complete.cases(tremuloides_1Mb_gene),]
cor.test(tremuloides_1Mb_gene_new$Gene_num,tremuloides_1Mb_gene_new$tremuloides_ldhat_new)
cor.test(tremuloides_1Mb_gene_new$Gene_num,tremuloides_1Mb_gene_new$tremuloides_ldhat_new,method="spearman")
pcor.test(tremuloides_1Mb_gene_new$Gene_num,tremuloides_1Mb_gene_new$tremuloides_ldhat_new,tremuloides_1Mb_gene_new[,-c(2,4)],method="spearman")
pcor.test(tremuloides_1Mb_gene_new$Gene_num,tremuloides_1Mb_gene_new$tremuloides_ldhat_new,tremuloides_1Mb_gene_new[,-c(2,4)])

###small and large gene number
tremuloides_1Mb_gene_new_small=tremuloides_1Mb_gene_new[which(tremuloides_1Mb_gene_new$Gene_num<85),]
tremuloides_1Mb_gene_new_large=tremuloides_1Mb_gene_new[which(tremuloides_1Mb_gene_new$Gene_num>=85),]

cor.test(tremuloides_1Mb_gene_new_small$Gene_num,tremuloides_1Mb_gene_new_small$tremuloides_ldhat_new,method="spearman")
pcor.test(tremuloides_1Mb_gene_new_small$Gene_num,tremuloides_1Mb_gene_new_small$tremuloides_ldhat_new,tremuloides_1Mb_gene_new_small[,-c(2,4)],method="spearman")

cor.test(tremuloides_1Mb_gene_new_large$Gene_num,tremuloides_1Mb_gene_new_large$tremuloides_ldhat_new,method="spearman")
pcor.test(tremuloides_1Mb_gene_new_large$Gene_num,tremuloides_1Mb_gene_new_large$tremuloides_ldhat_new,tremuloides_1Mb_gene_new_large[,-c(2,4)],method="spearman")


cor.test(tremuloides_1Mb_gene_new$Gene_num,tremuloides_1Mb_gene_new$tremuloides_ldhat_new,method="spearman")
cor.test(tremuloides_1Mb_gene_new$GC,tremuloides_1Mb_gene_new$tremuloides_ldhat_new,method="spearman")
cor.test(tremuloides_1Mb_gene_new$Thetas_tremuloides_numsites,tremuloides_1Mb_gene_new$tremuloides_ldhat_new,method="spearman")

cor.test(tremuloides_1Mb_gene_new$Gene_num,tremuloides_1Mb_gene_new$GC,method="spearman")
cor.test(tremuloides_1Mb_gene_new$Gene_num,tremuloides_1Mb_gene_new$Thetas_tremuloides_numsites,method="spearman")

cor.test(tremuloides_1Mb_gene_new$GC,tremuloides_1Mb_gene_new$Thetas_tremuloides_numsites,method="spearman")

###trichocarpa
trichocarpa_1Mb_gene=total_1Mb[,c("Thetas_trichocarpa_numsites","trichocarpa_ldhat_new","GC","Gene_num")]
trichocarpa_1Mb_gene_new=trichocarpa_1Mb_gene[complete.cases(trichocarpa_1Mb_gene),]
cor.test(trichocarpa_1Mb_gene_new$Gene_num,trichocarpa_1Mb_gene_new$trichocarpa_ldhat_new)
cor.test(trichocarpa_1Mb_gene_new$Gene_num,trichocarpa_1Mb_gene_new$trichocarpa_ldhat_new,method="spearman")
pcor.test(trichocarpa_1Mb_gene_new$Gene_num,trichocarpa_1Mb_gene_new$trichocarpa_ldhat_new,trichocarpa_1Mb_gene_new[,-c(2,4)],method="spearman")
pcor.test(trichocarpa_1Mb_gene_new$Gene_num,trichocarpa_1Mb_gene_new$trichocarpa_ldhat_new,trichocarpa_1Mb_gene_new[,-c(2,4)])

###small and large gene number
trichocarpa_1Mb_gene_new_small=trichocarpa_1Mb_gene_new[which(trichocarpa_1Mb_gene_new$Gene_num<85),]
trichocarpa_1Mb_gene_new_large=trichocarpa_1Mb_gene_new[which(trichocarpa_1Mb_gene_new$Gene_num>=85),]

cor.test(trichocarpa_1Mb_gene_new_small$Gene_num,trichocarpa_1Mb_gene_new_small$trichocarpa_ldhat_new,method="spearman")
pcor.test(trichocarpa_1Mb_gene_new_small$Gene_num,trichocarpa_1Mb_gene_new_small$trichocarpa_ldhat_new,trichocarpa_1Mb_gene_new_small[,-c(2,4)],method="spearman")

cor.test(trichocarpa_1Mb_gene_new_large$Gene_num,trichocarpa_1Mb_gene_new_large$trichocarpa_ldhat_new,method="spearman")
pcor.test(trichocarpa_1Mb_gene_new_large$Gene_num,trichocarpa_1Mb_gene_new_large$trichocarpa_ldhat_new,trichocarpa_1Mb_gene_new_large[,-c(2,4)],method="spearman")


cor.test(trichocarpa_1Mb_gene_new$Gene_num,trichocarpa_1Mb_gene_new$trichocarpa_ldhat_new,method="spearman")
cor.test(trichocarpa_1Mb_gene_new$GC,trichocarpa_1Mb_gene_new$trichocarpa_ldhat_new,method="spearman")
cor.test(trichocarpa_1Mb_gene_new$Thetas_trichocarpa_numsites,trichocarpa_1Mb_gene_new$trichocarpa_ldhat_new,method="spearman")

cor.test(trichocarpa_1Mb_gene_new$Gene_num,trichocarpa_1Mb_gene_new$GC,method="spearman")
cor.test(trichocarpa_1Mb_gene_new$Gene_num,trichocarpa_1Mb_gene_new$Thetas_trichocarpa_numsites,method="spearman")

cor.test(trichocarpa_1Mb_gene_new$GC,trichocarpa_1Mb_gene_new$Thetas_trichocarpa_numsites,method="spearman")



####correlation between d0-fold and pintergenic
#tremula
cor.test(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremula_four_fold_tP,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremula_four_fold_tajD,method="spearman")

cor.test(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremula_intergenic_tP,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremula_intergenic_tajD,method="spearman")
#tremuloides
cor.test(total_1Mb$tremuloides_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremuloides_four_fold_tP,method="spearman")
cor.test(total_1Mb$tremuloides_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremuloides_four_fold_tajD,method="spearman")

cor.test(total_1Mb$tremuloides_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremuloides_intergenic_tP,method="spearman")
cor.test(total_1Mb$tremuloides_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremuloides_intergenic_tajD,method="spearman")

#trichocarpa
cor.test(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_trichocarpa_four_fold_tP,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_trichocarpa_four_fold_tajD,method="spearman")

cor.test(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_trichocarpa_intergenic_tP,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_trichocarpa_intergenic_tajD,method="spearman")




#####correlation between d0-fold and pi4-fold
###tremula
tremula_1Mb_zero_fold=total_1Mb[,c("Thetas_tremula_four_fold_tP","Thetas_tremula_four_fold_numsites","tremula_ldhat_new","GC","Gene_num","tremula_trichocarpa_four_fold_fixed","tremula_trichocarpa_zero_fold_fixed","Thetas_tremula_zero_fold_numsites")]
tremula_1Mb_zero_fold_new=tremula_1Mb_zero_fold[complete.cases(tremula_1Mb_zero_fold),]

cor.test(tremula_1Mb_zero_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,method="spearman")
pcor.test(tremula_1Mb_zero_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,tremula_1Mb_zero_fold_new[,-c(1,7)],method="spearman")

###tremuloides
tremuloides_1Mb_zero_fold=total_1Mb[,c("Thetas_tremuloides_four_fold_tP","Thetas_tremuloides_four_fold_numsites","tremuloides_ldhat_new","GC","Gene_num","tremuloides_trichocarpa_four_fold_fixed","tremuloides_trichocarpa_zero_fold_fixed","Thetas_tremuloides_zero_fold_numsites")]
tremuloides_1Mb_zero_fold_new=tremuloides_1Mb_zero_fold[complete.cases(tremuloides_1Mb_zero_fold),]

cor.test(tremuloides_1Mb_zero_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_zero_fold_new$tremuloides_trichocarpa_zero_fold_fixed,method="spearman")
pcor.test(tremuloides_1Mb_zero_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_zero_fold_new$tremuloides_trichocarpa_zero_fold_fixed,tremuloides_1Mb_zero_fold_new[,-c(1,7)],method="spearman")

###trichocarpa
trichocarpa_1Mb_zero_fold=total_1Mb[,c("Thetas_trichocarpa_four_fold_tP","Thetas_trichocarpa_four_fold_numsites","trichocarpa_ldhat_new","GC","Gene_num","tremula_trichocarpa_four_fold_fixed","tremula_trichocarpa_zero_fold_fixed","Thetas_trichocarpa_zero_fold_numsites")]
trichocarpa_1Mb_zero_fold_new=trichocarpa_1Mb_zero_fold[complete.cases(trichocarpa_1Mb_zero_fold),]

cor.test(trichocarpa_1Mb_zero_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,method="spearman")
pcor.test(trichocarpa_1Mb_zero_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,trichocarpa_1Mb_zero_fold_new[,-c(1,7)],method="spearman")


#####correlation between d0-fold and piintergenic
###tremula
tremula_1Mb_zero_fold=total_1Mb[,c("Thetas_tremula_intergenic_tP","Thetas_tremula_intergenic_numsites","tremula_ldhat_new","GC","Gene_num","tremula_trichocarpa_intergenic_fixed","tremula_trichocarpa_zero_fold_fixed","Thetas_tremula_zero_fold_numsites")]
tremula_1Mb_zero_fold_new=tremula_1Mb_zero_fold[complete.cases(tremula_1Mb_zero_fold),]

cor.test(tremula_1Mb_zero_fold_new$Thetas_tremula_intergenic_tP,tremula_1Mb_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,method="spearman")
pcor.test(tremula_1Mb_zero_fold_new$Thetas_tremula_intergenic_tP,tremula_1Mb_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,tremula_1Mb_zero_fold_new[,-c(1,7)],method="spearman")

###tremuloides
tremuloides_1Mb_zero_fold=total_1Mb[,c("Thetas_tremuloides_intergenic_tP","Thetas_tremuloides_intergenic_numsites","tremuloides_ldhat_new","GC","Gene_num","tremuloides_trichocarpa_intergenic_fixed","tremuloides_trichocarpa_zero_fold_fixed","Thetas_tremuloides_zero_fold_numsites")]
tremuloides_1Mb_zero_fold_new=tremuloides_1Mb_zero_fold[complete.cases(tremuloides_1Mb_zero_fold),]

cor.test(tremuloides_1Mb_zero_fold_new$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_zero_fold_new$tremuloides_trichocarpa_zero_fold_fixed,method="spearman")
pcor.test(tremuloides_1Mb_zero_fold_new$Thetas_tremuloides_intergenic_tP,tremuloides_1Mb_zero_fold_new$tremuloides_trichocarpa_zero_fold_fixed,tremuloides_1Mb_zero_fold_new[,-c(1,7)],method="spearman")

###tremuloides
trichocarpa_1Mb_zero_fold=total_1Mb[,c("Thetas_trichocarpa_intergenic_tP","Thetas_trichocarpa_intergenic_numsites","trichocarpa_ldhat_new","GC","Gene_num","tremula_trichocarpa_intergenic_fixed","tremula_trichocarpa_zero_fold_fixed","Thetas_trichocarpa_zero_fold_numsites")]
trichocarpa_1Mb_zero_fold_new=trichocarpa_1Mb_zero_fold[complete.cases(trichocarpa_1Mb_zero_fold),]

cor.test(trichocarpa_1Mb_zero_fold_new$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,method="spearman")
pcor.test(trichocarpa_1Mb_zero_fold_new$Thetas_trichocarpa_intergenic_tP,trichocarpa_1Mb_zero_fold_new$tremula_trichocarpa_zero_fold_fixed,trichocarpa_1Mb_zero_fold_new[,-c(1,7)],method="spearman")


###correlation between recombiantion rate adn pi_0-fold/pi_4-fold

###tremula
tremula_1Mb_zero_four_fold=total_1Mb[,c("Thetas_tremula_four_fold_tP","Thetas_tremula_zero_fold_tP","Thetas_tremula_four_fold_numsites","tremula_ldhat_new","GC","Gene_num","Thetas_tremula_zero_fold_numsites")]
tremula_1Mb_zero_four_fold_new=tremula_1Mb_zero_four_fold[complete.cases(tremula_1Mb_zero_four_fold),]

cor.test(tremula_1Mb_zero_four_fold_new$Thetas_tremula_zero_fold_tP/tremula_1Mb_zero_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_zero_four_fold_new$tremula_ldhat_new,method="spearman")
cor.test(tremula_1Mb_zero_four_fold_new$Thetas_tremula_zero_fold_tP/tremula_1Mb_zero_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_zero_four_fold_new$tremula_ldhat_new)

pcor.test(tremula_1Mb_zero_four_fold_new$Thetas_tremula_zero_fold_tP/tremula_1Mb_zero_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_zero_four_fold_new$tremula_ldhat_new,tremula_1Mb_zero_four_fold_new[,-c(1,2,4)],method="spearman")


###tremuloides
tremuloides_1Mb_zero_four_fold=total_1Mb[,c("Thetas_tremuloides_four_fold_tP","Thetas_tremuloides_zero_fold_tP","Thetas_tremuloides_four_fold_numsites","tremuloides_ldhat_new","GC","Gene_num","Thetas_tremuloides_zero_fold_numsites")]
tremuloides_1Mb_zero_four_fold_new=tremuloides_1Mb_zero_four_fold[complete.cases(tremuloides_1Mb_zero_four_fold),]

cor.test(tremuloides_1Mb_zero_four_fold_new$Thetas_tremuloides_zero_fold_tP/tremuloides_1Mb_zero_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_zero_four_fold_new$tremuloides_ldhat_new,method="spearman")
cor.test(tremuloides_1Mb_zero_four_fold_new$Thetas_tremuloides_zero_fold_tP/tremuloides_1Mb_zero_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_zero_four_fold_new$tremuloides_ldhat_new)
pcor.test(tremuloides_1Mb_zero_four_fold_new$Thetas_tremuloides_zero_fold_tP/tremuloides_1Mb_zero_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_zero_four_fold_new$tremuloides_ldhat_new,tremuloides_1Mb_zero_four_fold_new[,-c(1,2,4)],method="spearman")

###trichocarpa
trichocarpa_1Mb_zero_four_fold=total_1Mb[,c("Thetas_trichocarpa_four_fold_tP","Thetas_trichocarpa_zero_fold_tP","Thetas_trichocarpa_four_fold_numsites","trichocarpa_ldhat_new","GC","Gene_num","Thetas_trichocarpa_zero_fold_numsites")]
trichocarpa_1Mb_zero_four_fold_new=trichocarpa_1Mb_zero_four_fold[complete.cases(trichocarpa_1Mb_zero_four_fold),]

cor.test(trichocarpa_1Mb_zero_four_fold_new$Thetas_trichocarpa_zero_fold_tP/trichocarpa_1Mb_zero_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_zero_four_fold_new$trichocarpa_ldhat_new,method="spearman")
cor.test(trichocarpa_1Mb_zero_four_fold_new$Thetas_trichocarpa_zero_fold_tP/trichocarpa_1Mb_zero_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_zero_four_fold_new$trichocarpa_ldhat_new)
pcor.test(trichocarpa_1Mb_zero_four_fold_new$Thetas_trichocarpa_zero_fold_tP/trichocarpa_1Mb_zero_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_zero_four_fold_new$trichocarpa_ldhat_new,trichocarpa_1Mb_zero_four_fold_new[,-c(1,2,4)],method="spearman")

###correlation between recombiantion rate adn d_0-fold/d_4-fold

###tremula
tremula_1Mb_zero_four_fold=total_1Mb[,c("tremula_trichocarpa_four_fold_fixed","tremula_trichocarpa_zero_fold_fixed","Thetas_tremula_four_fold_numsites","tremula_ldhat_new","GC","Gene_num","Thetas_tremula_zero_fold_numsites")]
tremula_1Mb_zero_four_fold_new=tremula_1Mb_zero_four_fold[complete.cases(tremula_1Mb_zero_four_fold),]

cor.test(tremula_1Mb_zero_four_fold_new$tremula_trichocarpa_zero_fold_fixed/tremula_1Mb_zero_four_fold_new$tremula_trichocarpa_four_fold_fixed,tremula_1Mb_zero_four_fold_new$tremula_ldhat_new,method="spearman")
cor.test(tremula_1Mb_zero_four_fold_new$tremula_trichocarpa_zero_fold_fixed/tremula_1Mb_zero_four_fold_new$tremula_trichocarpa_four_fold_fixed,tremula_1Mb_zero_four_fold_new$tremula_ldhat_new)

pcor.test(tremula_1Mb_zero_four_fold_new$tremula_trichocarpa_zero_fold_fixed/tremula_1Mb_zero_four_fold_new$tremula_trichocarpa_four_fold_fixed,tremula_1Mb_zero_four_fold_new$tremula_ldhat_new,tremula_1Mb_zero_four_fold_new[,-c(1,2,4)],method="spearman")


###tremuloides
tremuloides_1Mb_zero_four_fold=total_1Mb[,c("tremuloides_trichocarpa_four_fold_fixed","tremuloides_trichocarpa_zero_fold_fixed","Thetas_tremuloides_four_fold_numsites","tremuloides_ldhat_new","GC","Gene_num","Thetas_tremuloides_zero_fold_numsites")]
tremuloides_1Mb_zero_four_fold_new=tremuloides_1Mb_zero_four_fold[complete.cases(tremuloides_1Mb_zero_four_fold),]

cor.test(tremuloides_1Mb_zero_four_fold_new$tremuloides_trichocarpa_zero_fold_fixed/tremuloides_1Mb_zero_four_fold_new$tremuloides_trichocarpa_four_fold_fixed,tremuloides_1Mb_zero_four_fold_new$tremuloides_ldhat_new,method="spearman")
cor.test(tremuloides_1Mb_zero_four_fold_new$Thetas_tremuloides_zero_fold_tP/tremuloides_1Mb_zero_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_zero_four_fold_new$tremuloides_ldhat_new)
pcor.test(tremuloides_1Mb_zero_four_fold_new$tremuloides_trichocarpa_zero_fold_fixed/tremuloides_1Mb_zero_four_fold_new$tremuloides_trichocarpa_four_fold_fixed,tremuloides_1Mb_zero_four_fold_new$tremuloides_ldhat_new,tremuloides_1Mb_zero_four_fold_new[,-c(1,2,4)],method="spearman")

###trichocarpa
trichocarpa_1Mb_zero_four_fold=total_1Mb[,c("tremula_trichocarpa_four_fold_fixed","tremula_trichocarpa_zero_fold_fixed","Thetas_trichocarpa_four_fold_numsites","trichocarpa_ldhat_new","GC","Gene_num","Thetas_trichocarpa_zero_fold_numsites")]
trichocarpa_1Mb_zero_four_fold_new=trichocarpa_1Mb_zero_four_fold[complete.cases(trichocarpa_1Mb_zero_four_fold),]

cor.test(trichocarpa_1Mb_zero_four_fold_new$tremula_trichocarpa_zero_fold_fixed/trichocarpa_1Mb_zero_four_fold_new$tremula_trichocarpa_four_fold_fixed,trichocarpa_1Mb_zero_four_fold_new$trichocarpa_ldhat_new,method="spearman")
cor.test(trichocarpa_1Mb_zero_four_fold_new$Thetas_trichocarpa_zero_fold_tP/trichocarpa_1Mb_zero_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_zero_four_fold_new$trichocarpa_ldhat_new)
pcor.test(trichocarpa_1Mb_zero_four_fold_new$tremula_trichocarpa_zero_fold_fixed/trichocarpa_1Mb_zero_four_fold_new$tremula_trichocarpa_four_fold_fixed,trichocarpa_1Mb_zero_four_fold_new$trichocarpa_ldhat_new,trichocarpa_1Mb_zero_four_fold_new[,-c(1,2,4)],method="spearman")


###correlation between recombiantion rate adn pi_intron/pi_4-fold

###tremula
tremula_1Mb_intron_four_fold=total_1Mb[,c("Thetas_tremula_four_fold_tP","Thetas_tremula_intron_tP","Thetas_tremula_four_fold_numsites","tremula_ldhat_new","GC","Gene_num","Thetas_tremula_intron_numsites")]
tremula_1Mb_intron_four_fold_new=tremula_1Mb_intron_four_fold[complete.cases(tremula_1Mb_intron_four_fold),]

cor.test(tremula_1Mb_intron_four_fold_new$Thetas_tremula_intron_tP/tremula_1Mb_intron_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_intron_four_fold_new$tremula_ldhat_new,method="spearman")
cor.test(tremula_1Mb_intron_four_fold_new$Thetas_tremula_intron_tP/tremula_1Mb_intron_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_intron_four_fold_new$tremula_ldhat_new)

pcor.test(tremula_1Mb_intron_four_fold_new$Thetas_tremula_intron_tP/tremula_1Mb_intron_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_intron_four_fold_new$tremula_ldhat_new,tremula_1Mb_intron_four_fold_new[,-c(1,2,4)],method="spearman")


###tremuloides
tremuloides_1Mb_intron_four_fold=total_1Mb[,c("Thetas_tremuloides_four_fold_tP","Thetas_tremuloides_intron_tP","Thetas_tremuloides_four_fold_numsites","tremuloides_ldhat_new","GC","Gene_num","Thetas_tremuloides_intron_numsites")]
tremuloides_1Mb_intron_four_fold_new=tremuloides_1Mb_intron_four_fold[complete.cases(tremuloides_1Mb_intron_four_fold),]

cor.test(tremuloides_1Mb_intron_four_fold_new$Thetas_tremuloides_intron_tP/tremuloides_1Mb_intron_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_intron_four_fold_new$tremuloides_ldhat_new,method="spearman")
cor.test(tremuloides_1Mb_intron_four_fold_new$Thetas_tremuloides_intron_tP/tremuloides_1Mb_intron_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_intron_four_fold_new$tremuloides_ldhat_new)
pcor.test(tremuloides_1Mb_intron_four_fold_new$Thetas_tremuloides_intron_tP/tremuloides_1Mb_intron_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_intron_four_fold_new$tremuloides_ldhat_new,tremuloides_1Mb_intron_four_fold_new[,-c(1,2,4)],method="spearman")

###trichocarpa
trichocarpa_1Mb_intron_four_fold=total_1Mb[,c("Thetas_trichocarpa_four_fold_tP","Thetas_trichocarpa_intron_tP","Thetas_trichocarpa_four_fold_numsites","trichocarpa_ldhat_new","GC","Gene_num","Thetas_trichocarpa_intron_numsites")]
trichocarpa_1Mb_intron_four_fold_new=trichocarpa_1Mb_intron_four_fold[complete.cases(trichocarpa_1Mb_intron_four_fold),]

cor.test(trichocarpa_1Mb_intron_four_fold_new$Thetas_trichocarpa_intron_tP/trichocarpa_1Mb_intron_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_intron_four_fold_new$trichocarpa_ldhat_new,method="spearman")
cor.test(trichocarpa_1Mb_intron_four_fold_new$Thetas_trichocarpa_intron_tP/trichocarpa_1Mb_intron_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_intron_four_fold_new$trichocarpa_ldhat_new)
pcor.test(trichocarpa_1Mb_intron_four_fold_new$Thetas_trichocarpa_intron_tP/trichocarpa_1Mb_intron_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_intron_four_fold_new$trichocarpa_ldhat_new,trichocarpa_1Mb_intron_four_fold_new[,-c(1,2,4)],method="spearman")


###correlation between recombiantion rate adn pi_3UTR/pi_4-fold

###tremula
tremula_1Mb_3UTR_four_fold=total_1Mb[,c("Thetas_tremula_four_fold_tP","Thetas_tremula_3UTR_tP","Thetas_tremula_four_fold_numsites","tremula_ldhat_new","GC","Gene_num","Thetas_tremula_3UTR_numsites")]
tremula_1Mb_3UTR_four_fold_new=tremula_1Mb_3UTR_four_fold[complete.cases(tremula_1Mb_3UTR_four_fold),]

cor.test(tremula_1Mb_3UTR_four_fold_new$Thetas_tremula_3UTR_tP/tremula_1Mb_3UTR_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_3UTR_four_fold_new$tremula_ldhat_new,method="spearman")
cor.test(tremula_1Mb_3UTR_four_fold_new$Thetas_tremula_3UTR_tP/tremula_1Mb_3UTR_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_3UTR_four_fold_new$tremula_ldhat_new)

pcor.test(tremula_1Mb_3UTR_four_fold_new$Thetas_tremula_3UTR_tP/tremula_1Mb_3UTR_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_3UTR_four_fold_new$tremula_ldhat_new,tremula_1Mb_3UTR_four_fold_new[,-c(1,2,4)],method="spearman")


###tremuloides
tremuloides_1Mb_3UTR_four_fold=total_1Mb[,c("Thetas_tremuloides_four_fold_tP","Thetas_tremuloides_3UTR_tP","Thetas_tremuloides_four_fold_numsites","tremuloides_ldhat_new","GC","Gene_num","Thetas_tremuloides_3UTR_numsites")]
tremuloides_1Mb_3UTR_four_fold_new=tremuloides_1Mb_3UTR_four_fold[complete.cases(tremuloides_1Mb_3UTR_four_fold),]

cor.test(tremuloides_1Mb_3UTR_four_fold_new$Thetas_tremuloides_3UTR_tP/tremuloides_1Mb_3UTR_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_3UTR_four_fold_new$tremuloides_ldhat_new,method="spearman")
cor.test(tremuloides_1Mb_3UTR_four_fold_new$Thetas_tremuloides_3UTR_tP/tremuloides_1Mb_3UTR_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_3UTR_four_fold_new$tremuloides_ldhat_new)
pcor.test(tremuloides_1Mb_3UTR_four_fold_new$Thetas_tremuloides_3UTR_tP/tremuloides_1Mb_3UTR_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_3UTR_four_fold_new$tremuloides_ldhat_new,tremuloides_1Mb_3UTR_four_fold_new[,-c(1,2,4)],method="spearman")

###trichocarpa
trichocarpa_1Mb_3UTR_four_fold=total_1Mb[,c("Thetas_trichocarpa_four_fold_tP","Thetas_trichocarpa_3UTR_tP","Thetas_trichocarpa_four_fold_numsites","trichocarpa_ldhat_new","GC","Gene_num","Thetas_trichocarpa_3UTR_numsites")]
trichocarpa_1Mb_3UTR_four_fold_new=trichocarpa_1Mb_3UTR_four_fold[complete.cases(trichocarpa_1Mb_3UTR_four_fold),]

cor.test(trichocarpa_1Mb_3UTR_four_fold_new$Thetas_trichocarpa_3UTR_tP/trichocarpa_1Mb_3UTR_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_3UTR_four_fold_new$trichocarpa_ldhat_new,method="spearman")
cor.test(trichocarpa_1Mb_3UTR_four_fold_new$Thetas_trichocarpa_3UTR_tP/trichocarpa_1Mb_3UTR_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_3UTR_four_fold_new$trichocarpa_ldhat_new)
pcor.test(trichocarpa_1Mb_3UTR_four_fold_new$Thetas_trichocarpa_3UTR_tP/trichocarpa_1Mb_3UTR_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_3UTR_four_fold_new$trichocarpa_ldhat_new,trichocarpa_1Mb_3UTR_four_fold_new[,-c(1,2,4)],method="spearman")


###correlation between recombiantion rate adn pi_5UTR/pi_4-fold

###tremula
tremula_1Mb_5UTR_four_fold=total_1Mb[,c("Thetas_tremula_four_fold_tP","Thetas_tremula_5UTR_tP","Thetas_tremula_four_fold_numsites","tremula_ldhat_new","GC","Gene_num","Thetas_tremula_5UTR_numsites")]
tremula_1Mb_5UTR_four_fold_new=tremula_1Mb_5UTR_four_fold[complete.cases(tremula_1Mb_5UTR_four_fold),]

cor.test(tremula_1Mb_5UTR_four_fold_new$Thetas_tremula_5UTR_tP/tremula_1Mb_5UTR_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_5UTR_four_fold_new$tremula_ldhat_new,method="spearman")
cor.test(tremula_1Mb_5UTR_four_fold_new$Thetas_tremula_5UTR_tP/tremula_1Mb_5UTR_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_5UTR_four_fold_new$tremula_ldhat_new)

pcor.test(tremula_1Mb_5UTR_four_fold_new$Thetas_tremula_5UTR_tP/tremula_1Mb_5UTR_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_5UTR_four_fold_new$tremula_ldhat_new,tremula_1Mb_5UTR_four_fold_new[,-c(1,2,4)],method="spearman")


###tremuloides
tremuloides_1Mb_5UTR_four_fold=total_1Mb[,c("Thetas_tremuloides_four_fold_tP","Thetas_tremuloides_5UTR_tP","Thetas_tremuloides_four_fold_numsites","tremuloides_ldhat_new","GC","Gene_num","Thetas_tremuloides_5UTR_numsites")]
tremuloides_1Mb_5UTR_four_fold_new=tremuloides_1Mb_5UTR_four_fold[complete.cases(tremuloides_1Mb_5UTR_four_fold),]

cor.test(tremuloides_1Mb_5UTR_four_fold_new$Thetas_tremuloides_5UTR_tP/tremuloides_1Mb_5UTR_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_5UTR_four_fold_new$tremuloides_ldhat_new,method="spearman")
cor.test(tremuloides_1Mb_5UTR_four_fold_new$Thetas_tremuloides_5UTR_tP/tremuloides_1Mb_5UTR_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_5UTR_four_fold_new$tremuloides_ldhat_new)
pcor.test(tremuloides_1Mb_5UTR_four_fold_new$Thetas_tremuloides_5UTR_tP/tremuloides_1Mb_5UTR_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_5UTR_four_fold_new$tremuloides_ldhat_new,tremuloides_1Mb_5UTR_four_fold_new[,-c(1,2,4)],method="spearman")

###trichocarpa
trichocarpa_1Mb_5UTR_four_fold=total_1Mb[,c("Thetas_trichocarpa_four_fold_tP","Thetas_trichocarpa_5UTR_tP","Thetas_trichocarpa_four_fold_numsites","trichocarpa_ldhat_new","GC","Gene_num","Thetas_trichocarpa_5UTR_numsites")]
trichocarpa_1Mb_5UTR_four_fold_new=trichocarpa_1Mb_5UTR_four_fold[complete.cases(trichocarpa_1Mb_5UTR_four_fold),]

cor.test(trichocarpa_1Mb_5UTR_four_fold_new$Thetas_trichocarpa_5UTR_tP/trichocarpa_1Mb_5UTR_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_5UTR_four_fold_new$trichocarpa_ldhat_new,method="spearman")
cor.test(trichocarpa_1Mb_5UTR_four_fold_new$Thetas_trichocarpa_5UTR_tP/trichocarpa_1Mb_5UTR_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_5UTR_four_fold_new$trichocarpa_ldhat_new)
pcor.test(trichocarpa_1Mb_5UTR_four_fold_new$Thetas_trichocarpa_5UTR_tP/trichocarpa_1Mb_5UTR_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_5UTR_four_fold_new$trichocarpa_ldhat_new,trichocarpa_1Mb_5UTR_four_fold_new[,-c(1,2,4)],method="spearman")


###correlation between recombiantion rate adn pi_Intergenic/pi_4-fold

###tremula
tremula_1Mb_intergenic_four_fold=total_1Mb[,c("Thetas_tremula_four_fold_tP","Thetas_tremula_intergenic_tP","Thetas_tremula_four_fold_numsites","tremula_ldhat_new","GC","Gene_num","Thetas_tremula_intergenic_numsites")]
tremula_1Mb_intergenic_four_fold_new=tremula_1Mb_intergenic_four_fold[complete.cases(tremula_1Mb_intergenic_four_fold),]

cor.test(tremula_1Mb_intergenic_four_fold_new$Thetas_tremula_intergenic_tP/tremula_1Mb_intergenic_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_intergenic_four_fold_new$tremula_ldhat_new,method="spearman")
cor.test(tremula_1Mb_intergenic_four_fold_new$Thetas_tremula_intergenic_tP/tremula_1Mb_intergenic_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_intergenic_four_fold_new$tremula_ldhat_new)

pcor.test(tremula_1Mb_intergenic_four_fold_new$Thetas_tremula_intergenic_tP/tremula_1Mb_intergenic_four_fold_new$Thetas_tremula_four_fold_tP,tremula_1Mb_intergenic_four_fold_new$tremula_ldhat_new,tremula_1Mb_intergenic_four_fold_new[,-c(1,2,4)],method="spearman")


###tremuloides
tremuloides_1Mb_intergenic_four_fold=total_1Mb[,c("Thetas_tremuloides_four_fold_tP","Thetas_tremuloides_intergenic_tP","Thetas_tremuloides_four_fold_numsites","tremuloides_ldhat_new","GC","Gene_num","Thetas_tremuloides_intergenic_numsites")]
tremuloides_1Mb_intergenic_four_fold_new=tremuloides_1Mb_intergenic_four_fold[complete.cases(tremuloides_1Mb_intergenic_four_fold),]

cor.test(tremuloides_1Mb_intergenic_four_fold_new$Thetas_tremuloides_intergenic_tP/tremuloides_1Mb_intergenic_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_intergenic_four_fold_new$tremuloides_ldhat_new,method="spearman")
cor.test(tremuloides_1Mb_intergenic_four_fold_new$Thetas_tremuloides_intergenic_tP/tremuloides_1Mb_intergenic_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_intergenic_four_fold_new$tremuloides_ldhat_new)
pcor.test(tremuloides_1Mb_intergenic_four_fold_new$Thetas_tremuloides_intergenic_tP/tremuloides_1Mb_intergenic_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_intergenic_four_fold_new$tremuloides_ldhat_new,tremuloides_1Mb_intergenic_four_fold_new[,-c(1,2,4)],method="spearman")
pcor.test(tremuloides_1Mb_intergenic_four_fold_new$Thetas_tremuloides_intergenic_tP/tremuloides_1Mb_intergenic_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_intergenic_four_fold_new$tremuloides_ldhat_new,tremuloides_1Mb_intergenic_four_fold_new$Gene_num,method="spearman")


###trichocarpa
trichocarpa_1Mb_intergenic_four_fold=total_1Mb[,c("Thetas_trichocarpa_four_fold_tP","Thetas_trichocarpa_intergenic_tP","Thetas_trichocarpa_four_fold_numsites","trichocarpa_ldhat_new","GC","Gene_num","Thetas_trichocarpa_intergenic_numsites")]
trichocarpa_1Mb_intergenic_four_fold_new=trichocarpa_1Mb_intergenic_four_fold[complete.cases(trichocarpa_1Mb_intergenic_four_fold),]

cor.test(trichocarpa_1Mb_intergenic_four_fold_new$Thetas_trichocarpa_intergenic_tP/trichocarpa_1Mb_intergenic_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_intergenic_four_fold_new$trichocarpa_ldhat_new,method="spearman")
cor.test(trichocarpa_1Mb_intergenic_four_fold_new$Thetas_trichocarpa_intergenic_tP/trichocarpa_1Mb_intergenic_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_intergenic_four_fold_new$trichocarpa_ldhat_new)
pcor.test(trichocarpa_1Mb_intergenic_four_fold_new$Thetas_trichocarpa_intergenic_tP/trichocarpa_1Mb_intergenic_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_intergenic_four_fold_new$trichocarpa_ldhat_new,trichocarpa_1Mb_intergenic_four_fold_new[,-c(1,2,4)],method="spearman")
pcor.test(trichocarpa_1Mb_intergenic_four_fold_new$Thetas_trichocarpa_intergenic_tP/trichocarpa_1Mb_intergenic_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_intergenic_four_fold_new$trichocarpa_ldhat_new,trichocarpa_1Mb_intergenic_four_fold_new$Gene_num,method="spearman")


###correlation between recombiantion rate adn d_Intergenic/d_4-fold

###tremula
tremula_1Mb_intergenic_four_fold=total_1Mb[,c("tremula_trichocarpa_intergenic_fixed","tremula_trichocarpa_four_fold_fixed","Thetas_tremula_four_fold_numsites","tremula_ldhat_new","GC","Gene_num","Thetas_tremula_intergenic_numsites")]
tremula_1Mb_intergenic_four_fold_new=tremula_1Mb_intergenic_four_fold[complete.cases(tremula_1Mb_intergenic_four_fold),]

cor.test(tremula_1Mb_intergenic_four_fold_new$tremula_trichocarpa_intergenic_fixed/tremula_1Mb_intergenic_four_fold_new$tremula_trichocarpa_four_fold_fixed,tremula_1Mb_intergenic_four_fold_new$tremula_ldhat_new,method="spearman")
cor.test(tremula_1Mb_intergenic_four_fold_new$tremula_trichocarpa_intergenic_fixed/tremula_1Mb_intergenic_four_fold_new$tremula_trichocarpa_four_fold_fixed,tremula_1Mb_intergenic_four_fold_new$tremula_ldhat_new)

pcor.test(tremula_1Mb_intergenic_four_fold_new$tremula_trichocarpa_intergenic_fixed/tremula_1Mb_intergenic_four_fold_new$tremula_trichocarpa_four_fold_fixed,tremula_1Mb_intergenic_four_fold_new$tremula_ldhat_new,tremula_1Mb_intergenic_four_fold_new[,-c(1,2,4)],method="spearman")


###tremuloides
tremuloides_1Mb_intergenic_four_fold=total_1Mb[,c("tremuloides_trichocarpa_intergenic_fixed","tremuloides_trichocarpa_four_fold_fixed","Thetas_tremuloides_four_fold_numsites","tremuloides_ldhat_new","GC","Gene_num","Thetas_tremuloides_intergenic_numsites")]
tremuloides_1Mb_intergenic_four_fold_new=tremuloides_1Mb_intergenic_four_fold[complete.cases(tremuloides_1Mb_intergenic_four_fold),]

cor.test(tremuloides_1Mb_intergenic_four_fold_new$tremuloides_trichocarpa_intergenic_fixed/tremuloides_1Mb_intergenic_four_fold_new$tremuloides_trichocarpa_four_fold_fixed,tremuloides_1Mb_intergenic_four_fold_new$tremuloides_ldhat_new,method="spearman")
cor.test(tremuloides_1Mb_intergenic_four_fold_new$tremuloides_trichocarpa_intergenic_fixed/tremuloides_1Mb_intergenic_four_fold_new$tremuloides_trichocarpa_four_fold_fixed,tremuloides_1Mb_intergenic_four_fold_new$tremuloides_ldhat_new)
pcor.test(tremuloides_1Mb_intergenic_four_fold_new$tremuloides_trichocarpa_intergenic_fixed/tremuloides_1Mb_intergenic_four_fold_new$tremuloides_trichocarpa_four_fold_fixed,tremuloides_1Mb_intergenic_four_fold_new$tremuloides_ldhat_new,tremuloides_1Mb_intergenic_four_fold_new[,-c(1,2,4)],method="spearman")
pcor.test(tremuloides_1Mb_intergenic_four_fold_new$Thetas_tremuloides_intergenic_tP/tremuloides_1Mb_intergenic_four_fold_new$Thetas_tremuloides_four_fold_tP,tremuloides_1Mb_intergenic_four_fold_new$tremuloides_ldhat_new,tremuloides_1Mb_intergenic_four_fold_new$Gene_num,method="spearman")


###trichocarpa
trichocarpa_1Mb_intergenic_four_fold=total_1Mb[,c("tremula_trichocarpa_intergenic_fixed","tremula_trichocarpa_four_fold_fixed","trichocarpa_ldhat_new","GC","Gene_num","Thetas_trichocarpa_intergenic_numsites")]
trichocarpa_1Mb_intergenic_four_fold_new=trichocarpa_1Mb_intergenic_four_fold[complete.cases(trichocarpa_1Mb_intergenic_four_fold),]

cor.test(trichocarpa_1Mb_intergenic_four_fold_new$tremula_trichocarpa_intergenic_fixed/trichocarpa_1Mb_intergenic_four_fold_new$tremula_trichocarpa_four_fold_fixed,trichocarpa_1Mb_intergenic_four_fold_new$trichocarpa_ldhat_new,method="spearman")
cor.test(trichocarpa_1Mb_intergenic_four_fold_new$Thetas_trichocarpa_intergenic_tP/trichocarpa_1Mb_intergenic_four_fold_new$Thetas_trichocarpa_four_fold_tP,trichocarpa_1Mb_intergenic_four_fold_new$trichocarpa_ldhat_new)
pcor.test(trichocarpa_1Mb_intergenic_four_fold_new$tremula_trichocarpa_intergenic_fixed/trichocarpa_1Mb_intergenic_four_fold_new$tremula_trichocarpa_four_fold_fixed,trichocarpa_1Mb_intergenic_four_fold_new$trichocarpa_ldhat_new,trichocarpa_1Mb_intergenic_four_fold_new[,-c(1,2,4)],method="spearman")
pcor.test(trichocarpa_1Mb_intergenic_four_fold_new$tremula_trichocarpa_intergenic_fixed/trichocarpa_1Mb_intergenic_four_fold_new$tremula_trichocarpa_four_fold_fixed,trichocarpa_1Mb_intergenic_four_fold_new$trichocarpa_ldhat_new,trichocarpa_1Mb_intergenic_four_fold_new$Gene_num,method="spearman")







