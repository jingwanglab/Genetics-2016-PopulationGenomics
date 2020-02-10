#install.packages("pastecs")
library(pastecs)

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

####measure of relative variability: the coefficient of variation
install.packages("raster")
library(raster)
CV=function(value){
  mean_value=mean(value,na.rm=T)
  sd_value=sd(value,na.rm=T)
  (sd_value/mean_value)*100
}

##diveristy
cv(total_1Mb$Thetas_tremula_tP)
cv(total_1Mb$Thetas_tremuloides_tP)
cv(total_1Mb$Thetas_trichocarpa_tP)
##divergence
cv(total_1Mb$tremula_trichocarpa_dxy)
cv(total_1Mb$tremuloides_trichocarpa_dxy)
cv(total_1Mb$tremula_tremuloides_dxy)
###rho
cv(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),na.rm=T)
cv(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),na.rm=T)
cv(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),na.rm=T)

#tremuloides_total
round(quantile(total_1Mb$Thetas_tremuloides_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_tajD,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),c(.025,.5,.975),na.rm=T),digits=4)
#tremuloides_zero_fold
round(quantile(total_1Mb$Thetas_tremuloides_zero_fold_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_zero_fold_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_zero_fold_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#tremuloides_four_fold
round(quantile(total_1Mb$Thetas_tremuloides_four_fold_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_four_fold_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_four_fold_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#tremuloides_intron_fold
round(quantile(total_1Mb$Thetas_tremuloides_intron_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_intron_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_intron_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#tremuloides_5UTR_fold
round(quantile(total_1Mb$Thetas_tremuloides_5UTR_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_5UTR_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_5UTR_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#tremuloides_3UTR_fold
round(quantile(total_1Mb$Thetas_tremuloides_3UTR_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_3UTR_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_3UTR_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#tremuloides_intergenic_fold
round(quantile(total_1Mb$Thetas_tremuloides_intergenic_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_intergenic_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremuloides_intergenic_tajD,c(.025,.5,.975),na.rm=T),digits=4)


#tremula_total
round(quantile(total_1Mb$Thetas_tremula_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_tajD,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),c(.025,.5,.975),na.rm=T),digits=4)
#tremula_zero_fold
round(quantile(total_1Mb$Thetas_tremula_zero_fold_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_zero_fold_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_zero_fold_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#tremula_four_fold
round(quantile(total_1Mb$Thetas_tremula_four_fold_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_four_fold_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_four_fold_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#tremula_intron_fold
round(quantile(total_1Mb$Thetas_tremula_intron_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_intron_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_intron_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#tremula_5UTR_fold
round(quantile(total_1Mb$Thetas_tremula_5UTR_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_5UTR_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_5UTR_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#tremula_3UTR_fold
round(quantile(total_1Mb$Thetas_tremula_3UTR_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_3UTR_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_3UTR_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#tremula_intergenic_fold
round(quantile(total_1Mb$Thetas_tremula_intergenic_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_intergenic_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_tremula_intergenic_tajD,c(.025,.5,.975),na.rm=T),digits=4)


#trichocarpa_total
round(quantile(total_1Mb$Thetas_trichocarpa_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_tajD,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),c(.025,.5,.975),na.rm=T),digits=4)
#trichocarpa_zero_fold
round(quantile(total_1Mb$Thetas_trichocarpa_zero_fold_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_zero_fold_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_zero_fold_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#trichocarpa_four_fold
round(quantile(total_1Mb$Thetas_trichocarpa_four_fold_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_four_fold_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_four_fold_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#trichocarpa_intron_fold
round(quantile(total_1Mb$Thetas_trichocarpa_intron_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_intron_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_intron_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#trichocarpa_5UTR_fold
round(quantile(total_1Mb$Thetas_trichocarpa_5UTR_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_5UTR_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_5UTR_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#trichocarpa_3UTR_fold
round(quantile(total_1Mb$Thetas_trichocarpa_3UTR_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_3UTR_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_3UTR_tajD,c(.025,.5,.975),na.rm=T),digits=4)
#trichocarpa_intergenic_fold
round(quantile(total_1Mb$Thetas_trichocarpa_intergenic_tP,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_intergenic_tW,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$Thetas_trichocarpa_intergenic_tajD,c(.025,.5,.975),na.rm=T),digits=4)


####dxy between P.tremula and P.tremuloides
round(quantile(total_1Mb$tremula_tremuloides_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremula_tremuloides_zero_fold_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremula_tremuloides_four_fold_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremula_tremuloides_intron_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremula_tremuloides_3UTR_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremula_tremuloides_5UTR_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremula_tremuloides_intergenic_dxy,c(.025,.5,.975),na.rm=T),digits=4)

####dxy between P.tremula and P.trichocarpa
round(quantile(total_1Mb$tremula_trichocarpa_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremula_trichocarpa_zero_fold_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremula_trichocarpa_four_fold_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremula_trichocarpa_intron_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremula_trichocarpa_5UTR_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremula_trichocarpa_3UTR_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremula_trichocarpa_intergenic_dxy,c(.025,.5,.975),na.rm=T),digits=4)

####dxy between P.tremuloides and P.trichocarpa
round(quantile(total_1Mb$tremuloides_trichocarpa_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremuloides_trichocarpa_zero_fold_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremuloides_trichocarpa_four_fold_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremuloides_trichocarpa_intron_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremuloides_trichocarpa_5UTR_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremuloides_trichocarpa_3UTR_dxy,c(.025,.5,.975),na.rm=T),digits=4)
round(quantile(total_1Mb$tremuloides_trichocarpa_intergenic_dxy,c(.025,.5,.975),na.rm=T),digits=4)



#####Wilcoxon Sign rank test
###
wilcox.test(total_1Mb$Thetas_tremula_four_fold_tP,total_1Mb$Thetas_tremula_3UTR_tP)
wilcox.test(total_1Mb$Thetas_tremuloides_four_fold_tP,total_1Mb$Thetas_tremuloides_3UTR_tP)
wilcox.test(total_1Mb$Thetas_trichocarpa_four_fold_tP,total_1Mb$Thetas_trichocarpa_3UTR_tP)

wilcox.test(total_1Mb$Thetas_tremula_four_fold_tP,total_1Mb$Thetas_tremula_5UTR_tP)
wilcox.test(total_1Mb$Thetas_tremuloides_four_fold_tP,total_1Mb$Thetas_tremuloides_5UTR_tP)
wilcox.test(total_1Mb$Thetas_trichocarpa_four_fold_tP,total_1Mb$Thetas_trichocarpa_5UTR_tP)


###correlation of tP and tajD with species

##tremula-4fold
cor.test(total_1Mb$Thetas_tremula_four_fold_tP,total_1Mb$Thetas_tremula_zero_fold_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_four_fold_tP,total_1Mb$Thetas_tremula_3UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_four_fold_tP,total_1Mb$Thetas_tremula_5UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_four_fold_tP,total_1Mb$Thetas_tremula_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_four_fold_tP,total_1Mb$Thetas_tremula_intergenic_tP,method="spearman")
##tremula_0-fold
cor.test(total_1Mb$Thetas_tremula_zero_fold_tP,total_1Mb$Thetas_tremula_3UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_zero_fold_tP,total_1Mb$Thetas_tremula_5UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_zero_fold_tP,total_1Mb$Thetas_tremula_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_zero_fold_tP,total_1Mb$Thetas_tremula_intergenic_tP,method="spearman")
###tremula-3UTR
cor.test(total_1Mb$Thetas_tremula_3UTR_tP,total_1Mb$Thetas_tremula_5UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_3UTR_tP,total_1Mb$Thetas_tremula_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_3UTR_tP,total_1Mb$Thetas_tremula_intergenic_tP,method="spearman")
###tremula-5UTR
cor.test(total_1Mb$Thetas_tremula_5UTR_tP,total_1Mb$Thetas_tremula_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_5UTR_tP,total_1Mb$Thetas_tremula_intergenic_tP,method="spearman")
###tremula-intron
cor.test(total_1Mb$Thetas_tremula_intron_tP,total_1Mb$Thetas_tremula_intergenic_tP,method="spearman")


###tremuloides
##tremuloides-4fold
cor.test(total_1Mb$Thetas_tremuloides_four_fold_tP,total_1Mb$Thetas_tremuloides_zero_fold_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_four_fold_tP,total_1Mb$Thetas_tremuloides_3UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_four_fold_tP,total_1Mb$Thetas_tremuloides_5UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_four_fold_tP,total_1Mb$Thetas_tremuloides_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_four_fold_tP,total_1Mb$Thetas_tremuloides_intergenic_tP,method="spearman")
##tremuloides_0-fold
cor.test(total_1Mb$Thetas_tremuloides_zero_fold_tP,total_1Mb$Thetas_tremuloides_3UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_zero_fold_tP,total_1Mb$Thetas_tremuloides_5UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_zero_fold_tP,total_1Mb$Thetas_tremuloides_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_zero_fold_tP,total_1Mb$Thetas_tremuloides_intergenic_tP,method="spearman")
###tremuloides-3UTR
cor.test(total_1Mb$Thetas_tremuloides_3UTR_tP,total_1Mb$Thetas_tremuloides_5UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_3UTR_tP,total_1Mb$Thetas_tremuloides_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_3UTR_tP,total_1Mb$Thetas_tremuloides_intergenic_tP,method="spearman")
###tremuloides-5UTR
cor.test(total_1Mb$Thetas_tremuloides_5UTR_tP,total_1Mb$Thetas_tremuloides_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_5UTR_tP,total_1Mb$Thetas_tremuloides_intergenic_tP,method="spearman")
###tremuloides-intron
cor.test(total_1Mb$Thetas_tremuloides_intron_tP,total_1Mb$Thetas_tremuloides_intergenic_tP,method="spearman")


###trichocarpa
##trichocarpa-4fold
cor.test(total_1Mb$Thetas_trichocarpa_four_fold_tP,total_1Mb$Thetas_trichocarpa_zero_fold_tP,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_four_fold_tP,total_1Mb$Thetas_trichocarpa_3UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_four_fold_tP,total_1Mb$Thetas_trichocarpa_5UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_four_fold_tP,total_1Mb$Thetas_trichocarpa_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_four_fold_tP,total_1Mb$Thetas_trichocarpa_intergenic_tP,method="spearman")
##trichocarpa_0-fold
cor.test(total_1Mb$Thetas_trichocarpa_zero_fold_tP,total_1Mb$Thetas_trichocarpa_3UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_zero_fold_tP,total_1Mb$Thetas_trichocarpa_5UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_zero_fold_tP,total_1Mb$Thetas_trichocarpa_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_zero_fold_tP,total_1Mb$Thetas_trichocarpa_intergenic_tP,method="spearman")
###trichocarpa-3UTR
cor.test(total_1Mb$Thetas_trichocarpa_3UTR_tP,total_1Mb$Thetas_trichocarpa_5UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_3UTR_tP,total_1Mb$Thetas_trichocarpa_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_3UTR_tP,total_1Mb$Thetas_trichocarpa_intergenic_tP,method="spearman")
###trichocarpa-5UTR
cor.test(total_1Mb$Thetas_trichocarpa_5UTR_tP,total_1Mb$Thetas_trichocarpa_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_5UTR_tP,total_1Mb$Thetas_trichocarpa_intergenic_tP,method="spearman")
###trichocarpa-intron
cor.test(total_1Mb$Thetas_trichocarpa_intron_tP,total_1Mb$Thetas_trichocarpa_intergenic_tP,method="spearman")


###tajD
##tremula-4fold
cor.test(total_1Mb$Thetas_tremula_four_fold_tajD,total_1Mb$Thetas_tremula_zero_fold_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_four_fold_tajD,total_1Mb$Thetas_tremula_3UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_four_fold_tajD,total_1Mb$Thetas_tremula_5UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_four_fold_tajD,total_1Mb$Thetas_tremula_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_four_fold_tajD,total_1Mb$Thetas_tremula_intergenic_tajD,method="spearman")
##tremula_0-fold
cor.test(total_1Mb$Thetas_tremula_zero_fold_tajD,total_1Mb$Thetas_tremula_3UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_zero_fold_tajD,total_1Mb$Thetas_tremula_5UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_zero_fold_tajD,total_1Mb$Thetas_tremula_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_zero_fold_tajD,total_1Mb$Thetas_tremula_intergenic_tajD,method="spearman")
###tremula-3UTR
cor.test(total_1Mb$Thetas_tremula_3UTR_tajD,total_1Mb$Thetas_tremula_5UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_3UTR_tajD,total_1Mb$Thetas_tremula_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_3UTR_tajD,total_1Mb$Thetas_tremula_intergenic_tajD,method="spearman")
###tremula-5UTR
cor.test(total_1Mb$Thetas_tremula_5UTR_tajD,total_1Mb$Thetas_tremula_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_5UTR_tajD,total_1Mb$Thetas_tremula_intergenic_tajD,method="spearman")
###tremula-intron
cor.test(total_1Mb$Thetas_tremula_intron_tajD,total_1Mb$Thetas_tremula_intergenic_tajD,method="spearman")


###tremuloides
##tremuloides-4fold
cor.test(total_1Mb$Thetas_tremuloides_four_fold_tajD,total_1Mb$Thetas_tremuloides_zero_fold_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_four_fold_tajD,total_1Mb$Thetas_tremuloides_3UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_four_fold_tajD,total_1Mb$Thetas_tremuloides_5UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_four_fold_tajD,total_1Mb$Thetas_tremuloides_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_four_fold_tajD,total_1Mb$Thetas_tremuloides_intergenic_tajD,method="spearman")
##tremuloides_0-fold
cor.test(total_1Mb$Thetas_tremuloides_zero_fold_tajD,total_1Mb$Thetas_tremuloides_3UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_zero_fold_tajD,total_1Mb$Thetas_tremuloides_5UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_zero_fold_tajD,total_1Mb$Thetas_tremuloides_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_zero_fold_tajD,total_1Mb$Thetas_tremuloides_intergenic_tajD,method="spearman")
###tremuloides-3UTR
cor.test(total_1Mb$Thetas_tremuloides_3UTR_tajD,total_1Mb$Thetas_tremuloides_5UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_3UTR_tajD,total_1Mb$Thetas_tremuloides_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_3UTR_tajD,total_1Mb$Thetas_tremuloides_intergenic_tajD,method="spearman")
###tremuloides-5UTR
cor.test(total_1Mb$Thetas_tremuloides_5UTR_tajD,total_1Mb$Thetas_tremuloides_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_5UTR_tajD,total_1Mb$Thetas_tremuloides_intergenic_tajD,method="spearman")
###tremuloides-intron
cor.test(total_1Mb$Thetas_tremuloides_intron_tajD,total_1Mb$Thetas_tremuloides_intergenic_tajD,method="spearman")


###trichocarpa
##trichocarpa-4fold
cor.test(total_1Mb$Thetas_trichocarpa_four_fold_tajD,total_1Mb$Thetas_trichocarpa_zero_fold_tajD,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_four_fold_tajD,total_1Mb$Thetas_trichocarpa_3UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_four_fold_tajD,total_1Mb$Thetas_trichocarpa_5UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_four_fold_tajD,total_1Mb$Thetas_trichocarpa_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_four_fold_tajD,total_1Mb$Thetas_trichocarpa_intergenic_tajD,method="spearman")
##trichocarpa_0-fold
cor.test(total_1Mb$Thetas_trichocarpa_zero_fold_tajD,total_1Mb$Thetas_trichocarpa_3UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_zero_fold_tajD,total_1Mb$Thetas_trichocarpa_5UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_zero_fold_tajD,total_1Mb$Thetas_trichocarpa_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_zero_fold_tajD,total_1Mb$Thetas_trichocarpa_intergenic_tajD,method="spearman")
###trichocarpa-3UTR
cor.test(total_1Mb$Thetas_trichocarpa_3UTR_tajD,total_1Mb$Thetas_trichocarpa_5UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_3UTR_tajD,total_1Mb$Thetas_trichocarpa_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_3UTR_tajD,total_1Mb$Thetas_trichocarpa_intergenic_tajD,method="spearman")
###trichocarpa-5UTR
cor.test(total_1Mb$Thetas_trichocarpa_5UTR_tajD,total_1Mb$Thetas_trichocarpa_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_trichocarpa_5UTR_tajD,total_1Mb$Thetas_trichocarpa_intergenic_tajD,method="spearman")
###trichocarpa-intron
cor.test(total_1Mb$Thetas_trichocarpa_intron_tajD,total_1Mb$Thetas_trichocarpa_intergenic_tajD,method="spearman")



####species-level comparison

###tp
#0-fold
cor.test(total_1Mb$Thetas_tremula_zero_fold_tP,total_1Mb$Thetas_tremuloides_zero_fold_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_zero_fold_tP,total_1Mb$Thetas_trichocarpa_zero_fold_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_zero_fold_tP,total_1Mb$Thetas_trichocarpa_zero_fold_tP,method="spearman")
#4-fold
cor.test(total_1Mb$Thetas_tremula_four_fold_tP,total_1Mb$Thetas_tremuloides_four_fold_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_four_fold_tP,total_1Mb$Thetas_trichocarpa_four_fold_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_four_fold_tP,total_1Mb$Thetas_trichocarpa_four_fold_tP,method="spearman")
#intron
cor.test(total_1Mb$Thetas_tremula_intron_tP,total_1Mb$Thetas_tremuloides_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_intron_tP,total_1Mb$Thetas_trichocarpa_intron_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_intron_tP,total_1Mb$Thetas_trichocarpa_intron_tP,method="spearman")
#3UTR
cor.test(total_1Mb$Thetas_tremula_3UTR_tP,total_1Mb$Thetas_tremuloides_3UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_3UTR_tP,total_1Mb$Thetas_trichocarpa_3UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_3UTR_tP,total_1Mb$Thetas_trichocarpa_3UTR_tP,method="spearman")
#5UTR
cor.test(total_1Mb$Thetas_tremula_5UTR_tP,total_1Mb$Thetas_tremuloides_5UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_5UTR_tP,total_1Mb$Thetas_trichocarpa_5UTR_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_5UTR_tP,total_1Mb$Thetas_trichocarpa_5UTR_tP,method="spearman")
#intergenic
cor.test(total_1Mb$Thetas_tremula_intergenic_tP,total_1Mb$Thetas_tremuloides_intergenic_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremula_intergenic_tP,total_1Mb$Thetas_trichocarpa_intergenic_tP,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_intergenic_tP,total_1Mb$Thetas_trichocarpa_intergenic_tP,method="spearman")




###tajD
#0-fold
cor.test(total_1Mb$Thetas_tremula_zero_fold_tajD,total_1Mb$Thetas_tremuloides_zero_fold_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_zero_fold_tajD,total_1Mb$Thetas_trichocarpa_zero_fold_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_zero_fold_tajD,total_1Mb$Thetas_trichocarpa_zero_fold_tajD,method="spearman")
#4-fold
cor.test(total_1Mb$Thetas_tremula_four_fold_tajD,total_1Mb$Thetas_tremuloides_four_fold_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_four_fold_tajD,total_1Mb$Thetas_trichocarpa_four_fold_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_four_fold_tajD,total_1Mb$Thetas_trichocarpa_four_fold_tajD,method="spearman")
#intron
cor.test(total_1Mb$Thetas_tremula_intron_tajD,total_1Mb$Thetas_tremuloides_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_intron_tajD,total_1Mb$Thetas_trichocarpa_intron_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_intron_tajD,total_1Mb$Thetas_trichocarpa_intron_tajD,method="spearman")
#3UTR
cor.test(total_1Mb$Thetas_tremula_3UTR_tajD,total_1Mb$Thetas_tremuloides_3UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_3UTR_tajD,total_1Mb$Thetas_trichocarpa_3UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_3UTR_tajD,total_1Mb$Thetas_trichocarpa_3UTR_tajD,method="spearman")
#5UTR
cor.test(total_1Mb$Thetas_tremula_5UTR_tajD,total_1Mb$Thetas_tremuloides_5UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_5UTR_tajD,total_1Mb$Thetas_trichocarpa_5UTR_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_5UTR_tajD,total_1Mb$Thetas_trichocarpa_5UTR_tajD,method="spearman")
#intergenic
cor.test(total_1Mb$Thetas_tremula_intergenic_tajD,total_1Mb$Thetas_tremuloides_intergenic_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremula_intergenic_tajD,total_1Mb$Thetas_trichocarpa_intergenic_tajD,method="spearman")
cor.test(total_1Mb$Thetas_tremuloides_intergenic_tajD,total_1Mb$Thetas_trichocarpa_intergenic_tajD,method="spearman")


###dxy
###tremula-tremuloides
##0-fold
cor.test(total_1Mb$tremula_tremuloides_zero_fold_dxy,total_1Mb$tremula_tremuloides_four_fold_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_zero_fold_dxy,total_1Mb$tremula_tremuloides_intron_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_zero_fold_dxy,total_1Mb$tremula_tremuloides_5UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_zero_fold_dxy,total_1Mb$tremula_tremuloides_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_zero_fold_dxy,total_1Mb$tremula_tremuloides_intergenic_dxy,method="spearman")
##4-fold
cor.test(total_1Mb$tremula_tremuloides_four_fold_dxy,total_1Mb$tremula_tremuloides_intron_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_four_fold_dxy,total_1Mb$tremula_tremuloides_5UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_four_fold_dxy,total_1Mb$tremula_tremuloides_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_four_fold_dxy,total_1Mb$tremula_tremuloides_intergenic_dxy,method="spearman")
###intron
cor.test(total_1Mb$tremula_tremuloides_intron_dxy,total_1Mb$tremula_tremuloides_5UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_intron_dxy,total_1Mb$tremula_tremuloides_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_intron_dxy,total_1Mb$tremula_tremuloides_intergenic_dxy,method="spearman")
##5UTR
cor.test(total_1Mb$tremula_tremuloides_5UTR_dxy,total_1Mb$tremula_tremuloides_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_5UTR_dxy,total_1Mb$tremula_tremuloides_intergenic_dxy,method="spearman")
##3UTR
cor.test(total_1Mb$tremula_tremuloides_3UTR_dxy,total_1Mb$tremula_tremuloides_intergenic_dxy,method="spearman")

###tremula-trichocarpa
##0-fold
cor.test(total_1Mb$tremula_trichocarpa_zero_fold_dxy,total_1Mb$tremula_trichocarpa_four_fold_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_zero_fold_dxy,total_1Mb$tremula_trichocarpa_intron_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_zero_fold_dxy,total_1Mb$tremula_trichocarpa_5UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_zero_fold_dxy,total_1Mb$tremula_trichocarpa_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_zero_fold_dxy,total_1Mb$tremula_trichocarpa_intergenic_dxy,method="spearman")
##4-fold
cor.test(total_1Mb$tremula_trichocarpa_four_fold_dxy,total_1Mb$tremula_trichocarpa_intron_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_four_fold_dxy,total_1Mb$tremula_trichocarpa_5UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_four_fold_dxy,total_1Mb$tremula_trichocarpa_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_four_fold_dxy,total_1Mb$tremula_trichocarpa_intergenic_dxy,method="spearman")
###intron
cor.test(total_1Mb$tremula_trichocarpa_intron_dxy,total_1Mb$tremula_trichocarpa_5UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_intron_dxy,total_1Mb$tremula_trichocarpa_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_intron_dxy,total_1Mb$tremula_trichocarpa_intergenic_dxy,method="spearman")
##5UTR
cor.test(total_1Mb$tremula_trichocarpa_5UTR_dxy,total_1Mb$tremula_trichocarpa_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_5UTR_dxy,total_1Mb$tremula_trichocarpa_intergenic_dxy,method="spearman")
##3UTR
cor.test(total_1Mb$tremula_trichocarpa_3UTR_dxy,total_1Mb$tremula_trichocarpa_intergenic_dxy,method="spearman")


###tremuloides_trichocarpa
##0-fold
cor.test(total_1Mb$tremuloides_trichocarpa_zero_fold_dxy,total_1Mb$tremuloides_trichocarpa_four_fold_dxy,method="spearman")
cor.test(total_1Mb$tremuloides_trichocarpa_zero_fold_dxy,total_1Mb$tremuloides_trichocarpa_intron_dxy,method="spearman")
cor.test(total_1Mb$tremuloides_trichocarpa_zero_fold_dxy,total_1Mb$tremuloides_trichocarpa_5UTR_dxy,method="spearman")
cor.test(total_1Mb$tremuloides_trichocarpa_zero_fold_dxy,total_1Mb$tremuloides_trichocarpa_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremuloides_trichocarpa_zero_fold_dxy,total_1Mb$tremuloides_trichocarpa_intergenic_dxy,method="spearman")
##4-fold
cor.test(total_1Mb$tremuloides_trichocarpa_four_fold_dxy,total_1Mb$tremuloides_trichocarpa_intron_dxy,method="spearman")
cor.test(total_1Mb$tremuloides_trichocarpa_four_fold_dxy,total_1Mb$tremuloides_trichocarpa_5UTR_dxy,method="spearman")
cor.test(total_1Mb$tremuloides_trichocarpa_four_fold_dxy,total_1Mb$tremuloides_trichocarpa_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremuloides_trichocarpa_four_fold_dxy,total_1Mb$tremuloides_trichocarpa_intergenic_dxy,method="spearman")
###intron
cor.test(total_1Mb$tremuloides_trichocarpa_intron_dxy,total_1Mb$tremuloides_trichocarpa_5UTR_dxy,method="spearman")
cor.test(total_1Mb$tremuloides_trichocarpa_intron_dxy,total_1Mb$tremuloides_trichocarpa_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremuloides_trichocarpa_intron_dxy,total_1Mb$tremuloides_trichocarpa_intergenic_dxy,method="spearman")
##5UTR
cor.test(total_1Mb$tremuloides_trichocarpa_5UTR_dxy,total_1Mb$tremuloides_trichocarpa_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremuloides_trichocarpa_5UTR_dxy,total_1Mb$tremuloides_trichocarpa_intergenic_dxy,method="spearman")
##3UTR
cor.test(total_1Mb$tremuloides_trichocarpa_3UTR_dxy,total_1Mb$tremuloides_trichocarpa_intergenic_dxy,method="spearman")


#####species-level comparison

##0-fold
cor.test(total_1Mb$tremula_tremuloides_zero_fold_dxy,total_1Mb$tremula_trichocarpa_zero_fold_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_zero_fold_dxy,total_1Mb$tremuloides_trichocarpa_zero_fold_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_zero_fold_dxy,total_1Mb$tremuloides_trichocarpa_zero_fold_dxy,method="spearman")

##4-fold
cor.test(total_1Mb$tremula_tremuloides_four_fold_dxy,total_1Mb$tremula_trichocarpa_four_fold_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_four_fold_dxy,total_1Mb$tremuloides_trichocarpa_four_fold_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_four_fold_dxy,total_1Mb$tremuloides_trichocarpa_four_fold_dxy,method="spearman")

##intron
cor.test(total_1Mb$tremula_tremuloides_intron_dxy,total_1Mb$tremula_trichocarpa_intron_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_intron_dxy,total_1Mb$tremuloides_trichocarpa_intron_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_intron_dxy,total_1Mb$tremuloides_trichocarpa_intron_dxy,method="spearman")

##3UTR
cor.test(total_1Mb$tremula_tremuloides_3UTR_dxy,total_1Mb$tremula_trichocarpa_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_3UTR_dxy,total_1Mb$tremuloides_trichocarpa_3UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_3UTR_dxy,total_1Mb$tremuloides_trichocarpa_3UTR_dxy,method="spearman")

##5UTR
cor.test(total_1Mb$tremula_tremuloides_5UTR_dxy,total_1Mb$tremula_trichocarpa_5UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_5UTR_dxy,total_1Mb$tremuloides_trichocarpa_5UTR_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_5UTR_dxy,total_1Mb$tremuloides_trichocarpa_5UTR_dxy,method="spearman")

##intergenic
cor.test(total_1Mb$tremula_tremuloides_intergenic_dxy,total_1Mb$tremula_trichocarpa_intergenic_dxy,method="spearman")
cor.test(total_1Mb$tremula_tremuloides_intergenic_dxy,total_1Mb$tremuloides_trichocarpa_intergenic_dxy,method="spearman")
cor.test(total_1Mb$tremula_trichocarpa_intergenic_dxy,total_1Mb$tremuloides_trichocarpa_intergenic_dxy,method="spearman")




























