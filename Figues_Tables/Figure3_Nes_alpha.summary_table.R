#! /usr/bin/Rscript --no-save --no-restore

setwd("/proj/b2011141/nobackup/genomic_selection_paper/DFE/plot")
library(gplots)
library(RColorBrewer)

colors <- brewer.pal(9,"Set1")[c(1:9)]

###read.table
tremula_zero_fold=read.table("/proj/b2010014/GenomePaper/population_genetics/pan_genome/dfe/tremula/zero_fold/summary/tremula.dfe.zero_fold.summary",header=T)
tremuloides_zero_fold=read.table("/proj/b2010014/GenomePaper/population_genetics/pan_genome/dfe/tremuloides/zero_fold/summary/tremuloides.dfe.zero_fold.summary",header=T)
trichocarpa_zero_fold=read.table("/proj/b2010014/GenomePaper/population_genetics/pan_genome/dfe/trichocarpa/zero_fold/summary/trichocarpa.dfe.zero_fold.summary",header=T)

tremuloides_zero_fold$N2_N1=tremuloides_zero_fold$N2/100
tremuloides_zero_fold$t_N2=tremuloides_zero_fold$t2/tremuloides_zero_fold$N2

tremula_zero_fold$N2_N1=tremula_zero_fold$N2/100
tremula_zero_fold$t_N2=tremula_zero_fold$t2/tremula_zero_fold$N2

trichocarpa_zero_fold$N2_N1=trichocarpa_zero_fold$N2/100
trichocarpa_zero_fold$t_N2=trichocarpa_zero_fold$t2/trichocarpa_zero_fold$N2

###real data
tremula=as.table(as.matrix(rbind(tremula_zero_fold[1,c(17,18,5,9,10,11,12,15,16)])))

tremula_025=as.table(rbind(apply(tremula_zero_fold[c(2:201),c(17,18,5,9,10,11,12,15,16)],2,quantile,probs=c(.025))))

tremula_975=as.table(rbind(apply(tremula_zero_fold[c(2:201),c(17,18,5,9,10,11,12,15,16)],2,quantile,probs=c(.975))))

rbind(tremula,tremula_025,tremula_975)

#tremuloides
tremuloides=as.table(as.matrix(rbind(tremuloides_zero_fold[1,c(17,18,5,9,10,11,12,15,16)])))

tremuloides_025=as.table(rbind(apply(tremuloides_zero_fold[c(2:201),c(17,18,5,9,10,11,12,15,16)],2,quantile,probs=c(.025))))

tremuloides_975=as.table(rbind(apply(tremuloides_zero_fold[c(2:201),c(17,18,5,9,10,11,12,15,16)],2,quantile,probs=c(.975))))

rbind(tremuloides,tremuloides_025,tremuloides_975)

#trichocarpa
trichocarpa=as.table(as.matrix(rbind(trichocarpa_zero_fold[1,c(17,18,5,9,10,11,12,15,16)])))

trichocarpa_025=as.table(rbind(apply(trichocarpa_zero_fold[c(2:201),c(17,18,5,9,10,11,12,15,16)],2,quantile,probs=c(.025))))

trichocarpa_975=as.table(rbind(apply(trichocarpa_zero_fold[c(2:201),c(17,18,5,9,10,11,12,15,16)],2,quantile,probs=c(.975))))

rbind(trichocarpa,trichocarpa_025,trichocarpa_975)



