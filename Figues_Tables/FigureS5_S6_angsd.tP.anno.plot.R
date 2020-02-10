#! /usr/bin/Rscript --no-save --no-restore

library(ggplot2)
library(RColorBrewer)

setwd("/proj/b2010014/GenomePaper/population_genetics/pan_genome/plots/anno")

args=(commandArgs(TRUE))
window <- args[1]

colors <- brewer.pal(9,"Set1")[c(5,2,3)]


tremula_wd=paste("/proj/b2010014/GenomePaper/population_genetics/pan_genome/summary/",window,"/anno/tremula/",sep="")
tremuloides_wd=paste("/proj/b2010014/GenomePaper/population_genetics/pan_genome/summary/",window,"/anno/tremuloides/",sep="")
trichocarpa_wd=paste("/proj/b2010014/GenomePaper/population_genetics/pan_genome/summary/",window,"/anno/trichocarpa/",sep="")

####tremula
tremula_3UTR=read.table(paste(tremula_wd,"tremula.thetas.summary.3UTR.",window,".txt",sep=""),header=T)
tremula_5UTR=read.table(paste(tremula_wd,"tremula.thetas.summary.5UTR.",window,".txt",sep=""),header=T)
tremula_Intergenic=read.table(paste(tremula_wd,"tremula.thetas.summary.Intergenic.",window,".txt",sep=""),header=T)
tremula_Intron=read.table(paste(tremula_wd,"tremula.thetas.summary.Intron.",window,".txt",sep=""),header=T)
tremula_four_fold=read.table(paste(tremula_wd,"tremula.thetas.summary.four_fold.",window,".txt",sep=""),header=T)
tremula_zero_fold=read.table(paste(tremula_wd,"tremula.thetas.summary.zero_fold.",window,".txt",sep=""),header=T)
tremula_3UTR$anno="3'UTR"
tremula_3UTR$species="P.tremula"
tremula_5UTR$anno="5'UTR"
tremula_5UTR$species="P.tremula"
tremula_Intergenic$anno="Intergene"
tremula_Intergenic$species="P.tremula"
tremula_Intron$anno="Intron"
tremula_Intron$species="P.tremula"
tremula_four_fold$anno="4-fold"
tremula_four_fold$species="P.tremula"
tremula_zero_fold$anno="0-fold"
tremula_zero_fold$species="P.tremula"

#tremuloides
tremuloides_3UTR=read.table(paste(tremuloides_wd,"tremuloides.thetas.summary.3UTR.",window,".txt",sep=""),header=T)
tremuloides_5UTR=read.table(paste(tremuloides_wd,"tremuloides.thetas.summary.5UTR.",window,".txt",sep=""),header=T)
tremuloides_Intergenic=read.table(paste(tremuloides_wd,"tremuloides.thetas.summary.Intergenic.",window,".txt",sep=""),header=T)
tremuloides_Intron=read.table(paste(tremuloides_wd,"tremuloides.thetas.summary.Intron.",window,".txt",sep=""),header=T)
tremuloides_four_fold=read.table(paste(tremuloides_wd,"tremuloides.thetas.summary.four_fold.",window,".txt",sep=""),header=T)
tremuloides_zero_fold=read.table(paste(tremuloides_wd,"tremuloides.thetas.summary.zero_fold.",window,".txt",sep=""),header=T)

tremuloides_3UTR$anno="3'UTR"
tremuloides_3UTR$species="P.tremuloides"
tremuloides_5UTR$anno="5'UTR"
tremuloides_5UTR$species="P.tremuloides"
tremuloides_Intergenic$anno="Intergene"
tremuloides_Intergenic$species="P.tremuloides"
tremuloides_Intron$anno="Intron"
tremuloides_Intron$species="P.tremuloides"
tremuloides_four_fold$anno="4-fold"
tremuloides_four_fold$species="P.tremuloides"
tremuloides_zero_fold$anno="0-fold"
tremuloides_zero_fold$species="P.tremuloides"

#trichocarpa
trichocarpa_3UTR=read.table(paste(trichocarpa_wd,"trichocarpa.thetas.summary.3UTR.",window,".txt",sep=""),header=T)
trichocarpa_5UTR=read.table(paste(trichocarpa_wd,"trichocarpa.thetas.summary.5UTR.",window,".txt",sep=""),header=T)
trichocarpa_Intergenic=read.table(paste(trichocarpa_wd,"trichocarpa.thetas.summary.Intergenic.",window,".txt",sep=""),header=T)
trichocarpa_Intron=read.table(paste(trichocarpa_wd,"trichocarpa.thetas.summary.Intron.",window,".txt",sep=""),header=T)
trichocarpa_four_fold=read.table(paste(trichocarpa_wd,"trichocarpa.thetas.summary.four_fold.",window,".txt",sep=""),header=T)
trichocarpa_zero_fold=read.table(paste(trichocarpa_wd,"trichocarpa.thetas.summary.zero_fold.",window,".txt",sep=""),header=T)

trichocarpa_3UTR$anno="3'UTR"
trichocarpa_3UTR$species="P.trichocarpa"
trichocarpa_5UTR$anno="5'UTR"
trichocarpa_5UTR$species="P.trichocarpa"
trichocarpa_Intergenic$anno="Intergene"
trichocarpa_Intergenic$species="P.trichocarpa"
trichocarpa_Intron$anno="Intron"
trichocarpa_Intron$species="P.trichocarpa"
trichocarpa_four_fold$anno="4-fold"
trichocarpa_four_fold$species="P.trichocarpa"
trichocarpa_zero_fold$anno="0-fold"
trichocarpa_zero_fold$species="P.trichocarpa"

#total 

total=rbind(tremula_3UTR,tremula_5UTR,tremula_Intergenic,tremula_Intron,tremula_four_fold,tremula_zero_fold,tremuloides_3UTR,tremuloides_5UTR,tremuloides_Intergenic,tremuloides_Intron,tremuloides_four_fold,tremuloides_zero_fold,trichocarpa_3UTR,trichocarpa_5UTR,trichocarpa_Intergenic,trichocarpa_Intron,trichocarpa_four_fold,trichocarpa_zero_fold)

total$anno=factor(total$anno,levels=c("0-fold","4-fold","3'UTR","5'UTR","Intron","Intergene"))
total$species=factor(total$species,levels=c("P.tremula","P.tremuloides","P.trichocarpa"))

#save the RData
#save(total,file="total..RData")

#tP
tP=ggplot(total,aes(x=species,y=tP.norm,fill=species))+geom_boxplot()+scale_fill_manual(values=c("#FF7F00","#377EB8","#4DAF4A"))+ylim(0,0.04)
tP+facet_grid(.~anno)+xlab("Species")+ylab(expression(theta[pi]))+theme(axis.text.x=element_blank(),axis.ticks=element_blank())

ggsave(paste("tP.3species.anno.",window,".png",sep=""),width=8,height=6,dpi=300)

#tW
tW=ggplot(total,aes(x=species,y=tW.norm,fill=species))+geom_boxplot()+scale_fill_manual(values=c("#FF7F00","#377EB8","#4DAF4A"))+ylim(0,0.04)
tW+facet_grid(.~anno)+xlab("Species")+ylab(expression(theta[W]))+theme(axis.text.x=element_blank(),axis.ticks=element_blank())

ggsave(paste("tW.3species.anno.",window,".png",sep=""),width=8,height=6,dpi=300)

#tajD
tajD=ggplot(total,aes(x=species,y=tajD,fill=species))+geom_boxplot()+scale_fill_manual(values=c("#FF7F00","#377EB8","#4DAF4A"))
tajD+facet_grid(.~anno)+xlab("Species")+ylab("Tajima's D")+theme(axis.text.x=element_blank(),axis.ticks=element_blank())

ggsave(paste("tajD.3species.anno.",window,".png",sep=""),width=8,height=6,dpi=300)








