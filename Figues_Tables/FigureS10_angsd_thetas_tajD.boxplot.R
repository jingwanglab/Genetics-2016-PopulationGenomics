#! /usr/bin/Rscript --no-save --no-restore

library(RColorBrewer)
library(ggplot2)
library(gridExtra)

colors <- brewer.pal(9,"Set1")[c(5,2,3)]

setwd("/proj/b2010014/GenomePaper/population_genetics/pan_genome/plots/thetas_tajD")


args=(commandArgs(TRUE))
window <- args[1]

wd=paste("/proj/b2010014/GenomePaper/population_genetics/pan_genome/summary/",window,sep="")

tremula=read.table(paste(wd,"/tremula.thetas.summary.",window,".txt",sep=""),header=T)
tremuloides=read.table(paste(wd,"/tremuloides.thetas.summary.",window,".txt",sep=""),header=T)
trichocarpa=read.table(paste(wd,"/trichocarpa.thetas.summary.",window,".txt",sep=""),header=T)

tremula_new=as.data.frame(cbind(tremula$tP.norm,tremula$tW.norm,tremula$tajD))
names(tremula_new)=c("tP","tW","tajD")
tremula_new$species="P.tremula"

tremuloides_new=as.data.frame(cbind(tremuloides$tP.norm,tremuloides$tW.norm,tremuloides$tajD))
names(tremuloides_new)=c("tP","tW","tajD")
tremuloides_new$species="P.tremuloides"

trichocarpa_new=as.data.frame(cbind(trichocarpa$tP.norm,trichocarpa$tW.norm,trichocarpa$tajD))
names(trichocarpa_new)=c("tP","tW","tajD")
trichocarpa_new$species="P.trichocarpa"

total=rbind(tremula_new,tremuloides_new,trichocarpa_new)

### for recombination rate rho
tremula_rho=read.table(paste(wd,"/ldhat.tremula.",window,".summary.txt",sep=""),header=T)
tremuloides_rho=read.table(paste(wd,"/ldhat.tremuloides.",window,".summary.txt",sep=""),header=T)
trichocarpa_rho=read.table(paste(wd,"/ldhat.trichocarpa.",window,".summary.txt",sep=""),header=T)

tremula_rho_new=tremula_rho[which(tremula_rho$Num>=10),]
tremula_rho_new$ldhat_new=tremula_rho_new$ldhat/1000
tremula_rho_new$species="P.tremula"

tremuloides_rho_new=tremuloides_rho[which(tremuloides_rho$Num>=10),]
tremuloides_rho_new$ldhat_new=tremuloides_rho_new$ldhat/1000
tremuloides_rho_new$species="P.tremuloides"

trichocarpa_rho_new=trichocarpa_rho[which(trichocarpa_rho$Num>=10),]
trichocarpa_rho_new$ldhat_new=trichocarpa_rho_new$ldhat/1000
trichocarpa_rho_new$species="P.trichocarpa"

total_rho=rbind(tremula_rho_new,tremuloides_rho_new,trichocarpa_rho_new)

p=ggplot(total,aes(x=species,y=tP,fill=species))+geom_boxplot(fill=c(colors[1],colors[2],colors[3]))+xlab("Species")+ylab(expression(theta[pi]))+ggtitle("(a)")+theme(plot.title=element_text(hjust=0))
w=ggplot(total,aes(x=species,y=tW,fill=species))+geom_boxplot(fill=c(colors[1],colors[2],colors[3]))+xlab("Species")+ylab(expression(theta[W]))+ggtitle("(b)")+theme(plot.title=element_text(hjust=0))
d=ggplot(total,aes(x=species,y=tajD,fill=species))+geom_boxplot(fill=c(colors[1],colors[2],colors[3]))+xlab("Species")+ylab("Tajima's D")+ggtitle("(c)")+theme(plot.title=element_text(hjust=0))
rho=ggplot(total_rho,aes(x=species,y=ldhat_new,fill=species))+geom_boxplot(fill=c(colors[1],colors[2],colors[3]))+xlab("Species")+ylab(expression(rho))+ggtitle("(d)")+theme(plot.title=element_text(hjust=0))
grid.arrange(p,w,d,rho,ncol=2)
g=arrangeGrob(p,w,d,rho,ncol=2)

#ggsave(paste("3species.thetas_tajD_rho.",window,".png",sep=""),width=7,height=5.5,dpi=300,g)
ggsave(paste("3species.tP_W_tajD_rho.",window,".tiff",sep=""),width=4,height=4,dpi=300,g)

