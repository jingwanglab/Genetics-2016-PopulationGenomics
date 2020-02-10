#! /usr/bin/Rscript --no-save --no-restore

setwd("/proj/b2010014/GenomePaper/population_genetics/pan_genome/plots/dfe")
library(gplots)
library(RColorBrewer)
library(lattice)

colors <- brewer.pal(9,"Set1")[c(5,2,3)]

###read.table
tremula_zero_fold=read.table("/proj/b2010014/GenomePaper/population_genetics/pan_genome/dfe/tremula/zero_fold/summary/tremula.dfe.zero_fold.summary",header=T)
tremuloides_zero_fold=read.table("/proj/b2010014/GenomePaper/population_genetics/pan_genome/dfe/tremuloides/zero_fold/summary/tremuloides.dfe.zero_fold.summary",header=T)
trichocarpa_zero_fold=read.table("/proj/b2010014/GenomePaper/population_genetics/pan_genome/dfe/trichocarpa/zero_fold/summary/trichocarpa.dfe.zero_fold.summary",header=T)


###real data
Nes=as.data.frame(as.matrix(rbind(tremula_zero_fold[1,c(9,10,11,12)],tremuloides_zero_fold[1,c(9,10,11,12)],trichocarpa_zero_fold[1,c(9,10,11,12)])))
Nes$Nes_large_10=rowSums(Nes[,c("Nes_10_100","Nes_100")])
Nes_new=as.table(as.matrix(Nes[,c(1,2,5)]))

Nes_025=as.data.frame(rbind(apply(tremula_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025)),apply(tremuloides_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025)),apply(trichocarpa_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.025))))
Nes_025$Nes_large_10=rowSums(Nes_025[,c("Nes_10_100","Nes_100")])
Nes_025_new=as.table(as.matrix(Nes_025[,c(1,2,5)]))


Nes_975=as.data.frame(rbind(apply(tremula_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975)),apply(tremuloides_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975)),apply(trichocarpa_zero_fold[c(2:201),c(9,10,11,12)],2,quantile,probs=c(.975))))
Nes_975$Nes_large_10=rowSums(Nes_975[,c("Nes_10_100","Nes_100")])
Nes_975_new=as.table(as.matrix(Nes_975[,c(1,2,5)]))


alpha=as.table(as.matrix(rbind(tremula_zero_fold[1,15],tremuloides_zero_fold[1,15],trichocarpa_zero_fold[1,15])))
alpha_025=as.table(rbind(quantile(tremula_zero_fold[c(2:201),15],probs=c(.025)),quantile(tremuloides_zero_fold[c(2:201),15],probs=c(.025)),quantile(trichocarpa_zero_fold[c(2:201),15],probs=c(.025))))
alpha_975=as.table(rbind(quantile(tremula_zero_fold[c(2:201),15],probs=c(.975)),quantile(tremuloides_zero_fold[c(2:201),15],probs=c(.975)),quantile(trichocarpa_zero_fold[c(2:201),15],probs=c(.975))))

omega=as.table(as.matrix(rbind(tremula_zero_fold[1,16],tremuloides_zero_fold[1,16],trichocarpa_zero_fold[1,16])))
omega_025=as.table(rbind(quantile(tremula_zero_fold[c(2:201),16],probs=c(.025)),quantile(tremuloides_zero_fold[c(2:201),16],probs=c(.025)),quantile(trichocarpa_zero_fold[c(2:201),16],probs=c(.025))))
omega_975=as.table(rbind(quantile(tremula_zero_fold[c(2:201),16],probs=c(.975)),quantile(tremuloides_zero_fold[c(2:201),16],probs=c(.975)),quantile(trichocarpa_zero_fold[c(2:201),16],probs=c(.975))))


png(filename="tremula_tremuloides.Nes.alpha.png",width=6.5,height=4,units='in',res=300)
layout(matrix(c(1,1,1,2,3), 1, 5, byrow = TRUE))
par(mar=c(3.5,4.5,2,1))
###Nes
Nes_plot=barplot2(Nes_new,beside=T,names.arg=c(expression(paste("0<",N[e],"s<1")),expression(paste("1<",N[e],"s<10")),expression(paste(N[e],"s>10"))),plot.ci=T,ci.u=Nes_975_new,ci.l=Nes_025_new,col=c(colors[1],colors[2],colors[3]),ylab="Fraction of sites",ylim=c(0,0.8))
legend("top",c(expression(italic(P.tremula)),expression(italic(P.tremuloides)),expression(italic(P.trichocarpa))),bty="n",fill=c(colors[1],colors[2],colors[3]),horiz=T)
mtext("(a)",side=3,line=0.05,adj=0.5,font=1.5,cex=1)

###alpha
par(mar=c(3.5,2.5,2,1))
alpha_plot=barplot2(alpha,beside=T,plot.ci=T,ci.u=alpha_975,ci.l=alpha_025,col=c(colors[1],colors[2],colors[3]),names.arg=expression(alpha),ylim=c(0,0.7))
mtext("(b)",side=3,line=0.05,adj=0.5,font=1.5,cex=1)

par(mar=c(3.5,2.5,2,1))
omega_plot=barplot2(omega,beside=T,plot.ci=T,ci.u=omega_975,ci.l=omega_025,col=c(colors[1],colors[2],colors[3]),names.arg=expression(omega),ylim=c(0,0.3))
mtext("(c)",side=3,line=0.05,adj=0.5,font=1.5,cex=1)
dev.off()


