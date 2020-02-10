#! /usr/bin/Rscript --no-save --no-restore

library(RColorBrewer)
colors <- brewer.pal(9,"Set1")[c(5,2,3)]

setwd("/proj/b2010014/GenomePaper/population_genetics/pan_genome/plots/thetas_tajD")

args=(commandArgs(TRUE))
window <- args[1]

wd=paste("/proj/b2010014/GenomePaper/population_genetics/pan_genome/summary/",window,sep="")

tremula=read.table(paste(wd,"/tremula.thetas.summary.",window,".txt",sep=""),header=T)
tremuloides=read.table(paste(wd,"/tremuloides.thetas.summary.",window,".txt",sep=""),header=T)
trichocarpa=read.table(paste(wd,"/trichocarpa.thetas.summary.",window,".txt",sep=""),header=T)

chr=as.list(matrix(,19))
chr_pos=as.list(matrix(,19))
tajD_tremula=as.list(matrix(,19))
tajD_tremuloides=as.list(matrix(,19))
tajD_trichocarpa=as.list(matrix(,19))

for (i in 1:19) { 
  if (i<10) {
    j=paste("Chr","0",i,sep="")
  }
  else {
    j=paste("Chr",i,sep="")
  }
  chr_pos[[i]]=tremula$Pos[tremula$Chr==j]
  tajD_tremula[[i]]=tremula$tajD[tremula$Chr==j]
  tajD_tremuloides[[i]]=tremuloides$tajD[tremuloides$Chr==j]
  tajD_trichocarpa[[i]]=trichocarpa$tajD[trichocarpa$Chr==j]
}

tiff("3species.tajD.sliding.window.tiff",width = 6, height = 4, units = 'in', res=350)
mat=matrix(c(rep(1,49),rep(2,25),rep(0,1),rep(3,21),rep(4,24),rep(5,26),rep(0,4),rep(6,28),rep(7,16),rep(8,20),rep(0,11),rep(9,13),rep(10,22),rep(11,14),rep(12,16),rep(0,10),rep(13,16),rep(14,13),rep(15,15),rep(16,15),rep(0,16),rep(17,16),rep(18,17),rep(19,16),rep(0,26)),byrow=T,nrow=6)
layout(mat)

par(mar=c(1.5,4,1,0.5))
plot(chr_pos[[1]]/1e6,tajD_tremula[[1]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[1]]/1e6,tajD_tremuloides[[1]],col=colors[2])
lines(chr_pos[[1]]/1e6,tajD_trichocarpa[[1]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
par(mgp=c(0,1,0))
axis(side=2,at=seq(-2,1.5,0.5),cex.axis=0.7)
mtext("Chr01 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[2]]/1e6,tajD_tremula[[2]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[2]]/1e6,tajD_tremuloides[[2]],col=colors[2])
lines(chr_pos[[2]]/1e6,tajD_trichocarpa[[2]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr02 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,4,1,0.5))
plot(chr_pos[[3]]/1e6,tajD_tremula[[3]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[3]]/1e6,tajD_tremuloides[[3]],col=colors[2])
lines(chr_pos[[3]]/1e6,tajD_trichocarpa[[3]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
par(mgp=c(0,1,0))
axis(side=2,at=seq(-2,1.5,0.5),cex.axis=0.7)
mtext("Chr03 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[4]]/1e6,tajD_tremula[[4]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[4]]/1e6,tajD_tremuloides[[4]],col=colors[2])
lines(chr_pos[[4]]/1e6,tajD_trichocarpa[[4]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr04 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
par(mgp=c(0,0.03,0))
plot(chr_pos[[5]]/1e6,tajD_tremula[[5]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[5]]/1e6,tajD_tremuloides[[5]],col=colors[2])
lines(chr_pos[[5]]/1e6,tajD_trichocarpa[[5]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr05 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,4,1,0.5))
plot(chr_pos[[6]]/1e6,tajD_tremula[[6]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[6]]/1e6,tajD_tremuloides[[6]],col=colors[2])
lines(chr_pos[[6]]/1e6,tajD_trichocarpa[[6]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
par(mgp=c(0,1,0))
axis(side=2,at=seq(-2,1.5,0.5),cex.axis=0.7)
mtext("Chr06 (Mb)",side=1,cex=0.5,line=0.2)
#mtext("Tajima's D",side=2,cex=0.8,line=-0.5,padj=6.3,adj=-0.22)
mtext("Tajima's D",side=2,line=2,cex=0.7,adj=0.5)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[7]]/1e6,tajD_tremula[[7]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[7]]/1e6,tajD_tremuloides[[7]],col=colors[2])
lines(chr_pos[[7]]/1e6,tajD_trichocarpa[[7]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr07 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[8]]/1e6,tajD_tremula[[8]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[8]]/1e6,tajD_tremuloides[[8]],col=colors[2])
lines(chr_pos[[8]]/1e6,tajD_trichocarpa[[8]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr08 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,4,1,0.5))
plot(chr_pos[[9]]/1e6,tajD_tremula[[9]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[9]]/1e6,tajD_tremuloides[[9]],col=colors[2])
lines(chr_pos[[9]]/1e6,tajD_trichocarpa[[9]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
par(mgp=c(0,1,0))
axis(side=2,at=seq(-2,1.5,0.5),cex.axis=0.7)
mtext("Chr09 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[10]]/1e6,tajD_tremula[[10]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[10]]/1e6,tajD_tremuloides[[10]],col=colors[2])
lines(chr_pos[[10]]/1e6,tajD_trichocarpa[[10]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr10 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[11]]/1e6,tajD_tremula[[11]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[11]]/1e6,tajD_tremuloides[[11]],col=colors[2])
lines(chr_pos[[11]]/1e6,tajD_trichocarpa[[11]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr11 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[12]]/1e6,tajD_tremula[[12]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[12]]/1e6,tajD_tremuloides[[12]],col=colors[2])
lines(chr_pos[[12]]/1e6,tajD_trichocarpa[[12]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr12 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,4,1,0.5))
plot(chr_pos[[13]]/1e6,tajD_tremula[[13]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[13]]/1e6,tajD_tremuloides[[13]],col=colors[2])
lines(chr_pos[[13]]/1e6,tajD_trichocarpa[[13]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
par(mgp=c(0,1,0))
axis(side=2,at=seq(-2,1.5,0.5),cex.axis=0.7)
mtext("Chr13 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[14]]/1e6,tajD_tremula[[14]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[14]]/1e6,tajD_tremuloides[[14]],col=colors[2])
lines(chr_pos[[14]]/1e6,tajD_trichocarpa[[14]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr14 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[15]]/1e6,tajD_tremula[[15]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[15]]/1e6,tajD_tremuloides[[15]],col=colors[2])
lines(chr_pos[[15]]/1e6,tajD_trichocarpa[[15]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr15 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[16]]/1e6,tajD_tremula[[16]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[16]]/1e6,tajD_tremuloides[[16]],col=colors[2])
lines(chr_pos[[16]]/1e6,tajD_trichocarpa[[16]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr16 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,4,1,0.5))
plot(chr_pos[[17]]/1e6,tajD_tremula[[17]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[17]]/1e6,tajD_tremuloides[[17]],col=colors[2])
lines(chr_pos[[17]]/1e6,tajD_trichocarpa[[17]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
par(mgp=c(0,1,0))
axis(side=2,at=seq(-2,1.5,0.5),cex.axis=0.7)
mtext("Chr17 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[18]]/1e6,tajD_tremula[[18]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[18]]/1e6,tajD_tremuloides[[18]],col=colors[2])
lines(chr_pos[[18]]/1e6,tajD_trichocarpa[[18]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr18 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[19]]/1e6,tajD_tremula[[19]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(-2,1.5))
lines(chr_pos[[19]]/1e6,tajD_tremuloides[[19]],col=colors[2])
lines(chr_pos[[19]]/1e6,tajD_trichocarpa[[19]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=-2,cex.axis=0.5)
mtext("Chr19 (Mb)",side=1,cex=0.5,line=0.2)

 add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }

add_legend("top", legend=c("P.tremula", "P.tremuloides","P.trichocarpa"), lwd=1,col=c(colors[1],colors[2],colors[3]),horiz=TRUE, bty='n', cex=0.8)
dev.off()



