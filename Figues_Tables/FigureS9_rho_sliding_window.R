#! /usr/bin/Rscript --no-save --no-restore

library(RColorBrewer)
colors <- brewer.pal(9,"Set1")[c(5,2,3)]

setwd("/proj/b2010014/GenomePaper/population_genetics/pan_genome/plots/thetas_tajD")

args=(commandArgs(TRUE))
window <- args[1]

wd=paste("/proj/b2010014/GenomePaper/population_genetics/pan_genome/summary/",window,sep="")

total_new=read.table(paste(wd,"/",window,".summary.txt",sep=""),header=T)

total_new$tremula_ldhat_new=total_new$Ldhat_tremula/1000
total_new$tremuloides_ldhat_new=total_new$Ldhat_tremuloides/1000
total_new$trichocarpa_ldhat_new=total_new$Ldhat_trichocarpa/1000


total_new[which(total_new$Ldhat_tremula_num<50),]$tremula_ldhat_new="NA"
total_new[which(total_new$Ldhat_tremuloides_num<50),]$tremuloides_ldhat_new="NA"
total_new[which(total_new$Ldhat_trichocarpa_num<50),]$trichocarpa_ldhat_new="NA"

chr=as.list(matrix(,19))
chr_pos=as.list(matrix(,19))
rho_tremula=as.list(matrix(,19))
rho_tremuloides=as.list(matrix(,19))
rho_trichocarpa=as.list(matrix(,19))

for (i in 1:19) { 
  if (i<10) {
    j=paste("Chr","0",i,sep="")
  }
  else {
    j=paste("Chr",i,sep="")
  }
  chr_pos[[i]]=total_new$Pos[total_new$Chr==j]
  rho_tremula[[i]]=total_new$tremula_ldhat_new[total_new$Chr==j]
  rho_tremuloides[[i]]=total_new$tremuloides_ldhat_new[total_new$Chr==j]
  rho_trichocarpa[[i]]=total_new$trichocarpa_ldhat_new[total_new$Chr==j]
}

tiff("3species.rho.sliding.window.tiff",width = 6, height = 4, units = 'in', res=350)
mat=matrix(c(rep(1,49),rep(2,25),rep(0,1),rep(3,21),rep(4,24),rep(5,26),rep(0,4),rep(6,28),rep(7,16),rep(8,20),rep(0,11),rep(9,13),rep(10,22),rep(11,14),rep(12,16),rep(0,10),rep(13,16),rep(14,13),rep(15,15),rep(16,15),rep(0,16),rep(17,16),rep(18,17),rep(19,16),rep(0,26)),byrow=T,nrow=6)
layout(mat)

par(mar=c(1.5,4,1,0.5))
plot(chr_pos[[1]]/1e6,rho_tremula[[1]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[1]]/1e6,rho_tremuloides[[1]],col=colors[2])
lines(chr_pos[[1]]/1e6,rho_trichocarpa[[1]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
par(mgp=c(0,1,0))
axis(side=2,at=seq(0,0.025,0.005),cex.axis=0.7)
mtext("Chr01 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[2]]/1e6,rho_tremula[[2]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[2]]/1e6,rho_tremuloides[[2]],col=colors[2])
lines(chr_pos[[2]]/1e6,rho_trichocarpa[[2]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
mtext("Chr02 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,4,1,0.5))
plot(chr_pos[[3]]/1e6,rho_tremula[[3]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[3]]/1e6,rho_tremuloides[[3]],col=colors[2])
lines(chr_pos[[3]]/1e6,rho_trichocarpa[[3]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
par(mgp=c(0,1,0))
axis(side=2,at=seq(0,0.025,0.005),cex.axis=0.7)
mtext("Chr03 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[4]]/1e6,rho_tremula[[4]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[4]]/1e6,rho_tremuloides[[4]],col=colors[2])
lines(chr_pos[[4]]/1e6,rho_trichocarpa[[4]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
mtext("Chr04 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
par(mgp=c(0,0.03,0))
plot(chr_pos[[5]]/1e6,rho_tremula[[5]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[5]]/1e6,rho_tremuloides[[5]],col=colors[2])
lines(chr_pos[[5]]/1e6,rho_trichocarpa[[5]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
mtext("Chr05 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,4,1,0.5))
plot(chr_pos[[6]]/1e6,rho_tremula[[6]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[6]]/1e6,rho_tremuloides[[6]],col=colors[2])
lines(chr_pos[[6]]/1e6,rho_trichocarpa[[6]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
par(mgp=c(0,1,0))
axis(side=2,at=seq(0,0.025,0.005),cex.axis=0.7)
mtext("Chr06 (Mb)",side=1,cex=0.5,line=0.2)
mtext(expression(rho),side=3,cex=0.8,line=-0.5,padj=6.3,adj=-0.22)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[7]]/1e6,rho_tremula[[7]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[7]]/1e6,rho_tremuloides[[7]],col=colors[2])
lines(chr_pos[[7]]/1e6,rho_trichocarpa[[7]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
mtext("Chr07 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[8]]/1e6,rho_tremula[[8]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[8]]/1e6,rho_tremuloides[[8]],col=colors[2])
lines(chr_pos[[8]]/1e6,rho_trichocarpa[[8]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
mtext("Chr08 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,4,1,0.5))
plot(chr_pos[[9]]/1e6,rho_tremula[[9]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[9]]/1e6,rho_tremuloides[[9]],col=colors[2])
lines(chr_pos[[9]]/1e6,rho_trichocarpa[[9]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
par(mgp=c(0,1,0))
axis(side=2,at=seq(0,0.025,0.005),cex.axis=0.7)
mtext("Chr09 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[10]]/1e6,rho_tremula[[10]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[10]]/1e6,rho_tremuloides[[10]],col=colors[2])
lines(chr_pos[[10]]/1e6,rho_trichocarpa[[10]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
mtext("Chr10 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[11]]/1e6,rho_tremula[[11]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[11]]/1e6,rho_tremuloides[[11]],col=colors[2])
lines(chr_pos[[11]]/1e6,rho_trichocarpa[[11]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
mtext("Chr11 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[12]]/1e6,rho_tremula[[12]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[12]]/1e6,rho_tremuloides[[12]],col=colors[2])
lines(chr_pos[[12]]/1e6,rho_trichocarpa[[12]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
mtext("Chr12 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,4,1,0.5))
plot(chr_pos[[13]]/1e6,rho_tremula[[13]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[13]]/1e6,rho_tremuloides[[13]],col=colors[2])
lines(chr_pos[[13]]/1e6,rho_trichocarpa[[13]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
par(mgp=c(0,1,0))
axis(side=2,at=seq(0,0.025,0.005),cex.axis=0.7)
mtext("Chr13 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[14]]/1e6,rho_tremula[[14]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[14]]/1e6,rho_tremuloides[[14]],col=colors[2])
lines(chr_pos[[14]]/1e6,rho_trichocarpa[[14]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
mtext("Chr14 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[15]]/1e6,rho_tremula[[15]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[15]]/1e6,rho_tremuloides[[15]],col=colors[2])
lines(chr_pos[[15]]/1e6,rho_trichocarpa[[15]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
mtext("Chr15 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[16]]/1e6,rho_tremula[[16]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[16]]/1e6,rho_tremuloides[[16]],col=colors[2])
lines(chr_pos[[16]]/1e6,rho_trichocarpa[[16]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
mtext("Chr16 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,4,1,0.5))
plot(chr_pos[[17]]/1e6,rho_tremula[[17]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[17]]/1e6,rho_tremuloides[[17]],col=colors[2])
lines(chr_pos[[17]]/1e6,rho_trichocarpa[[17]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
par(mgp=c(0,1,0))
axis(side=2,at=seq(0,0.025,0.005),cex.axis=0.7)
mtext("Chr17 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[18]]/1e6,rho_tremula[[18]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[18]]/1e6,rho_tremuloides[[18]],col=colors[2])
lines(chr_pos[[18]]/1e6,rho_trichocarpa[[18]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
mtext("Chr18 (Mb)",side=1,cex=0.5,line=0.2)

par(mar=c(1.5,1,1,0.5))
plot(chr_pos[[19]]/1e6,rho_tremula[[19]],axes=F,col=colors[1],xlab="",ylab="",type="l",ylim=c(0,0.025))
lines(chr_pos[[19]]/1e6,rho_tremuloides[[19]],col=colors[2])
lines(chr_pos[[19]]/1e6,rho_trichocarpa[[19]],col=colors[3])
par(mgp=c(0,-0.015,0))
axis(side=1,tck=-0.015,pos=0,cex.axis=0.5)
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



