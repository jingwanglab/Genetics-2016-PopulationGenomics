library(RColorBrewer)
colors <- brewer.pal(9,"Set1")[c(5,2,3)]

setwd("/Users/Jing/Dropbox/paperI/data/1Mb")

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

###plots of correlation between gene number and recombination
tiff(filename="3species.1Mb.gene_density.rho.tiff",width = 7, height = 4.5, units = 'in', res=300)
par(mfrow=c(1,3))
par(mar=c(4,4.5,2,1))

plot(total_1Mb$Gene_num,total_1Mb$tremula_ldhat_new,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab=expression(paste("P.tremula","(",rho,")",sep=" ")),xlab="Gene number")
z=lm(total_1Mb$tremula_ldhat_new~total_1Mb$Gene_num)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$Gene_num,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(50,0.005,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(50,0.0047,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(50,0.0047,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

plot(total_1Mb$Gene_num,total_1Mb$tremuloides_ldhat_new,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab=expression(paste("P.tremuloides","(",rho,")",sep=" ")),xlab="Gene number")
z=lm(total_1Mb$tremuloides_ldhat_new~total_1Mb$Gene_num)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$Gene_num,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(50,0.0133,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(50,0.0125,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(50,0.0125,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

plot(total_1Mb$Gene_num,total_1Mb$trichocarpa_ldhat_new,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab=expression(paste("P.trichocarpa","(",rho,")",sep=" ")),xlab="Gene number")
z=lm(total_1Mb$trichocarpa_ldhat_new~total_1Mb$Gene_num)
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$Gene_num,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(50,0.0033,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(50,0.0031,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(50,0.0031,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(c)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

dev.off()



###plots
png(filename="3species.1Mb.rho_pi_df.png",width = 5, height = 5, units = 'in', res=300)
par(mfrow=c(3,2))
par(mar=c(4,4,1,1))
par(tcl=-0.3, mgp=c(2, 0.6, 0))

#pi
plot(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[4-fold]))
z=lm(total_1Mb$Thetas_tremula_four_fold_tP~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$Thetas_tremula_four_fold_tP,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.001,0.017,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.001,0.015,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.002,0.015,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)
legend("bottomright",lty=NA,lwd=2,col=colors[1],pch=NA, bty="n", legend="P. tremula", text.font=3)


#fixed
plot(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[4-fold]))
z=lm(total_1Mb$tremula_trichocarpa_four_fold_fixed~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$tremula_trichocarpa_four_fold_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.0056,0.0275,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0056,0.025,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0056,0.025,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))


###tremuloides
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$Thetas_tremuloides_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[4-fold]))
z=lm(total_1Mb$Thetas_tremuloides_four_fold_tP~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$Thetas_tremuloides_four_fold_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.0025,0.02,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0025,0.018,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.018,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)
legend("bottomright",lty=NA,lwd=2,col=colors[2],pch=NA, bty="n", legend="P. tremuloides", text.font=3)

#fixed
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremuloides_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[4-fold]))
z=lm(total_1Mb$tremuloides_trichocarpa_four_fold_fixed~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$tremuloides_trichocarpa_four_fold_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.0148,0.024,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0025,0.022,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0148,0.022,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

#P.trichocarpa
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$Thetas_trichocarpa_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[4-fold]))
z=lm(total_1Mb$Thetas_trichocarpa_four_fold_tP~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)))
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$Thetas_trichocarpa_four_fold_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.00085,0.0092,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.00085,0.008,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0014,0.008,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(c)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)
legend("bottomright",lty=NA,lwd=2,col=colors[3],pch=NA, bty="n", legend="P. trichocarpa", text.font=3)

#fixed
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$tremula_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[4-fold]))
z=lm(total_1Mb$tremula_trichocarpa_four_fold_fixed~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),use="complete.obs")
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$tremula_trichocarpa_four_fold_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0036,0.0265,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.00085,0.0243,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0036,0.0243,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

dev.off()


###plots
png(filename="3species.1Mb.rho_intergenic_df.png",width = 6, height = 6, units = 'in', res=300)
par(mfrow=c(3,2))
par(mar=c(4,4.5,1,1))

#pi
plot(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_intergenic_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[Intergenic]))
z=lm(total_1Mb$Thetas_tremula_intergenic_tP~total_1Mb$tremula_ldhat_new)
#r=round(cor(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_intergenic_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.002,0.03,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.002,0.027,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.002,0.027,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_trichocarpa_intergenic_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[Intergenic]))
z=lm(total_1Mb$tremula_trichocarpa_intergenic_fixed~total_1Mb$tremula_ldhat_new)
#r=round(cor(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_trichocarpa_intergenic_fixed),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.0018,0.0305,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0018,0.028,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0018,0.028,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

###tremuloides
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$Thetas_tremuloides_intergenic_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[Intergenic]))
z=lm(total_1Mb$Thetas_tremuloides_intergenic_tP~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(total_1Mb$tremuloides_ldhat_new,total_1Mb$Thetas_tremuloides_intergenic_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.027,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.0245,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.0245,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremuloides_trichocarpa_intergenic_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[Intergenic]))
z=lm(total_1Mb$tremuloides_trichocarpa_intergenic_fixed~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremuloides_trichocarpa_intergenic_fixed),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.029,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.06,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.026,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

#P.trichocarpa
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$Thetas_trichocarpa_intergenic_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[Intergenic]))
z=lm(total_1Mb$Thetas_trichocarpa_intergenic_tP~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)))
#r=round(cor(total_1Mb$trichocarpa_ldhat_new,total_1Mb$Thetas_trichocarpa_intergenic_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0014,0.017,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0014,0.015,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0014,0.015,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(c)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$tremula_trichocarpa_intergenic_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[Intergenic]))
z=lm(total_1Mb$tremula_trichocarpa_intergenic_fixed~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)))
#r=round(cor(total_1Mb$trichocarpa_ldhat_new,total_1Mb$tremula_trichocarpa_intergenic_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0014,0.031,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0014,0.0285,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0014,0.0285,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()

######test for Tajima's D
png(filename="3species.1Mb.rho_tajD_df.png",width = 6, height = 6, units = 'in', res=300)
par(mfrow=c(3,2))
par(mar=c(4,4.5,1,1))

#tajD
plot(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_four_fold_tajD,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(TajD[4-fold]))
z=lm(total_1Mb$Thetas_tremula_four_fold_tajD~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$Thetas_tremula_four_fold_tP,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.002,0.4,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.002,0.1,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.002,0.1,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[4-fold]))
z=lm(total_1Mb$tremula_trichocarpa_four_fold_fixed~total_1Mb$tremula_ldhat_new)
r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$tremula_trichocarpa_four_fold_fixed,use="complete"),2)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.0018,0.0275,bquote(paste(r,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0018,0.025,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0018,0.025,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))


###tremuloides
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$Thetas_tremuloides_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[4-fold]))
z=lm(total_1Mb$Thetas_tremuloides_four_fold_tP~total_1Mb$tremuloides_ldhat_new)
r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$Thetas_tremuloides_four_fold_tP,use="complete"),2)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.02,bquote(paste(r,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.018,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.018,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremuloides_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[4-fold]))
z=lm(total_1Mb$tremuloides_trichocarpa_four_fold_fixed~total_1Mb$tremuloides_ldhat_new)
r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$tremuloides_trichocarpa_four_fold_fixed,use="complete"),2)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.024,bquote(paste(r,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.022,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.022,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

#P.trichocarpa
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$Thetas_trichocarpa_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[4-fold]))
z=lm(total_1Mb$Thetas_trichocarpa_four_fold_tP~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)))
r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$Thetas_trichocarpa_four_fold_tP,use="complete"),2)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0014,0.0092,bquote(paste(r,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0014,0.008,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0014,0.008,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(c)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$tremula_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[4-fold]))
z=lm(total_1Mb$tremula_trichocarpa_four_fold_fixed~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),use="complete.obs")
r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$tremula_trichocarpa_four_fold_fixed,use="complete"),2)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0015,0.0265,bquote(paste(r,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0015,0.0245,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0015,0.0245,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()


######intron
###plots
png(filename="3species.1Mb.rho_intron_pi_df.png",width = 6, height = 6, units = 'in', res=300)
par(mfrow=c(3,2))
par(mar=c(4,4.5,1,1))

#pi
plot(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_intron_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[intron]))
z=lm(total_1Mb$Thetas_tremula_intron_tP~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$Thetas_tremula_intron_tP,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.002,0.015,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.002,0.013,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.002,0.013,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_trichocarpa_intron_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[intron]))
z=lm(total_1Mb$tremula_trichocarpa_intron_fixed~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$tremula_trichocarpa_intron_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.0018,0.02,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0018,0.0185,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0018,0.0185,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

###tremuloides
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$Thetas_tremuloides_intron_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[intron]))
z=lm(total_1Mb$Thetas_tremuloides_intron_tP~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$Thetas_tremuloides_intron_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.015,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.013,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.013,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremuloides_trichocarpa_intron_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[intron]))
z=lm(total_1Mb$tremuloides_trichocarpa_intron_fixed~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$tremuloides_trichocarpa_intron_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.018,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.0165,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.0165,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

#P.trichocarpa
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$Thetas_trichocarpa_intron_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[intron]))
z=lm(total_1Mb$Thetas_trichocarpa_intron_tP~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)))
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$Thetas_trichocarpa_intron_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0014,0.0078,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0014,0.0065,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0014,0.0065,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(c)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$tremula_trichocarpa_intron_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[intron]))
z=lm(total_1Mb$tremula_trichocarpa_intron_fixed~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),use="complete.obs")
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$tremula_trichocarpa_intron_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0015,0.02,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0015,0.0185,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0015,0.0185,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()


######0-fold
###plots
png(filename="3species.1Mb.rho_zero_fold_pi_df.png",width = 6, height = 6, units = 'in', res=300)
par(mfrow=c(3,2))
par(mar=c(4,4.5,1,1))

#pi
plot(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_zero_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[0-fold]))
z=lm(total_1Mb$Thetas_tremula_zero_fold_tP~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$Thetas_tremula_zero_fold_tP,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.002,0.0065,bquote(paste(r,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.002,0.0058,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.002,0.0058,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_trichocarpa_zero_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[0-fold]))
z=lm(total_1Mb$tremula_trichocarpa_zero_fold_fixed~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$tremula_trichocarpa_zero_fold_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.0018,0.0095,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0018,0.0087,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0018,0.0087,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))


###tremuloides
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$Thetas_tremuloides_zero_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[0-fold]))
z=lm(total_1Mb$Thetas_tremuloides_zero_fold_tP~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$Thetas_tremuloides_zero_fold_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.007,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.0062,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.0062,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremuloides_trichocarpa_zero_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[0-fold]))
z=lm(total_1Mb$tremuloides_trichocarpa_zero_fold_fixed~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$tremuloides_trichocarpa_zero_fold_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.0085,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.0076,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.0076,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

#P.trichocarpa
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$Thetas_trichocarpa_zero_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[0-fold]))
z=lm(total_1Mb$Thetas_trichocarpa_zero_fold_tP~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)))
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$Thetas_trichocarpa_zero_fold_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0014,0.0033,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0014,0.0028,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0014,0.0028,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(c)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$tremula_trichocarpa_zero_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[0-fold]))
z=lm(total_1Mb$tremula_trichocarpa_zero_fold_fixed~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),use="complete.obs")
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$tremula_trichocarpa_zero_fold_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0015,0.009,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0015,0.0081,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0015,0.0081,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()


######3UTR
###plots
png(filename="3species.1Mb.rho_3UTR_pi_df.png",width = 6, height = 6, units = 'in', res=300)
par(mfrow=c(3,2))
par(mar=c(4,4.5,1,1))

#pi
plot(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_3UTR_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[UTR3]))
z=lm(total_1Mb$Thetas_tremula_3UTR_tP~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$Thetas_tremula_3UTR_tP,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.002,0.015,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.002,0.013,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.002,0.013,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_trichocarpa_3UTR_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[UTR3]))
z=lm(total_1Mb$tremula_trichocarpa_3UTR_fixed~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$tremula_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.0018,0.0245,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0018,0.022,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0018,0.022,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

###tremuloides
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$Thetas_tremuloides_3UTR_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[UTR3]))
z=lm(total_1Mb$Thetas_tremuloides_3UTR_tP~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$Thetas_tremuloides_3UTR_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.018,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.0155,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.0155,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremuloides_trichocarpa_3UTR_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[UTR3]))
z=lm(total_1Mb$tremuloides_trichocarpa_3UTR_fixed~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$tremuloides_trichocarpa_3UTR_fixed,use="complete"),2)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.021,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.019,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.019,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

#P.trichocarpa
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$Thetas_trichocarpa_3UTR_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[UTR3]))
z=lm(total_1Mb$Thetas_trichocarpa_3UTR_tP~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)))
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$Thetas_trichocarpa_3UTR_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0014,0.0085,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0014,0.0075,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0014,0.0075,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(c)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$tremula_trichocarpa_3UTR_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[UTR3]))
z=lm(total_1Mb$tremula_trichocarpa_3UTR_fixed~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),use="complete.obs")
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$tremula_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0015,0.0235,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0015,0.0215,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0015,0.0215,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()


######5UTR
###plots
png(filename="3species.1Mb.rho_5UTR_pi_df.png",width = 6, height = 6, units = 'in', res=300)
par(mfrow=c(3,2))
par(mar=c(4,4.5,1,1))

#pi
plot(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_5UTR_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[UTR5]))
z=lm(total_1Mb$Thetas_tremula_5UTR_tP~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$Thetas_tremula_5UTR_tP,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.002,0.015,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.002,0.013,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.002,0.013,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_trichocarpa_5UTR_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[UTR5]))
z=lm(total_1Mb$tremula_trichocarpa_5UTR_fixed~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$tremula_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.0018,0.0245,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0018,0.022,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0018,0.022,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

###tremuloides
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$Thetas_tremuloides_5UTR_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[UTR5]))
z=lm(total_1Mb$Thetas_tremuloides_5UTR_tP~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$Thetas_tremuloides_3UTR_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.015,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.013,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.013,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremuloides_trichocarpa_5UTR_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[UTR5]))
z=lm(total_1Mb$tremuloides_trichocarpa_5UTR_fixed~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$tremuloides_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.022,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.02,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.02,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

#P.trichocarpa
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$Thetas_trichocarpa_5UTR_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[UTR5]))
z=lm(total_1Mb$Thetas_trichocarpa_5UTR_tP~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)))
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$Thetas_trichocarpa_3UTR_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0014,0.01,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0014,0.0085,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0014,0.0085,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(c)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$tremula_trichocarpa_5UTR_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[UTR5]))
z=lm(total_1Mb$tremula_trichocarpa_5UTR_fixed~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),use="complete.obs")
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$tremula_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0015,0.0235,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0015,0.021,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0015,0.021,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()


####pN/pS
png(filename="3species.1Mb.rho_HRI.png",width = 7, height = 6, units = 'in', res=300)
par(mfrow=c(3,4))
par(mar=c(2.5,4.5,1,0.5))

#pi
plot(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_zero_fold_tP/total_1Mb$Thetas_tremula_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[0-fold]/theta[4-fold]))
z=lm(total_1Mb$Thetas_tremula_zero_fold_tP/total_1Mb$Thetas_tremula_four_fold_tP~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$Thetas_tremula_5UTR_tP,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.002,0.6,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.002,0.53,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.002,0.53,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.4,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_trichocarpa_zero_fold_fixed/total_1Mb$tremula_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[0-fold]/d[4-fold]))
z=lm(total_1Mb$tremula_trichocarpa_zero_fold_fixed/total_1Mb$tremula_trichocarpa_four_fold_fixed~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$tremula_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.0018,0.49,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0018,0.45,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0018,0.45,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

#intergenic
plot(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_intergenic_tP/total_1Mb$Thetas_tremula_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[Intergenic]/theta[4-fold]))
z=lm(total_1Mb$Thetas_tremula_intergenic_tP/total_1Mb$Thetas_tremula_four_fold_tP~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$Thetas_tremula_5UTR_tP,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.002,6.5,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.002,5.6,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.002,0.53,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
#mtext("(a)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_trichocarpa_intergenic_fixed/total_1Mb$tremula_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[Intergenic]/d[4-fold]))
z=lm(total_1Mb$tremula_trichocarpa_intergenic_fixed/total_1Mb$tremula_trichocarpa_four_fold_fixed~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$tremula_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.0018,2,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0018,1.8,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0018,1.83,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))


###tremuloides
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$Thetas_tremuloides_zero_fold_tP/total_1Mb$Thetas_tremuloides_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[0-fold]/theta[4-fold]))
z=lm(total_1Mb$Thetas_tremuloides_zero_fold_tP/total_1Mb$Thetas_tremuloides_four_fold_tP~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$Thetas_tremuloides_3UTR_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.53,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.48,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.48,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.4,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremuloides_trichocarpa_zero_fold_fixed/total_1Mb$tremuloides_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[0-fold]/d[4-fold]))
z=lm(total_1Mb$tremuloides_trichocarpa_zero_fold_fixed/total_1Mb$tremuloides_trichocarpa_four_fold_fixed~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$tremuloides_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.48,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.44,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.44,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

#intergenic
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$Thetas_tremuloides_intergenic_tP/total_1Mb$Thetas_tremuloides_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[Intergenic]/theta[4-fold]))
z=lm(total_1Mb$Thetas_tremuloides_intergenic_tP/total_1Mb$Thetas_tremuloides_four_fold_tP~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$Thetas_tremuloides_3UTR_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,5,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,4.35,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.48,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
#mtext("(b)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremuloides_trichocarpa_intergenic_fixed/total_1Mb$tremuloides_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[Intergenic]/d[4-fold]))
z=lm(total_1Mb$tremuloides_trichocarpa_intergenic_fixed/total_1Mb$tremuloides_trichocarpa_four_fold_fixed~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$tremuloides_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,1.8,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,1.6,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,1.65,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

#P.trichocarpa
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$Thetas_trichocarpa_zero_fold_tP/total_1Mb$Thetas_trichocarpa_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[0-fold]/theta[4-fold]))
z=lm(total_1Mb$Thetas_trichocarpa_zero_fold_tP/total_1Mb$Thetas_trichocarpa_four_fold_tP~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)))
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$Thetas_trichocarpa_3UTR_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0014,1,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0014,0.83,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0014,0.87,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(c)",side=3,line=0.05,adj=-0.4,font=1.5,cex=1)

#fixed
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$tremula_trichocarpa_zero_fold_fixed/total_1Mb$tremula_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[0-fold]/d[4-fold]))
z=lm(total_1Mb$tremula_trichocarpa_zero_fold_fixed/total_1Mb$tremula_trichocarpa_four_fold_fixed~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),use="complete.obs")
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$tremula_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0015,0.49,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0015,0.45,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0015,0.45,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))


#intergenic
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$Thetas_trichocarpa_intergenic_tP/total_1Mb$Thetas_trichocarpa_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[Intergenic]/theta[4-fold]))
z=lm(total_1Mb$Thetas_trichocarpa_intergenic_tP/total_1Mb$Thetas_trichocarpa_four_fold_tP~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)))
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$Thetas_trichocarpa_3UTR_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0014,6.5,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0014,5.7,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0014,0.87,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
#mtext("(c)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$tremula_trichocarpa_intergenic_fixed/total_1Mb$tremula_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[Intergenic]/d[4-fold]))
z=lm(total_1Mb$tremula_trichocarpa_intergenic_fixed/total_1Mb$tremula_trichocarpa_four_fold_fixed~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),use="complete.obs")
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$tremula_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0015,2,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0015,0.45,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0015,1.85,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

mtext(expression(rho),side=1,line=1,adj=-2.1,font=1.5,cex=1)
dev.off()



####pIntergenic/pS
png(filename="3species.1Mb.rho_pNpS_dNdS.png",width = 6, height = 6, units = 'in', res=300)
par(mfrow=c(3,2))
par(mar=c(4,4.5,1,1))

#pi
plot(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_intergenic_tP/total_1Mb$Thetas_tremula_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[0-fold]/theta[4-fold]))
z=lm(total_1Mb$Thetas_tremula_intergenic_tP/total_1Mb$Thetas_tremula_four_fold_tP~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$Thetas_tremula_5UTR_tP,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.002,7,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.002,6,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.002,0.53,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_trichocarpa_intergenic_fixed/total_1Mb$tremula_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[0-fold]/d[4-fold]))
z=lm(total_1Mb$tremula_trichocarpa_intergenic_fixed/total_1Mb$tremula_trichocarpa_four_fold_fixed~total_1Mb$tremula_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$tremula_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.0018,2.1,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0018,1.9,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0018,1.9,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

###tremuloides
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$Thetas_tremuloides_intergenic_tP/total_1Mb$Thetas_tremuloides_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[0-fold]/theta[4-fold]))
z=lm(total_1Mb$Thetas_tremuloides_intergenic_tP/total_1Mb$Thetas_tremuloides_four_fold_tP~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$Thetas_tremuloides_3UTR_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,5.5,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,4.9,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.48,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremuloides_trichocarpa_intergenic_fixed/total_1Mb$tremuloides_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[0-fold]/d[4-fold]))
z=lm(total_1Mb$tremuloides_trichocarpa_intergenic_fixed/total_1Mb$tremuloides_trichocarpa_four_fold_fixed~total_1Mb$tremuloides_ldhat_new)
#r=round(cor(as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)),total_1Mb$tremuloides_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.005,1.9,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,1.8,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,1.7,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

#P.trichocarpa
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$Thetas_trichocarpa_intergenic_tP/total_1Mb$Thetas_trichocarpa_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(theta[0-fold]/theta[4-fold]))
z=lm(total_1Mb$Thetas_trichocarpa_intergenic_tP/total_1Mb$Thetas_trichocarpa_four_fold_tP~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)))
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$Thetas_trichocarpa_3UTR_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0014,7,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0014,6,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0014,0.87,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(c)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#fixed
plot(total_1Mb$trichocarpa_ldhat_new,total_1Mb$tremula_trichocarpa_intergenic_fixed/total_1Mb$tremula_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(rho),ylab=expression(d[0-fold]/d[4-fold]))
z=lm(total_1Mb$tremula_trichocarpa_intergenic_fixed/total_1Mb$tremula_trichocarpa_four_fold_fixed~as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),use="complete.obs")
#r=round(cor(as.numeric(as.character(total_1Mb$trichocarpa_ldhat_new)),total_1Mb$tremula_trichocarpa_3UTR_fixed,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0015,2.1,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0015,0.45,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0015,1.9,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()


####Relationship between dxy in 0-fold and diveristy in 4-fold or intergenic
###plots
png(filename="3species.1Mb.0_4_intergenic.df.png",width = 6, height = 6, units = 'in', res=300)
par(mfrow=c(3,2))
par(mar=c(4,4.5,1,1))

#pi,4-fold
plot(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremula_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(d[0-fold]),ylab=expression(theta[4-fold]))
z=lm(total_1Mb$Thetas_tremula_four_fold_tP~total_1Mb$tremula_trichocarpa_zero_fold_fixed)
#r=round(cor(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremula_four_fold_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.0051,0.015,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0051,0.013,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0051,0.013,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#pi,intergenic
plot(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremula_intergenic_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(d[0-fold]),ylab=expression(theta[Intergenic]))
z=lm(total_1Mb$Thetas_tremula_intergenic_tP~total_1Mb$tremula_trichocarpa_zero_fold_fixed)
#r=round(cor(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremula_intergenic_tP),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(0.0051,0.029,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0051,0.026,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0051,0.026,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))


###tremuloides
plot(total_1Mb$tremuloides_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremuloides_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(d[0-fold]),ylab=expression(theta[4-fold]))
z=lm(total_1Mb$Thetas_tremuloides_four_fold_tP~total_1Mb$tremuloides_trichocarpa_zero_fold_fixed)
#r=round(cor(total_1Mb$tremuloides_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremuloides_four_fold_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.0047,0.018,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0047,0.016,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0047,0.016,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#pi,intergenic
plot(total_1Mb$tremuloides_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremuloides_intergenic_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(d[0-fold]),ylab=expression(theta[Intergenic]))
z=lm(total_1Mb$Thetas_tremuloides_intergenic_tP~total_1Mb$tremuloides_trichocarpa_zero_fold_fixed)
#r=round(cor(total_1Mb$tremuloides_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_tremuloides_intergenic_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.0047,0.026,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0047,0.024,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0047,0.024,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

###trichocarpa
plot(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_trichocarpa_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(d[0-fold]),ylab=expression(theta[4-fold]))
z=lm(total_1Mb$Thetas_trichocarpa_four_fold_tP~total_1Mb$tremula_trichocarpa_zero_fold_fixed)
#r=round(cor(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_trichocarpa_four_fold_tP,use="complete"),4)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0051,0.008,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0051,0.007,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0051,0.007,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(c)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#pi,intergenic
plot(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_trichocarpa_intergenic_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(d[0-fold]),ylab=expression(theta[Intergenic]))
z=lm(total_1Mb$Thetas_trichocarpa_intergenic_tP~total_1Mb$tremula_trichocarpa_zero_fold_fixed)
#r=round(cor(total_1Mb$tremula_trichocarpa_zero_fold_fixed,total_1Mb$Thetas_trichocarpa_intergenic_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(0.0051,0.0155,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.0051,0.014,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.0051,0.014,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()
 

###Correlation betweeen gene denisty and dIntergenic and d4-fold

###plots
png(filename="3species.1Mb.gene_density_divergence.png",width = 6, height = 6, units = 'in', res=300)
par(mfrow=c(3,2))
par(mar=c(4,4.5,1,1))

#tremula
#d4-fold
plot(total_1Mb$Gene_num,total_1Mb$tremula_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab="Gene num",ylab=expression(d[4-fold]))
z=lm(total_1Mb$tremula_trichocarpa_four_fold_fixed~total_1Mb$Gene_num)
#r=round(cor(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_intergenic_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(50,0.025,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(50,0.022,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(50,0.022,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#dIntergenic
plot(total_1Mb$Gene_num,total_1Mb$tremula_trichocarpa_intergenic_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab="Gene num",ylab=expression(d[Intergenic]))
z=lm(total_1Mb$tremula_trichocarpa_intergenic_fixed~total_1Mb$Gene_num)
#r=round(cor(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_intergenic_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[1],lwd=2) # equivalent to abline(reg = z) or
text(50,0.03,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(50,0.027,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(50,0.027,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

###tremuloides
#d4-fold
plot(total_1Mb$Gene_num,total_1Mb$tremuloides_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab="Gene num",ylab=expression(d[4-fold]))
z=lm(total_1Mb$tremuloides_trichocarpa_four_fold_fixed~total_1Mb$Gene_num)
#r=round(cor(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_intergenic_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(50,0.023,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(50,0.02,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(50,0.021,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#dIntergenic
plot(total_1Mb$Gene_num,total_1Mb$tremuloides_trichocarpa_intergenic_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab="Gene num",ylab=expression(d[Intergenic]))
z=lm(total_1Mb$tremuloides_trichocarpa_intergenic_fixed~total_1Mb$Gene_num)
#r=round(cor(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_intergenic_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(50,0.029,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(50,0.027,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(50,0.027,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))

#trichocarpa
#d4-fold
plot(total_1Mb$Gene_num,total_1Mb$tremula_trichocarpa_four_fold_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab="Gene num",ylab=expression(d[4-fold]))
z=lm(total_1Mb$tremula_trichocarpa_four_fold_fixed~total_1Mb$Gene_num)
#r=round(cor(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_intergenic_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(50,0.025,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(50,0.022,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(50,0.022,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(c)",side=3,line=0.05,adj=-0.28,font=1.5,cex=1)

#dIntergenic
plot(total_1Mb$Gene_num,total_1Mb$tremula_trichocarpa_intergenic_fixed,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab="Gene num",ylab=expression(d[Intergenic]))
z=lm(total_1Mb$tremula_trichocarpa_intergenic_fixed~total_1Mb$Gene_num)
#r=round(cor(total_1Mb$tremula_ldhat_new,total_1Mb$Thetas_tremula_intergenic_tP,use="complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col=colors[3],lwd=2) # equivalent to abline(reg = z) or
text(50,0.03,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(50,0.027,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(50,0.027,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()







#####for speciation paper
library(RColorBrewer)
colors <- brewer.pal(10,"Paired")[c(1,2,3,4,5,6,7,8,9,10,11,12)]


setwd("~/Dropbox/combined_paper/paper/version3/plots")
png(filename="1Mb.rho_fst.png",width = 6, height = 4, units = 'in', res=300)

par(mfrow=c(1,2))
par(mar=c(4,4,1.5,1))
plot(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_tremuloides_fst,pch=19,col="grey",cex=0.7,xlab=expression(paste(rho[(P.tremula)])),ylab=expression(F[ST]))
z=lm(total_1Mb$tremula_tremuloides_fst~as.numeric(as.character(total_1Mb$tremula_ldhat_new)))
r=round(cor(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_tremuloides_fst),4)
#r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,4)
abline(z,col=colors[6],lwd=2) # equivalent to abline(reg = z) or
#text(0.005,0.55,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
text(0.005,0.55,bquote(paste(r,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.5,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.5,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.3,font=1.5,cex=1.5)

plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremula_tremuloides_fst,pch=19,col="grey",cex=0.7,xlab=expression(paste(rho[(P.tremuloides)])),ylab=expression(F[ST]))
z=lm(total_1Mb$tremula_tremuloides_fst~as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)))
r=round(cor(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremula_tremuloides_fst),4)
#r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,4)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.013,0.55,bquote(paste(r,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.013,0.5,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.013,0.5,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.3,font=1.5,cex=1.5)
dev.off()

####dxy
png(filename="1Mb.rho_dxy.png",width = 5.5, height = 7, units = 'in', res=300)

par(mfrow=c(2,1))
par(mar=c(4,4.5,1.5,1))
plot(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_tremuloides_dxy,pch=19,col="grey",cex=0.8,xlab=expression(paste(rho[(P.tremula)])),ylab=expression(d[xy]))
z=lm(total_1Mb$tremula_tremuloides_dxy~as.numeric(as.character(total_1Mb$tremula_ldhat_new)))
#r=round(cor(total_1Mb$tremula_ldhat_new,total_1Mb$tremula_tremuloides_fst),4)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,4)
abline(z,col=colors[6],lwd=2) # equivalent to abline(reg = z) or
text(0.005,0.03,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.005,0.028,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.005,0.028,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(a)",side=3,line=0.05,adj=-0.15,font=1.5,cex=1.5)

plot(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremula_tremuloides_dxy,pch=19,col="grey",cex=0.8,xlab=expression(paste(rho[(P.tremuloides)])),ylab=expression(d[xy]))
z=lm(total_1Mb$tremula_tremuloides_dxy~as.numeric(as.character(total_1Mb$tremuloides_ldhat_new)))
#r=round(cor(total_1Mb$tremuloides_ldhat_new,total_1Mb$tremula_tremuloides_fst),4)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,4)
abline(z,col=colors[2],lwd=2) # equivalent to abline(reg = z) or
text(0.013,0.03,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(0.013,0.028,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(0.013,0.028,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
mtext("(b)",side=3,line=0.05,adj=-0.15,font=1.5,cex=1.5)
dev.off()


####linear model of 4-fold diveristy~4-fold divergence+0-fold divergence
tremula_0_4_fold=lm(total_1Mb$Thetas_tremula_four_fold_tP~total_1Mb$tremula_trichocarpa_four_fold_fixed+total_1Mb$tremula_trichocarpa_zero_fold_fixed)
summary(tremula_0_4_fold)
tremuloides_0_4_fold=lm(total_1Mb$Thetas_tremuloides_four_fold_tP~total_1Mb$tremuloides_trichocarpa_four_fold_fixed+total_1Mb$tremuloides_trichocarpa_zero_fold_fixed)
summary(tremuloides_0_4_fold)
trichocarpa_0_4_fold=lm(total_1Mb$Thetas_trichocarpa_four_fold_tP~total_1Mb$tremula_trichocarpa_four_fold_fixed+total_1Mb$tremula_trichocarpa_zero_fold_fixed)
summary(trichocarpa_0_4_fold)





