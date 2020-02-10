library(RColorBrewer)
colors <- brewer.pal(9,"Set1")[c(5,2,3)]

setwd("~/Dropbox/paperI/data/100kb")
total_new=read.table("100kb.summary.txt",header=T)

total_new$tremula_ldhat_new=total_new$Ldhat_tremula/1000
total_new$tremuloides_ldhat_new=total_new$Ldhat_tremuloides/1000
total_new$trichocarpa_ldhat_new=total_new$Ldhat_trichocarpa/1000

total_new[which(total_new$Ldhat_tremula_num<50),]$tremula_ldhat_new="NA"
total_new[which(total_new$Ldhat_tremuloides_num<50),]$tremuloides_ldhat_new="NA"
total_new[which(total_new$Ldhat_trichocarpa_num<50),]$trichocarpa_ldhat_new="NA"

total_new[which(total_new$Ldhat_tremula_num<50),]$Ld_tremula="NA"
total_new[which(total_new$Ldhat_tremuloides_num<50),]$Ld_tremuloides="NA"
total_new[which(total_new$Ldhat_trichocarpa_num<50),]$Ld_trichocarpa="NA"

pal <- colorRampPalette(c("light blue", "yellow", "red"))

mat=matrix(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3),rep(4,3),rep(5,3),rep(6,3),rep(4,3),rep(5,3),rep(6,3),rep(7,3),rep(8,3),rep(9,3),rep(10,3),rep(11,3),rep(12,3),rep(10,3),rep(11,3),rep(12,3),rep(10,3),rep(11,3),rep(12,3),rep(13,3),rep(14,3),rep(15,3),rep(16,3),rep(17,3),rep(18,3),rep(16,3),rep(17,3),rep(18,3),rep(16,3),rep(17,3),rep(18,3)),byrow=T,nrow=12)

layout(mat)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}



opa = par(tcl=-0.3, mgp=c(2, 0.6, 0))
opa
par(opa)


png("3species.tP.tajD.rho.cor.100kb.png",width = 6, height = 6, units = 'in',res=500)
par(mfrow=c(3,3))
par(mar=c(4,4,1,1))
par(tcl=-0.3,mgp=c(2,0.4,0))
par(asp=1)


colors_pi=densCols(total_new$Thetas_tremula_tP,total_new$Thetas_tremuloides_tP,colramp=pal)
plot(total_new$Thetas_tremula_tP,total_new$Thetas_tremuloides_tP,pch=19,col=colors_pi,cex=.5,cex.lab=1,cex.axis=0.7,xlim=c(0,0.08),ylim=c(0,0.08),xlab=expression(paste(italic(P.tremula)," (",theta[pi],")")),ylab=expression(paste(italic(P.tremuloides)," (",theta[pi],")")))
abline(0,1,lty=3,lwd=0.5)
cor.test(total_new$Thetas_tremula_tP,total_new$Thetas_tremuloides_tP,method="spearman")
#text(0.0175,0.06,expression(paste(rho,"=0.829")^"***"))
#legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression(paste(rho,"=0.829")^"***"))
legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression("rho=0.829"^"***"))
mtext("(a)",side=3,line=0.05,adj=-0.4,font=2,cex=1)

par(mar=c(4,4,1,1))
colors_pi=densCols(total_new$Thetas_tremula_tP,total_new$Thetas_trichocarpa_tP,colramp=pal)
plot(total_new$Thetas_tremula_tP,total_new$Thetas_trichocarpa_tP,pch=19,col=colors_pi,cex=.5,cex.lab=1,cex.axis=0.7,xlim=c(0,0.08),ylim=c(0,0.08),xlab=expression(paste(italic(P.tremula)," (",theta[pi],")")),ylab=expression(paste(italic(P.trichocarpa)," (",theta[pi],")")))
abline(0,1,lty=3,lwd=0.5)
cor.test(total_new$Thetas_tremula_tP,total_new$Thetas_trichocarpa_tP,method="spearman")
#text(0.0175,0.06,expression(paste(rho,"=0.675")^"***"))
#legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression(paste(rho,"=0.675")^"***"))
legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression("rho=0.675"^"***"))


par(mar=c(4,4,1,1))
colors_pi=densCols(total_new$Thetas_tremuloides_tP,total_new$Thetas_trichocarpa_tP,colramp=pal)
plot(total_new$Thetas_tremuloides_tP,total_new$Thetas_trichocarpa_tP,pch=19,col=colors_pi,cex=.5,cex.lab=1,cex.axis=0.7,xlim=c(0,0.08),ylim=c(0,0.08),xlab=expression(paste(italic(P.tremuloides)," (",theta[pi],")")),ylab=expression(paste(italic(P.trichocarpa)," (",theta[pi],")")))
abline(0,1,lty=3,lwd=0.5)
cor.test(total_new$Thetas_tremuloides_tP,total_new$Thetas_trichocarpa_tP,method="spearman")
#text(0.0175,0.06,expression(paste(rho,"=0.658")^"***"))
#legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression(paste(rho,"=0.658")^"***"))
legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression("rho=0.658"^"***"))


par(mar=c(4,4,1,1))
colors_pi=densCols(total_new$Thetas_tremula_tajD,total_new$Thetas_tremuloides_tajD,colramp=pal)
plot(total_new$Thetas_tremula_tajD,total_new$Thetas_tremuloides_tajD,pch=19,col=colors_pi,cex=.5,cex.lab=1,cex.axis=0.7,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),xlab=expression(paste(italic(P.tremula)," (TajD)")),ylab=expression(paste(italic(P.tremuloides)," (TajD)")))
abline(0,1,lty=3,lwd=0.5)
cor.test(total_new$Thetas_tremula_tajD,total_new$Thetas_tremuloides_tajD,method="spearman")
#text(-1.55,1.2,expression(paste(rho,"=0.319")^"***"))
mtext("(b)",side=3,line=0.05,adj=-0.4,font=2,cex=1)
#legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression(paste(rho,"=0.319")^"***"))
legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression("rho=0.319"^"***"))


par(mar=c(4,4,1,1))
colors_pi=densCols(total_new$Thetas_tremula_tajD,total_new$Thetas_trichocarpa_tajD,colramp=pal)
plot(total_new$Thetas_tremula_tajD,total_new$Thetas_trichocarpa_tajD,pch=19,col=colors_pi,cex=.5,cex.lab=1,cex.axis=0.7,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),xlab=expression(paste(italic(P.tremula)," (TajD)")),ylab=expression(paste(italic(P.trichocarpa)," (TajD)")))
abline(0,1,lty=3,lwd=0.5)
cor.test(total_new$Thetas_tremula_tajD,total_new$Thetas_trichocarpa_tajD,method="spearman")
#text(-1.55,1.2,expression(paste(rho,"=0.099")^"**"))
#legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression(paste(rho,"=0.099")^"**"))
legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression("rho=0.099"^"**"))


par(mar=c(4,4,1,1))
colors_pi=densCols(total_new$Thetas_tremuloides_tajD,total_new$Thetas_trichocarpa_tajD,colramp=pal)
plot(total_new$Thetas_tremuloides_tajD,total_new$Thetas_trichocarpa_tajD,pch=19,col=colors_pi,cex=.5,cex.lab=1,cex.axis=0.7,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),xlab=expression(paste(italic(P.tremuloides)," (TajD)")),ylab=expression(paste(italic(P.trichocarpa)," (TajD)")))
abline(0,1,lty=3,lwd=0.5)
cor.test(total_new$Thetas_tremuloides_tajD,total_new$Thetas_trichocarpa_tajD,method="spearman")
#text(-1.55,1.2,expression(paste(rho,"=0.164")^"***"))
#legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression(paste(rho,"=0.164")^"***"))
legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression("rho=0.164"^"***"))


total_new$tremula_ldhat_new=as.numeric(as.character(total_new$tremula_ldhat_new))
total_new$tremuloides_ldhat_new=as.numeric(as.character(total_new$tremuloides_ldhat_new))
total_new$trichocarpa_ldhat_new=as.numeric(as.character(total_new$trichocarpa_ldhat_new))
total_new$Ld_tremula=as.numeric(as.character(total_new$Ld_tremula))
total_new$Ld_tremuloides=as.numeric(as.character(total_new$Ld_tremuloides))
total_new$Ld_trichocarpa=as.numeric(as.character(total_new$Ld_trichocarpa))

par(mar=c(4,4,1,1))
colors_pi=densCols(total_new$tremula_ldhat_new,total_new$tremuloides_ldhat_new,colramp=pal)
plot(total_new$tremula_ldhat_new,total_new$tremuloides_ldhat_new,pch=19,col=colors_pi,cex=.5,cex.lab=1,cex.axis=0.7,xlim=c(0,0.035),ylim=c(0,0.035),xlab=expression(paste(italic(P.tremula)," (",rho,")")),ylab=expression(paste(italic(P.tremuloides)," (",rho,")")))
abline(0,1,lty=3,lwd=0.5)
cor.test(total_new$tremula_ldhat_new,total_new$tremuloides_ldhat_new,method="spearman")
#text(0.0075,0.026,expression(paste(rho,"=0.514")^"***"))
mtext("(c)",side=3,line=0.05,adj=-0.4,font=2,cex=1)
#legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression(paste(rho,"=0.514")^"***"))
legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression("rho=0.514"^"***"))


par(mar=c(4,4,1,1))
colors_pi=densCols(total_new$tremula_ldhat_new,total_new$trichocarpa_ldhat_new,colramp=pal)
plot(total_new$tremula_ldhat_new,total_new$trichocarpa_ldhat_new,pch=19,col=colors_pi,cex=.5,cex.lab=1,cex.axis=0.7,xlim=c(0,0.035),ylim=c(0,0.035),xlab=expression(paste(italic(P.tremula)," (",rho,")")),ylab=expression(paste(italic(P.trichocarpa)," (",rho,")")))
abline(0,1,lty=3,lwd=0.5)
cor.test(total_new$tremula_ldhat_new,total_new$trichocarpa_ldhat_new,method="spearman")
#text(0.0075,0.026,expression(paste(rho,"=0.317")^"***"))
#legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression(paste(rho,"=0.317")^"***"))
legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression("rho=0.317"^"***"))

par(mar=c(4,4,1,1))
colors_pi=densCols(total_new$tremuloides_ldhat_new,total_new$trichocarpa_ldhat_new,colramp=pal)
plot(total_new$tremuloides_ldhat_new,total_new$trichocarpa_ldhat_new,pch=19,col=colors_pi,cex=.5,cex.lab=1,cex.axis=0.7,xlim=c(0,0.035),ylim=c(0,0.035),xlab=expression(paste(italic(P.tremuloides)," (",rho,")")),ylab=expression(paste(italic(P.trichocarpa)," (",rho,")")))
abline(0,1,lty=3,lwd=0.5)
cor.test(total_new$tremuloides_ldhat_new,total_new$trichocarpa_ldhat_new,method="spearman")
#text(0.0075,0.026,expression(paste(rho,"=0.306")^"***"))
#legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression(paste(rho,"=0.306")^"***"))
legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression("rho=0.306"^"***"))

dev.off()

