library(RColorBrewer)
colors <- brewer.pal(9,"Set1")[c(5,2,3)]

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


total_1Mb_new=total_1Mb[complete.cases(total_1Mb),]

###plot1: correlation between gene number and coding density
tiff(filename="Gene_num_vs_Gene_density.tiff",width = 4, height = 4, units = 'in', res=300)
par(mfcol=c(1,1))
plot(total_1Mb$Gene_num,total_1Mb$Coding_prop,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab="Gene number",ylab="Coding density")
z=lm(total_1Mb$Coding_prop~total_1Mb$Gene_num)
#r=round(cor(as.numeric(as.character(total_1Mb$tremula_ldhat_new)),total_1Mb$Gene_num,use = "complete"),2)
r=round(summary(z)$adj.r.squared,3)
p=round(summary(z)$coefficients[,4][2] ,3)
abline(z,col="black",lwd=1) # equivalent to abline(reg = z) or
text(50,0.2,bquote(paste(R^2,"=",.(r),sep="")),cex=1)
ifelse(p<0.001,text(50,0.18,expression(paste(italic(P),"<0.001",sep="")),cex=1),text(50,0.0047,bquote(paste(italic(P),"=",.(p),sep="")),cex=1))
dev.off()


###tremula
tremula_cor_parameters=total_1Mb_new[c("GC","Gene_num","tremula_ldhat_new","Thetas_tremula_four_fold_tP","tremula_trichocarpa_four_fold_fixed","Thetas_tremula_intergenic_tP","tremula_trichocarpa_intergenic_fixed","tremula_trichocarpa_zero_fold_fixed")]


###tremuloides
tremuloides_cor_parameters=total_1Mb_new[c("GC","Gene_num","tremuloides_ldhat_new","Thetas_tremuloides_four_fold_tP","tremuloides_trichocarpa_four_fold_fixed","Thetas_tremuloides_intergenic_tP","tremuloides_trichocarpa_intergenic_fixed","tremuloides_trichocarpa_zero_fold_fixed")]

###trichocarpa
trichocarpa_cor_parameters=total_1Mb_new[c("GC","Gene_num","trichocarpa_ldhat_new","Thetas_trichocarpa_four_fold_tP","tremula_trichocarpa_four_fold_fixed","Thetas_trichocarpa_intergenic_tP","tremula_trichocarpa_intergenic_fixed","tremula_trichocarpa_zero_fold_fixed")]
trichocarpa_cor_parameters$trichocarpa_ldhat_new=as.numeric(as.character(trichocarpa_cor_parameters$trichocarpa_ldhat_new))


###plots
panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y,method="spearman")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y,method="spearman")
  Signif <- ifelse(round(test$p.value,3)<0.001,"P<0.001",paste("P=",round(test$p.value,3)))  
  text(0.5, 0.5,cex=0.8,paste("r=",txt))
  text(.5, .75,cex=0.8,Signif)
}
panel.smooth<-function (x, y, col = "gray", bg = NA, pch = 18, 
                        cex = 0.6, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE, breaks=15)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="blue", ...)
}

tremula_cor_parameters_new=tremula_cor_parameters[,c(2,3,4,6)]
tremula_cor_parameters_new_small=tremula_cor_parameters_new[which(tremula_cor_parameters_new$Gene_num<85),]
tremula_cor_parameters_new_large=tremula_cor_parameters_new[which(tremula_cor_parameters_new$Gene_num>=85),]

cor.test(tremula_cor_parameters_new_small$Gene_num,tremula_cor_parameters_new_small$tremula_ldhat_new,method="spearman")
cor.test(tremula_cor_parameters_new_large$Gene_num,tremula_cor_parameters_new_large$tremula_ldhat_new,method="spearman")

cor.test(tremula_cor_parameters_new_small$Gene_num,tremula_cor_parameters_new_small$Thetas_tremula_four_fold_tP,method="spearman")
cor.test(tremula_cor_parameters_new_large$Gene_num,tremula_cor_parameters_new_large$Thetas_tremula_four_fold_tP,method="spearman")

cor.test(tremula_cor_parameters_new_small$Gene_num,tremula_cor_parameters_new_small$Thetas_tremula_intergenic_tP,method="spearman")
cor.test(tremula_cor_parameters_new_large$Gene_num,tremula_cor_parameters_new_large$Thetas_tremula_intergenic_tP,method="spearman")


tremuloides_cor_parameters_new=tremuloides_cor_parameters[,c(2,3,4,6)]
tremuloides_cor_parameters_new_small=tremuloides_cor_parameters_new[which(tremuloides_cor_parameters_new$Gene_num<85),]
tremuloides_cor_parameters_new_large=tremuloides_cor_parameters_new[which(tremuloides_cor_parameters_new$Gene_num>=85),]

cor.test(tremuloides_cor_parameters_new_small$Gene_num,tremuloides_cor_parameters_new_small$tremuloides_ldhat_new,method="spearman")
cor.test(tremuloides_cor_parameters_new_large$Gene_num,tremuloides_cor_parameters_new_large$tremuloides_ldhat_new,method="spearman")

cor.test(tremuloides_cor_parameters_new_small$Gene_num,tremuloides_cor_parameters_new_small$Thetas_tremuloides_four_fold_tP,method="spearman")
cor.test(tremuloides_cor_parameters_new_large$Gene_num,tremuloides_cor_parameters_new_large$Thetas_tremuloides_four_fold_tP,method="spearman")

cor.test(tremuloides_cor_parameters_new_small$Gene_num,tremuloides_cor_parameters_new_small$Thetas_tremuloides_intergenic_tP,method="spearman")
cor.test(tremuloides_cor_parameters_new_large$Gene_num,tremuloides_cor_parameters_new_large$Thetas_tremuloides_intergenic_tP,method="spearman")

trichocarpa_cor_parameters_new=trichocarpa_cor_parameters[,c(2,3,4,6)]
trichocarpa_cor_parameters_new_small=trichocarpa_cor_parameters_new[which(trichocarpa_cor_parameters_new$Gene_num<85),]
trichocarpa_cor_parameters_new_large=trichocarpa_cor_parameters_new[which(trichocarpa_cor_parameters_new$Gene_num>=85),]

cor.test(trichocarpa_cor_parameters_new_small$Gene_num,trichocarpa_cor_parameters_new_small$trichocarpa_ldhat_new,method="spearman")
cor.test(trichocarpa_cor_parameters_new_large$Gene_num,trichocarpa_cor_parameters_new_large$trichocarpa_ldhat_new,method="spearman")

cor.test(trichocarpa_cor_parameters_new_small$Gene_num,trichocarpa_cor_parameters_new_small$Thetas_trichocarpa_four_fold_tP,method="spearman")
cor.test(trichocarpa_cor_parameters_new_large$Gene_num,trichocarpa_cor_parameters_new_large$Thetas_trichocarpa_four_fold_tP,method="spearman")

cor.test(trichocarpa_cor_parameters_new_small$Gene_num,trichocarpa_cor_parameters_new_small$Thetas_trichocarpa_intergenic_tP,method="spearman")
cor.test(trichocarpa_cor_parameters_new_large$Gene_num,trichocarpa_cor_parameters_new_large$Thetas_trichocarpa_intergenic_tP,method="spearman")


png(filename="3species.1Mb.gene_density.correlation.png",width = 7, height = 7, units = 'in', res=300)


par(mfcol=c(3,3))
par(mar=c(2,4.5,2,1))
par(tcl=-0.3, mgp=c(2, 0.6, 0))
#plot(tremula_cor_parameters_new$Gene_num,tremula_cor_parameters_new$tremula_ldhat_new,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab=expression(paste("P.tremula","(",rho,")",sep=" ")),xlab="Gene number")
plot(tremula_cor_parameters_new$Gene_num,tremula_cor_parameters_new$tremula_ldhat_new,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,xlab="",ylab="")
abline(v=85,lty=2)
panel.smooth(tremula_cor_parameters_new$Gene_num,tremula_cor_parameters_new$tremula_ldhat_new,col.smooth = colors[1], span = 2/3,iter = 10)
mtext("(a)",side=3,line=0.05,adj=-0.35,font=1.5,cex=1.2)
mtext("P.tremula",side=3,line=0.05,adj=.45,font=1.5,cex=1.2)
mtext(expression(rho),side=3,line=-5.5,adj=-0.25,font=1,cex=1)

#plot(tremula_cor_parameters_new$Gene_num,tremula_cor_parameters_new$Thetas_tremula_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab=expression(paste("P.tremula","(",theta[4-fold],")",sep=" ")),xlab="Gene number")
plot(tremula_cor_parameters_new$Gene_num,tremula_cor_parameters_new$Thetas_tremula_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab="",xlab="")
abline(v=85,lty=2)
panel.smooth(tremula_cor_parameters_new$Gene_num,tremula_cor_parameters_new$Thetas_tremula_four_fold_tP,col.smooth = colors[1], span = 2/3,iter = 10)
mtext("(d)",side=3,line=0.05,adj=-0.35,font=1.5,cex=1.2)
mtext(expression(theta[4-fold]),side=3,line=-5.5,adj=-0.4,font=1,cex=1)


plot(tremula_cor_parameters_new$Gene_num,tremula_cor_parameters_new$Thetas_tremula_intergenic_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab=expression(paste("P.tremula","(",theta[Intergenic],")",sep=" ")),xlab="Gene number")
abline(v=85,lty=2)
panel.smooth(tremula_cor_parameters_new$Gene_num,tremula_cor_parameters_new$Thetas_tremula_intergenic_tP,col.smooth = colors[1], span = 2/3,iter = 10)
mtext("(g)",side=3,line=0.05,adj=-0.35,font=1.5,cex=1.2)


###tremuloides
plot(tremuloides_cor_parameters_new$Gene_num,tremuloides_cor_parameters_new$tremuloides_ldhat_new,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab=expression(paste("P.tremuloides","(",rho,")",sep=" ")),xlab="Gene number")
abline(v=85,lty=2)
panel.smooth(tremuloides_cor_parameters_new$Gene_num,tremuloides_cor_parameters_new$tremuloides_ldhat_new,col.smooth = colors[2], span = 2/3,iter = 10)
mtext("(b)",side=3,line=0.05,adj=-0.35,font=1.5,cex=1.2)


plot(tremuloides_cor_parameters_new$Gene_num,tremuloides_cor_parameters_new$Thetas_tremuloides_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab=expression(paste("P.tremuloides","(",theta[4-fold],")",sep=" ")),xlab="Gene number")
abline(v=85,lty=2)
panel.smooth(tremuloides_cor_parameters_new$Gene_num,tremuloides_cor_parameters_new$Thetas_tremuloides_four_fold_tP,col.smooth = colors[2], span = 2/3,iter = 10)
mtext("(e)",side=3,line=0.05,adj=-0.35,font=1.5,cex=1.2)

plot(tremuloides_cor_parameters_new$Gene_num,tremuloides_cor_parameters_new$Thetas_tremuloides_intergenic_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab=expression(paste("P.tremuloides","(",theta[Intergenic],")",sep=" ")),xlab="Gene number")
abline(v=85,lty=2)
panel.smooth(tremuloides_cor_parameters_new$Gene_num,tremuloides_cor_parameters_new$Thetas_tremuloides_intergenic_tP,col.smooth = colors[2], span = 2/3,iter = 10)
mtext("(h)",side=3,line=0.05,adj=-0.35,font=1.5,cex=1.2)

###trichocarpa
plot(trichocarpa_cor_parameters_new$Gene_num,trichocarpa_cor_parameters_new$trichocarpa_ldhat_new,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab=expression(paste("P.trichocarpa","(",rho,")",sep=" ")),xlab="Gene number")
abline(v=85,lty=2)
panel.smooth(trichocarpa_cor_parameters_new$Gene_num,trichocarpa_cor_parameters_new$trichocarpa_ldhat_new,col.smooth = colors[3], span = 2/3,iter = 10)
mtext("(c)",side=3,line=0.05,adj=-0.35,font=1.5,cex=1.2)

plot(trichocarpa_cor_parameters_new$Gene_num,trichocarpa_cor_parameters_new$Thetas_trichocarpa_four_fold_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab=expression(paste("P.trichocarpa","(",theta[4-fold],")",sep=" ")),xlab="Gene number")
abline(v=85,lty=2)
panel.smooth(trichocarpa_cor_parameters_new$Gene_num,trichocarpa_cor_parameters_new$Thetas_trichocarpa_four_fold_tP,col.smooth = colors[3], span = 2/3,iter = 10)
mtext("(f)",side=3,line=0.05,adj=-0.35,font=1.5,cex=1.2)

plot(trichocarpa_cor_parameters_new$Gene_num,trichocarpa_cor_parameters_new$Thetas_trichocarpa_intergenic_tP,pch=19,col="grey",cex=.5,cex.lab=1,cex.axis=0.7,ylab=expression(paste("P.trichocarpa","(",theta[Intergenic],")",sep=" ")),xlab="Gene number")
abline(v=85,lty=2)
panel.smooth(trichocarpa_cor_parameters_new$Gene_num,trichocarpa_cor_parameters_new$Thetas_trichocarpa_intergenic_tP,col.smooth = colors[3], span = 2/3,iter = 10)
mtext("(i)",side=3,line=0.05,adj=-0.35,font=1.5,cex=1.2)
dev.off()



###plots for tremula
png(filename="tremula.1Mb.cor.gene.png",width = 7, height = 7, units = 'in', res=300)
pairs(tremula_cor_parameters,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("GC","Gene num",expression(rho),expression(theta[4-fold]),expression(d[4-fold]),expression(theta[Intergenic]),expression(d[Intergenic]),expression(d[0-fold])))
#pairs(tremula_cor_parameters_3,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("GC","Gene num",expression(rho),expression(theta[4-fold]),expression(TajD[4-fold]),expression(d[4-fold]),expression(theta[Intergenic]),expression(TajD[Intergenic]),expression(d[Intergenic]),expression(d[0-fold])))
dev.off()

pairs(tremula_cor_parameters_new,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("Gene num",expression(rho),expression(theta[4-fold]),expression(theta[Intergenic])))
pairs(tremuloides_cor_parameters_new,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("Gene num",expression(rho),expression(theta[4-fold]),expression(theta[Intergenic])))
pairs(trichocarpa_cor_parameters_new,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("Gene num",expression(rho),expression(theta[4-fold]),expression(theta[Intergenic])))




png(filename="tremula.1Mb.cor.coding.png",width = 6, height = 6, units = 'in', res=300)
pairs(tremula_cor_parameters_2,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("GC","Coding %",expression(rho),expression(theta[4-fold]),expression(d[4-fold]),expression(theta[Intergenic]),expression(d[Intergenic]),expression(d[0-fold])))
dev.off()

pairs(tremula_cor_parameters_3,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("GC","Gene num",expression(rho),expression(theta[4-fold]),expression(TajD[4-fold]),expression(d[4-fold]),expression(theta[Intergenic]),expression(TajD[Intergenic]),expression(d[Intergenic]),expression(d[0-fold])))
pairs(tremuloides_cor_parameters_3,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("GC","Gene num",expression(rho),expression(theta[4-fold]),expression(TajD[4-fold]),expression(d[4-fold]),expression(theta[Intergenic]),expression(TajD[Intergenic]),expression(d[Intergenic]),expression(d[0-fold])))
pairs(trichocarpa_cor_parameters_3,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("GC","Gene num",expression(rho),expression(theta[4-fold]),expression(TajD[4-fold]),expression(d[4-fold]),expression(theta[Intergenic]),expression(TajD[Intergenic]),expression(d[Intergenic]),expression(d[0-fold])))






###plots for tremuloides
png(filename="tremuloides.1Mb.cor.gene.png",width = 7, height = 7, units = 'in', res=300)
pairs(tremuloides_cor_parameters,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("GC","Gene num",expression(rho),expression(theta[4-fold]),expression(d[4-fold]),expression(theta[Intergenic]),expression(d[Intergenic]),expression(d[0-fold])))
#pairs(tremuloides_cor_parameters_3,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("GC","Gene num",expression(rho),expression(theta[4-fold]),expression(TajD[4-fold]),expression(d[4-fold]),expression(theta[Intergenic]),expression(TajD[Intergenic]),expression(d[Intergenic]),expression(d[0-fold])))
dev.off()

png(filename="tremuloides.1Mb.cor.coding.png",width = 6, height = 6, units = 'in', res=300)
pairs(tremuloides_cor_parameters_2,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("GC","Coding %",expression(rho),expression(theta[4-fold]),expression(d[4-fold]),expression(theta[Intergenic]),expression(d[Intergenic]),expression(d[0-fold])))
dev.off()

###plots for trichocarpa
png(filename="trichocarpa.1Mb.cor.gene.png",width = 7, height = 7, units = 'in', res=300)
pairs(trichocarpa_cor_parameters,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("GC","Gene num",expression(rho),expression(theta[4-fold]),expression(d[4-fold]),expression(theta[Intergenic]),expression(d[Intergenic]),expression(d[0-fold])))
#pairs(trichocarpa_cor_parameters_3,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("GC","Gene num",expression(rho),expression(theta[4-fold]),expression(TajD[4-fold]),expression(d[4-fold]),expression(theta[Intergenic]),expression(TajD[Intergenic]),expression(d[Intergenic]),expression(d[0-fold])))

dev.off()

png(filename="trichocarpa.1Mb.cor.coding.png",width = 6, height = 6, units = 'in', res=300)
pairs(trichocarpa_cor_parameters_2,lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,labels=c("GC","Coding %",expression(rho),expression(theta[4-fold]),expression(d[4-fold]),expression(theta[Intergenic]),expression(d[Intergenic]),expression(d[0-fold])))
dev.off()






