library(pastecs)
library(RColorBrewer)
colors <- brewer.pal(9,"Set1")[c(5,2,3)]
setwd("~/Dropbox/paper1_genomic_summary/data/100kb")

total_100kb=read.table("100kb.summary.txt",header=T)
total_100kb=total_100kb[which(total_100kb$Thetas_tremula_numsites>10000),]
total_100kb$tremula_ldhat_new=total_100kb$Ldhat_tremula/1000
total_100kb$tremuloides_ldhat_new=total_100kb$Ldhat_tremuloides/1000
total_100kb$trichocarpa_ldhat_new=total_100kb$Ldhat_trichocarpa/1000
#total_100kb$tremula_alpha=(total_100kb$tremula_trichocarpa_zero_fold_dxy/total_100kb$Thetas_tremula_zero_fold_tP)/(total_100kb$tremula_trichocarpa_four_fold_dxy/total_100kb$Thetas_tremula_four_fold_tP)

total_100kb[which(total_100kb$Ldhat_tremula_num<50),]$tremula_ldhat_new="NA"
total_100kb[which(total_100kb$Ldhat_tremuloides_num<50),]$tremuloides_ldhat_new="NA"
total_100kb[which(total_100kb$Ldhat_trichocarpa_num<50),]$trichocarpa_ldhat_new="NA"

###rho/theta

#install.packages("sm")
library(sm)
tremula=as.numeric(as.character(total_100kb$tremula_ldhat_new))/total_100kb$Thetas_tremula_tW
tremula_table=as.data.frame(cbind(tremula,"P.tremula"))
names(tremula_table)=c("rho_tW","species")
tremuloides=as.numeric(as.character(total_100kb$tremuloides_ldhat_new))/total_100kb$Thetas_tremuloides_tW
tremuloides_table=as.data.frame(cbind(tremuloides,"P.tremuloides"))
names(tremuloides_table)=c("rho_tW","species")
trichocarpa=as.numeric(as.character(total_100kb$trichocarpa_ldhat_new))/total_100kb$Thetas_trichocarpa_tW
trichocarpa_table=as.data.frame(cbind(trichocarpa,"P.trichocarpa"))
names(trichocarpa_table)=c("rho_tW","species")
total_table=rbind(tremula_table,tremuloides_table,trichocarpa_table)
total_table_new=total_table[which(total_table$rho_tW!="NA"),]
total_table_new$rho_tW=as.numeric(as.character(total_table_new$rho_tW))
species_f=factor(total_table$species)
with(total_table_new,sm.density.compare(rho_tW,species,col=c(colors[1],colors[2],colors[3]),panel=T,lwd=c(2,2,2),lty=c(1,1,1),xlab=expression(rho/theta[W]),xlim=c(0,2)))





plot.multi.dens <- function(s)
{
  junk.x = NULL
  junk.y = NULL
  for(i in 1:length(s))
  {
    junk.x = c(junk.x, density(s[[i]])$x)
    junk.y = c(junk.y, density(s[[i]])$y)
  }
  xr <- range(junk.x)
  yr <- range(junk.y)
  plot(density(s[[1]]), xlim = xr, ylim = yr, main = "",xlab=expression(rho/theta[W]))
  for(i in 1:length(s))
  {
    lines(density(s[[i]]), xlim = xr, ylim = yr, col = colors[i],lwd=2)
  }
}
#usage:

# the input of the following function MUST be a numeric list
png(filename="3species.rho.tW.png",width = 6, height = 6, units = 'in', res=300)
plot.multi.dens(list(tremula[which(tremula!="NA")],tremuloides[which(tremuloides!="NA")],trichocarpa[which(trichocarpa!="NA")]))
#install.packages("Hmisc")
legend("topright",c("P.tremula","P.tremuloides","P.trichocarpa"),lty=c(1,1,1),lwd=c(2,2,2),col=c(colors[1],colors[2],colors[3]))
dev.off()

mean(tremula,na.rm=T)
mean(tremuloides,na.rm=T)
mean(trichocarpa,na.rm=T)



