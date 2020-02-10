#! /usr/bin/Rscript --no-save --no-restore

#Here, this R script is used to make LD decay plot, which based on the genotype_r2 produced by vcftools

library(RColorBrewer)
library(ggplot2)
library(gridExtra)

colors <- brewer.pal(9,"Set1")[c(5,2,3)]

setwd("/proj/b2010014/GenomePaper/population_genetics/pan_genome/plots/ld")

#tremula
tremula=read.table("/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremula/plink/thin_ld/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.tremula.plink.thin.r2.ld",header=TRUE)
#tremuloides
tremuloides=read.table("/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremuloides/plink/thin_ld/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.tremuloides.plink.thin.r2.ld",header=TRUE)
#trichocarpa
trichocarpa=read.table("/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/trichocarpa/plink/thin_ld/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.trichocarpa.plink.thin.r2.ld",header=TRUE)

tremula_n=48
tremuloides_n=44
trichocarpa_n=48


#LD Function
tremula_Er<-function(C_,d){
length(d)
res<-((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2)/(tremula_n*(2+C_*d)*(11+C_*d))))
return(res)
}

tremuloides_Er<-function(C_,d){
length(d)
res<-((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2)/(tremuloides_n*(2+C_*d)*(11+C_*d))))
return(res)
}

trichocarpa_Er<-function(C_,d){
length(d)
res<-((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2)/(trichocarpa_n*(2+C_*d)*(11+C_*d))))
return(res)
}


#for P.tremula 
tremula_d=(tremula$BP_B-tremula$BP_A)
tremula_d_order=tremula[order(tremula_d),]
tremula_ld=tremula_d_order$R2

tremula_nlm=nls(tremula_ld~tremula_Er(C_,tremula_d[order(tremula_d)]),start=list(C_=0.01))
tremula_C_<-summary(tremula_nlm)$coefficients[1]

#for P.tremuloides
tremuloides_d=(tremuloides$BP_B-tremuloides$BP_A)
tremuloides_d_order=tremuloides[order(tremuloides_d),]
tremuloides_ld=tremuloides_d_order$R2

tremuloides_nlm=nls(tremuloides_ld~tremuloides_Er(C_,tremuloides_d[order(tremuloides_d)]),start=list(C_=0.01))
tremuloides_C_<-summary(tremuloides_nlm)$coefficients[1]

#for P.trichocarpa
trichocarpa_d=(trichocarpa$BP_B-trichocarpa$BP_A)
trichocarpa_d_order=trichocarpa[order(trichocarpa_d),]
trichocarpa_ld=trichocarpa_d_order$R2

trichocarpa_nlm=nls(trichocarpa_ld~trichocarpa_Er(C_,trichocarpa_d[order(trichocarpa_d)]),start=list(C_=0.01))
trichocarpa_C_<-summary(trichocarpa_nlm)$coefficients[1]

tremula_table=data.frame(Distance=tremula_d[order(tremula_d)],r2=tremula_Er(tremula_C_,tremula_d[order(tremula_d)]),species="P.tremula")
tremuloides_table=data.frame(Distance=tremuloides_d[order(tremuloides_d)],r2=tremuloides_Er(tremuloides_C_,tremuloides_d[order(tremuloides_d)]),species="P.tremuloides")
trichocarpa_table=data.frame(Distance=trichocarpa_d[order(trichocarpa_d)],r2=trichocarpa_Er(trichocarpa_C_,trichocarpa_d[order(trichocarpa_d)]),species="P.trichocarpa")

total_table=rbind(tremula_table,tremuloides_table,trichocarpa_table)

#save(total_table,file="total_table.RData")
total_plot=ggplot(total_table,aes(x=Distance,y=r2,colour=species))+geom_line(aes(colour = species),size=1.2)+ggtitle("(a)")+theme(plot.title=element_text(hjust=0))+ylab(expression(r^2))+xlab("Distance(bp)")+theme(legend.position="none")+scale_color_manual(values=c(colors[1],colors[2],colors[3]))

#ggsave("LD.3species.png",width=5,height=5,dpi=300)

###tremula_plot
tremula_d_order$distance=tremula_d_order$BP_B-tremula_d_order$BP_A
tremula_quantile_95=rep(0,50000)
tremula_quantile_05=rep(0,50000)
for(i in 1:50000)
{
tremula_quantile_05[i]=quantile(tremula_d_order[which(tremula_d_order$distance==i),"R2"],0.05)
tremula_quantile_95[i]=quantile(tremula_d_order[which(tremula_d_order$distance==i),"R2"],0.95)
}
tremula_new=as.data.frame(cbind(c(1:50000),unique(tremula_table$r2),tremula_quantile_05,tremula_quantile_95))
names(tremula_new)=c("distance","r2","quantile_05","quantile_95")
tremula_plot=ggplot(tremula_new,aes(x=distance,y=r2))+geom_ribbon(aes(ymin=quantile_05,ymax=quantile_95,x=distance),fill=colors[1],alpha=0.2)+geom_line(colour=colors[1])+ggtitle("(b)")+theme(plot.title=element_text(hjust=0))+ylab(expression(r^2))+xlab("Distance(bp)")

###tremuloides_plot
tremuloides_d_order$distance=tremuloides_d_order$BP_B-tremuloides_d_order$BP_A
tremuloides_quantile_95=rep(0,50000)
tremuloides_quantile_05=rep(0,50000)
for(i in 1:50000)
{
tremuloides_quantile_05[i]=quantile(tremuloides_d_order[which(tremuloides_d_order$distance==i),"R2"],0.05)
tremuloides_quantile_95[i]=quantile(tremuloides_d_order[which(tremuloides_d_order$distance==i),"R2"],0.95)
}
tremuloides_new=as.data.frame(cbind(c(1:50000),unique(tremuloides_table$r2),tremuloides_quantile_05,tremuloides_quantile_95))
names(tremuloides_new)=c("distance","r2","quantile_05","quantile_95")
tremuloides_plot=ggplot(tremuloides_new,aes(x=distance,y=r2))+geom_ribbon(aes(ymin=quantile_05,ymax=quantile_95,x=distance),fill=colors[2],alpha=0.2)+geom_line(colour=colors[2])+ggtitle("(c)")+theme(plot.title=element_text(hjust=0))+ylab(expression(r^2))+xlab("Distance(bp)")

###trichocarpa_plot
trichocarpa_d_order$distance=trichocarpa_d_order$BP_B-trichocarpa_d_order$BP_A
trichocarpa_quantile_95=rep(0,50000)
trichocarpa_quantile_05=rep(0,50000)
for(i in 1:50000)
{
trichocarpa_quantile_05[i]=quantile(trichocarpa_d_order[which(trichocarpa_d_order$distance==i),"R2"],0.05)
trichocarpa_quantile_95[i]=quantile(trichocarpa_d_order[which(trichocarpa_d_order$distance==i),"R2"],0.95)
}
trichocarpa_new=as.data.frame(cbind(c(1:50000),unique(trichocarpa_table$r2),trichocarpa_quantile_05,trichocarpa_quantile_95))
names(trichocarpa_new)=c("distance","r2","quantile_05","quantile_95")
trichocarpa_plot=ggplot(trichocarpa_new,aes(x=distance,y=r2))+geom_ribbon(aes(ymin=quantile_05,ymax=quantile_95,x=distance),fill=colors[3],alpha=0.2)+geom_line(colour=colors[3])+ggtitle("(d)")+theme(plot.title=element_text(hjust=0))+ylab(expression(r^2))+xlab("Distance(bp)")


grid.arrange(total_plot,tremula_plot,tremuloides_plot,trichocarpa_plot,nrow=2)
#g=arrangeGrob(total_plot,tremula_plot,tremuloides_plot,trichocarpa_plot,nrow=2)
g=arrangeGrob(tremula_plot,tremuloides_plot,trichocarpa_plot,total_plot,nrow=2)
print(g)
#ggsave(paste("3species.thetas_tajD_rho.",window,".png",sep=""),width=7,height=5.5,dpi=300,g)
ggsave(filename="3species.ld.quantile90.png",width=8,height=8,dpi=300)



