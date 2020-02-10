#! /usr/bin/Rscript --no-save --no-restore

#Here, this R script is used to make LD decay plot, which based on the genotype_r2 produced by vcftools

library(RColorBrewer)
library(ggplot2)
colors <- brewer.pal(9,"Paired")[c(2,1,3)]

setwd("/proj/b2010014/GenomePaper/population_genetics/pan_genome/plots/ld")

#tremuloides
tremuloides=read.table("/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremuloides/plink/thin_ld/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.tremuloides.plink.thin.r2.ld",header=TRUE)
#alb
alb=read.table("/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremuloides/alb/plink/thin_ld/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.tremuloides.recode.alb.plink.thin.r2.ld",header=TRUE)
#wil
wil=read.table("/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremuloides/wil/plink/thin_ld/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.tremuloides.recode.wil.plink.thin.r2.ld",header=TRUE)

tremuloides_n=44
alb_n=24
wil_n=20

#LD Function
tremuloides_Er<-function(C_,d){
length(d)
res<-((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2)/(tremuloides_n*(2+C_*d)*(11+C_*d))))
return(res)
}

alb_Er<-function(C_,d){
length(d)
res<-((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2)/(alb_n*(2+C_*d)*(11+C_*d))))
return(res)
}

wil_Er<-function(C_,d){
length(d)
res<-((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2)/(wil_n*(2+C_*d)*(11+C_*d))))
return(res)
}

#for P.tremuloides
tremuloides_d=(tremuloides$BP_B-tremuloides$BP_A)
tremuloides_d_order=tremuloides[order(tremuloides_d),]
tremuloides_ld=tremuloides_d_order$R2

tremuloides_nlm=nls(tremuloides_ld~tremuloides_Er(C_,tremuloides_d[order(tremuloides_d)]),start=list(C_=0.01))
tremuloides_C_<-summary(tremuloides_nlm)$coefficients[1]

#for alb
alb_d=(alb$BP_B-alb$BP_A)
alb_d_order=alb[order(alb_d),]
alb_ld=alb_d_order$R2

alb_nlm=nls(alb_ld~alb_Er(C_,alb_d[order(alb_d)]),start=list(C_=0.01))
alb_C_<-summary(alb_nlm)$coefficients[1]

#for wil
wil_d=(wil$BP_B-wil$BP_A)
wil_d_order=wil[order(wil_d),]
wil_ld=wil_d_order$R2

wil_nlm=nls(wil_ld~wil_Er(C_,wil_d[order(wil_d)]),start=list(C_=0.01))
wil_C_<-summary(wil_nlm)$coefficients[1]

tremuloides_table=data.frame(Distance=tremuloides_d[order(tremuloides_d)],r2=tremuloides_Er(tremuloides_C_,tremuloides_d[order(tremuloides_d)]),species="P.tremuloides")
alb_table=data.frame(Distance=alb_d[order(alb_d)],r2=alb_Er(alb_C_,alb_d[order(alb_d)]),species="Alberta")
wil_table=data.frame(Distance=wil_d[order(wil_d)],r2=wil_Er(wil_C_,wil_d[order(wil_d)]),species="Wisconsin")

total_table=rbind(tremuloides_table,alb_table,wil_table)

save(total_table,file="alb_wil_table.RData")
d=ggplot(total_table,aes(x=Distance,y=r2,colour=species))+geom_line(aes(colour = species),size=1.5)
d+scale_color_manual(values=c(colors[1],colors[2],colors[3]))+xlab("Distance(bp)")+ylab(expression(r^2))+theme(legend.title=element_blank())

ggsave("LD.alb_wil.png",width=5,height=5,dpi=300)




