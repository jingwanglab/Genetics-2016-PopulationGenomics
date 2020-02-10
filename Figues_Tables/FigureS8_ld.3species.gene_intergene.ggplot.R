#! /usr/bin/Rscript --no-save --no-restore

#Here, this R script is used to make LD decay plot, which based on the genotype_r2 produced by vcftools

library(gridExtra)
library(RColorBrewer)
library(ggplot2)
colors <- brewer.pal(12,"Paired")[c(7,8,1,2,3,4)]

setwd("/proj/b2010014/GenomePaper/population_genetics/pan_genome/plots/ld")

#tremula
tremula_gene=read.table("/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremula/plink_gene_intergene/gene/thin_ld/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.tremula.gene.plink.thin.r2.ld",header=T)
tremula_intergene=read.table("/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremula/plink_gene_intergene/intergene/thin_ld/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.tremula.intergene.plink.thin.r2.ld",header=T)
#tremuloides
tremuloides_gene=read.table("/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremuloides/plink_gene_intergene/gene/thin_ld/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.tremuloides.gene.plink.thin.r2.ld",header=T)
tremuloides_intergene=read.table("/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremuloides/plink_gene_intergene/intergene/thin_ld/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.tremuloides.intergene.plink.thin.r2.ld",header=T)
#trichocarpa
trichocarpa_gene=read.table("/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/trichocarpa/plink_gene_intergene/gene/thin_ld/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.trichocarpa.gene.plink.thin.r2.ld",header=T)
trichocarpa_intergene=read.table("/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/trichocarpa/plink_gene_intergene/intergene/thin_ld/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.trichocarpa.intergene.plink.thin.r2.ld",header=T)

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


#for P.tremula_gene 
tremula_gene_d=(tremula_gene$BP_B-tremula_gene$BP_A)
tremula_gene_d_order=tremula_gene[order(tremula_gene_d),]
tremula_gene_ld=tremula_gene_d_order$R2

tremula_gene_nlm=nls(tremula_gene_ld~tremula_Er(C_,tremula_gene_d[order(tremula_gene_d)]),start=list(C_=0.01))
tremula_gene_C_<-summary(tremula_gene_nlm)$coefficients[1]

#for P.tremula_intergene
tremula_intergene_d=(tremula_intergene$BP_B-tremula_intergene$BP_A)
tremula_intergene_d_order=tremula_intergene[order(tremula_intergene_d),]
tremula_intergene_ld=tremula_intergene_d_order$R2

tremula_intergene_nlm=nls(tremula_intergene_ld~tremula_Er(C_,tremula_intergene_d[order(tremula_intergene_d)]),start=list(C_=0.01))
tremula_intergene_C_<-summary(tremula_intergene_nlm)$coefficients[1]

#for P.tremuloides_gene
tremuloides_gene_d=(tremuloides_gene$BP_B-tremuloides_gene$BP_A)
tremuloides_gene_d_order=tremuloides_gene[order(tremuloides_gene_d),]
tremuloides_gene_ld=tremuloides_gene_d_order$R2

tremuloides_gene_nlm=nls(tremuloides_gene_ld~tremuloides_Er(C_,tremuloides_gene_d[order(tremuloides_gene_d)]),start=list(C_=0.01))
tremuloides_gene_C_<-summary(tremuloides_gene_nlm)$coefficients[1]

#for P.tremuloides_intergene
tremuloides_intergene_d=(tremuloides_intergene$BP_B-tremuloides_intergene$BP_A)
tremuloides_intergene_d_order=tremuloides_intergene[order(tremuloides_intergene_d),]
tremuloides_intergene_ld=tremuloides_intergene_d_order$R2

tremuloides_intergene_nlm=nls(tremuloides_intergene_ld~tremuloides_Er(C_,tremuloides_intergene_d[order(tremuloides_intergene_d)]),start=list(C_=0.01))
tremuloides_intergene_C_<-summary(tremuloides_intergene_nlm)$coefficients[1]

#for P.trichocarpa_gene
trichocarpa_gene_d=(trichocarpa_gene$BP_B-trichocarpa_gene$BP_A)
trichocarpa_gene_d_order=trichocarpa_gene[order(trichocarpa_gene_d),]
trichocarpa_gene_ld=trichocarpa_gene_d_order$R2

trichocarpa_gene_nlm=nls(trichocarpa_gene_ld~trichocarpa_Er(C_,trichocarpa_gene_d[order(trichocarpa_gene_d)]),start=list(C_=0.01))
trichocarpa_gene_C_<-summary(trichocarpa_gene_nlm)$coefficients[1]

#for P.trichocarpa_intergene
trichocarpa_intergene_d=(trichocarpa_intergene$BP_B-trichocarpa_intergene$BP_A)
trichocarpa_intergene_d_order=trichocarpa_intergene[order(trichocarpa_intergene_d),]
trichocarpa_intergene_ld=trichocarpa_intergene_d_order$R2

trichocarpa_intergene_nlm=nls(trichocarpa_intergene_ld~trichocarpa_Er(C_,trichocarpa_intergene_d[order(trichocarpa_intergene_d)]),start=list(C_=0.01))
trichocarpa_intergene_C_<-summary(trichocarpa_intergene_nlm)$coefficients[1]

tremula_gene_table=data.frame(Distance=tremula_gene_d[order(tremula_gene_d)],r2=tremula_Er(tremula_gene_C_,tremula_gene_d[order(tremula_gene_d)]),species="P.tremula",anno="Gene")
tremula_intergene_table=data.frame(Distance=tremula_intergene_d[order(tremula_intergene_d)],r2=tremula_Er(tremula_intergene_C_,tremula_intergene_d[order(tremula_intergene_d)]),species="P.tremula",anno="Intergene")
tremuloides_gene_table=data.frame(Distance=tremuloides_gene_d[order(tremuloides_gene_d)],r2=tremuloides_Er(tremuloides_gene_C_,tremuloides_gene_d[order(tremuloides_gene_d)]),species="P.tremuloides",anno="Gene")
tremuloides_intergene_table=data.frame(Distance=tremuloides_intergene_d[order(tremuloides_intergene_d)],r2=tremuloides_Er(tremuloides_intergene_C_,tremuloides_intergene_d[order(tremuloides_intergene_d)]),species="P.tremuloides",anno="Intergene")
trichocarpa_gene_table=data.frame(Distance=trichocarpa_gene_d[order(trichocarpa_gene_d)],r2=trichocarpa_Er(trichocarpa_gene_C_,trichocarpa_gene_d[order(trichocarpa_gene_d)]),species="P.trichocarpa",anno="Gene")
trichocarpa_intergene_table=data.frame(Distance=trichocarpa_intergene_d[order(trichocarpa_intergene_d)],r2=trichocarpa_Er(trichocarpa_intergene_C_,trichocarpa_intergene_d[order(trichocarpa_intergene_d)]),species="P.trichocarpa",anno="Intergene")

total_table=rbind(tremula_gene_table,tremula_intergene_table,tremuloides_gene_table,tremuloides_intergene_table,trichocarpa_gene_table,trichocarpa_intergene_table)

#save(total_table,file="total_table.gene_intergene.RData")
tremula_p=ggplot(total_table[which(total_table$species=="P.tremula"),],aes(x=Distance,y=r2,colour=anno))+geom_line(aes(colour=anno),size=1.3)+scale_color_manual(values=c(colors[1],colors[2]))+xlab("Distance(bp)")+ylab(expression(r^2))+theme(legend.title=element_blank())+ggtitle("P.tremula")
tremula_p
ggsave("LD.tremula.gene_intergene.png",width=5,height=5,dpi=300)
tremuloides_p=ggplot(total_table[which(total_table$species=="P.tremuloides"),],aes(x=Distance,y=r2,colour=anno))+geom_line(aes(colour=anno),size=1.3)+scale_color_manual(values=c(colors[3],colors[4]))+xlab("Distance(bp)")+ylab(expression(r^2))+theme(legend.title=element_blank())+ggtitle("P.tremuloides")
tremuloides_p
ggsave("LD.tremuloides.gene_intergene.png",width=5,height=5,dpi=300)
trichocarpa_p=ggplot(total_table[which(total_table$species=="P.trichocarpa"),],aes(x=Distance,y=r2,colour=anno))+geom_line(aes(colour=anno),size=1.3)+scale_color_manual(values=c(colors[5],colors[6]))+xlab("Distance(bp)")+ylab(expression(r^2))+theme(legend.title=element_blank())+ggtitle("P.trichocarpa")
trichocarpa_p
ggsave("LD.trichocarpa.gene_intergene.png",width=5,height=5,dpi=300)




