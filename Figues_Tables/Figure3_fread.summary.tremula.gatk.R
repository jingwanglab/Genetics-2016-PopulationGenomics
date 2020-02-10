#! /usr/bin/Rscript --no-save --no-restore

setwd("/gulo/proj_nobackup/b2011141/genomic_selection_paper/tremula/count2/annotations")
##install.packages("data.table")
library(data.table)
#library(boot)

##read table
four_fold=fread("/proj/b2010014/GenomePaper/population_genetics/pan_genome/shared/three_species/four_fold/three_species_ref_allele.four_fold.bed")
zero_fold=fread("/proj/b2010014/GenomePaper/population_genetics/pan_genome/shared/three_species/zero_fold/three_species_ref_allele.zero_fold.bed")
#zero_fold=fread("/proj/b2010014/GenomePaper/population_genetics/pan_genome/shared/three_species/zero_fold/three_species_ref_allele.zero_fold.bed")

n_four_fold=3398334
n_zero_fold=16518838
#n_intergenic=73462757
#n_intron=

###----1------real data for SFS
##1.4-fold

four_fold$tremula_maf=four_fold$V4*48
four_fold$tremula_maf=as.integer(four_fold$tremula_maf)

four_fold_count_table=as.numeric(table(four_fold$tremula_maf))

four_fold_sfs=rep(0,49)
for (i in 1:25) {four_fold_sfs[i]=four_fold_count_table[i]+four_fold_count_table[50-i]}
four_fold_sfs[25]=four_fold_sfs[25]/2
sum=sum(four_fold_sfs)
four_fold_sfs[1]=n_four_fold-sum+four_fold_sfs[1]


##2.zero_fold
zero_fold$tremula_maf=zero_fold$V4*48
zero_fold$tremula_maf=as.integer(zero_fold$tremula_maf)

zero_fold_count_table=as.numeric(table(zero_fold$tremula_maf))

zero_fold_sfs=rep(0,49)
for (i in 1:25) {zero_fold_sfs[i]=zero_fold_count_table[i]+zero_fold_count_table[50-i]}
zero_fold_sfs[25]=zero_fold_sfs[25]/2
sum=sum(zero_fold_sfs)
zero_fold_sfs[1]=n_zero_fold-sum+zero_fold_sfs[1]


###---2-----creating output file for est-dfe
###1.0-fold
tremula="/proj/b2010014/GenomePaper/population_genetics/pan_genome/shared/three_species/zero_fold/dfe/tremula/real"
cat("1",file=paste(tremula,"/tremula.zero_fold.real.txt",sep=""))
cat("\n",file=paste(tremula,"/tremula.zero_fold.real.txt",sep=""),append=T)
cat("48",file=paste(tremula,"/tremula.zero_fold.real.txt",sep=""),append=T)
cat("\n",file=paste(tremula,"/tremula.zero_fold.real.txt",sep=""),append=T)
cat(zero_fold_sfs,sep=" ",file=paste(tremula,"/tremula.zero_fold.real.txt",sep=""),append=T)
cat("\n",file=paste(tremula,"/tremula.zero_fold.real.txt",sep=""),append=T)
cat(four_fold_sfs,sep=" ",file=paste(tremula,"/tremula.zero_fold.real.txt",sep=""),append=T)


###----3-----creating the output file for est_alpha_omega for when running DFE
##1.4-fold (neutral)
fixed_four_fold=length(which(four_fold$V7=="fixed"))

##2.zero-fold 
fixed_zero_fold=length(which(zero_fold$V7=="fixed"))

###1.-fold
cat("1",n_zero_fold,fixed_zero_fold,file=paste(tremula,"/tremula.zero_fold.omega.txt",sep=""),sep=" ")
cat("\n",file=paste(tremula,"/tremula.zero_fold.omega.txt",sep=""),append=T)
cat("0",n_four_fold,fixed_four_fold,file=paste(tremula,"/tremula.zero_fold.omega.txt",sep=""),sep=" ",append=T)


