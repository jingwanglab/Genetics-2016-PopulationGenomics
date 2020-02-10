#! /usr/bin/Rscript --no-save --no-restore

##install.packages("data.table")
library(data.table)
#library(boot)

##read table
four_fold=fread("/proj/b2010014/GenomePaper/population_genetics/pan_genome/shared/three_species/four_fold/three_species_ref_allele.four_fold.bed")
zero_fold=fread("/proj/b2010014/GenomePaper/population_genetics/pan_genome/shared/three_species/zero_fold/three_species_ref_allele.zero_fold.bed")

n_four_fold=3398334
n_zero_fold=16518838

###----1------real data for SFS
##1.4-fold

four_fold$tremuloides_maf=four_fold$V5*44
four_fold$tremuloides_maf=as.integer(four_fold$tremuloides_maf)

four_fold_count_table=as.numeric(table(four_fold$tremuloides_maf))

four_fold_sfs=rep(0,45)
for (i in 1:23) {four_fold_sfs[i]=four_fold_count_table[i]+four_fold_count_table[46-i]}
four_fold_sfs[23]=four_fold_sfs[23]/2
sum=sum(four_fold_sfs)
four_fold_sfs[1]=n_four_fold-sum+four_fold_sfs[1]


##2.zero_fold
zero_fold$tremuloides_maf=zero_fold$V5*44
zero_fold$tremuloides_maf=as.integer(zero_fold$tremuloides_maf)

zero_fold_count_table=as.numeric(table(zero_fold$tremuloides_maf))

zero_fold_sfs=rep(0,45)
for (i in 1:23) {zero_fold_sfs[i]=zero_fold_count_table[i]+zero_fold_count_table[46-i]}
zero_fold_sfs[23]=zero_fold_sfs[23]/2
sum=sum(zero_fold_sfs)
zero_fold_sfs[1]=n_zero_fold-sum+zero_fold_sfs[1]


###---2-----creating output file for est-dfe
###1.0-fold
tremuloides="/proj/b2010014/GenomePaper/population_genetics/pan_genome/shared/three_species/zero_fold/dfe/tremuloides/real"
cat("1",file=paste(tremuloides,"/tremuloides.zero_fold.real.txt",sep=""))
cat("\n",file=paste(tremuloides,"/tremuloides.zero_fold.real.txt",sep=""),append=T)
cat("44",file=paste(tremuloides,"/tremuloides.zero_fold.real.txt",sep=""),append=T)
cat("\n",file=paste(tremuloides,"/tremuloides.zero_fold.real.txt",sep=""),append=T)
cat(zero_fold_sfs,sep=" ",file=paste(tremuloides,"/tremuloides.zero_fold.real.txt",sep=""),append=T)
cat("\n",file=paste(tremuloides,"/tremuloides.zero_fold.real.txt",sep=""),append=T)
cat(four_fold_sfs,sep=" ",file=paste(tremuloides,"/tremuloides.zero_fold.real.txt",sep=""),append=T)


###----3-----creating the output file for est_alpha_omega for when running DFE
##1.4-fold (neutral)
fixed_four_fold=length(which(four_fold$V8=="fixed"))

##2.zero-fold 
fixed_zero_fold=length(which(zero_fold$V8=="fixed"))

###1.-fold
cat("1",n_zero_fold,fixed_zero_fold,file=paste(tremuloides,"/tremuloides.zero_fold.omega.txt",sep=""),sep=" ")
cat("\n",file=paste(tremuloides,"/tremuloides.zero_fold.omega.txt",sep=""),append=T)
cat("0",n_four_fold,fixed_four_fold,file=paste(tremuloides,"/tremuloides.zero_fold.omega.txt",sep=""),sep=" ",append=T)


