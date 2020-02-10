#! /usr/bin/Rscript --no-save --no-restore

##install.packages("data.table")
library(data.table)
#library(boot)

###read the argument of the number of bootstrap values
args=(commandArgs(TRUE))
bootstrap_n=args[1]

##read table
four_fold=fread("/proj/b2010014/GenomePaper/population_genetics/pan_genome/shared/three_species/four_fold/three_species_ref_allele.four_fold.bed")
zero_fold=fread("/proj/b2010014/GenomePaper/population_genetics/pan_genome/shared/three_species/zero_fold/three_species_ref_allele.zero_fold.bed")

n_four_fold=3398334
n_zero_fold=16518838

f=function(count,number){  ###count is the table and n is the number of sites for either 4-fold or 0-fold sites
        n=c(1:nrow(count))  ###number of rows of count
        n_boot=sample(n,size=length(n),replace=T)
        count_sample=count[n_boot,]
        maf=as.integer(count_sample$V5*44)
	sfs_table=as.numeric(table(maf))
        sfs=rep(0,45)
        for (i in 1:23){sfs[i]=sfs_table[i]+sfs_table[46-i]}
	sfs[23]=sfs[23]/2
	sum=sum(sfs)
	sfs[1]=number-sum+sfs[1]
        fixed_table=length(which(count_sample$V8=="fixed"))
        return(list(sfs,fixed_table))
}


###----1------real data for SFS
##1.4-fold

four_fold_sfs=f(four_fold,n_four_fold)[[1]]

##2.zero_fold

zero_fold_sfs=f(zero_fold,n_zero_fold)[[1]]

###---2-----creating output file for est-dfe
###1.0-fold
Outfile=paste("/proj/b2010014/GenomePaper/population_genetics/pan_genome/shared/three_species/zero_fold/dfe/tremuloides/bootstrap/bootstrap",bootstrap_n,"/tremuloides.zero_fold.bootstrap",bootstrap_n,".txt",sep="")
cat("1",file=Outfile)
cat("\n",file=Outfile,append=T)
cat("44",file=Outfile,append=T)
cat("\n",file=Outfile,append=T)
cat(zero_fold_sfs,sep=" ",file=Outfile,append=T)
cat("\n",file=Outfile,append=T)
cat(four_fold_sfs,sep=" ",file=Outfile,append=T)


###----3-----creating the output file for est_alpha_omega for when running DFE
##1.4-fold (neutral)
fixed_four_fold=f(four_fold,n_four_fold)[[2]]

##2.zero-fold 
fixed_zero_fold=f(zero_fold,n_zero_fold)[[2]]


Outfile_omega=paste("/proj/b2010014/GenomePaper/population_genetics/pan_genome/shared/three_species/zero_fold/dfe/tremuloides/bootstrap/bootstrap",bootstrap_n,"/tremuloides.zero_fold.bootstrap",bootstrap_n,".omega.txt",sep="")

###1.-fold
cat("1",n_zero_fold,fixed_zero_fold,file=Outfile_omega,sep=" ")
cat("\n",file=Outfile_omega,append=T)
cat("0",n_four_fold,fixed_four_fold,file=Outfile_omega,sep=" ",append=T)


