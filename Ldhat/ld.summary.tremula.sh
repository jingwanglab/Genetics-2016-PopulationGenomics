#! /bin/bash -l

# This is a small script using vcftools to remove indels

set +e
set -x

#SBATCH -A b2010014
#SBATCH -p core
#SBATCH -o ld.summary.tremula.out
#SBATCH -e ld.summary.tremula.err
#SBATCH -J ld.summary.tremula.job
#SBATCH -t 1-00:00:00
#SBATCH --mail-user jing.wang@emg.umu.se
#SBATCH --mail-type=ALL

ld_Dir="/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremula/ld/$1"
ld_bed="/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremula/ldhat/bed/$1"
Output="/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremula/ld/summary/$1"

if [ ! -d "$Output" ]; then
mkdir -p $Output
fi

mean="/proj/b2011141/tools/mean"

echo -e "Chr\tWin\tScaffold\tPair_Num\tr2" > $Output/ld.tremula.$1.summary.txt

for file in $ld_Dir/*ld
do
input=${file##*/}
chr_out=${input%%.geno.ld}
number=${chr_out##*.}
chr_n=${chr_out#tremula_}
chr_number=${chr_n%.*}
echo $chr_number
chr="Chr"$chr_number
echo $chr

line=$(cat $file | wc -l)
echo $line
###n, rather than represent the number of SNPs, it represents the pairwise comparison number 
n=$(echo "$line-1" | bc)
echo $n
ld_mean=$($mean col=5 header=1 $file | sed '1d')

scaffold=$(sed '1d' $ld_bed/${chr_out}.bed |cut -f 2 )
echo $scaffold

if [ $1 == "10kb" ] ; then
middle_scaffold=$(echo $scaffold+5001 |bc)
fi
if [ $1 == "100kb" ] ; then
middle_scaffold=$(echo $scaffold+50001 |bc)
fi
if [ $1 == "500kb" ] ; then
middle_scaffold=$(echo $scaffold+250001 |bc)
fi
if [ $1 == "1Mb" ] ; then
middle_scaffold=$(echo $scaffold+500001 |bc)
fi

echo $middle_scaffold
echo -e "$chr\t$number\t$middle_scaffold\t$n\t$ld_mean" >> $Output/ld.tremula.$1.summary.txt
done

sort -k1,1V -k2,2n $Output/ld.tremula.$1.summary.txt > $Output/temp && mv $Output/temp $Output/ld.tremula.$1.summary.txt


