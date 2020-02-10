#! /bin/bash -l

# This is a small script using vcftools to remove indels

set -e
set -x

#SBATCH -A b2011141
#SBATCH -p core
#SBATCH -o Vcf_ldhat.tremula.anno.out
#SBATCH -e Vcf_ldhat.tremula.anno.err
#SBATCH -J Vcf_ldhat.tremula.anno.job
#SBATCH -t 2-00:00:00
#SBATCH --mail-user jing.wang@emg.umu.se
#SBATCH --mail-type=ALL

set +e

module load bioinfo-tools
#module load vcftools

vcftools="/proj/b2011141/tools/vcftools/bin/vcftools"
bgzip="/proj/b2011141/tools/bgzip"

Inputvcf=$1
VCFDir=`dirname $1`
anno=$2 ###3UTR 5UTR intron zero_fold four_fold intergenic
window=$3 ##100kb 500kb 1Mb
InputBed="/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremula/ldhat/bed/$3"
chr=$4

interval="/proj/b2011141/tools/ldhat/interval"
stat="/proj/b2011141/tools/ldhat/stat"
lkgen="/proj/b2011141/tools/ldhat/lkgen"

lookup="/proj/b2011141/tools/ldhat/lookup/3species/tremula_trichocarpa.lk"

Out=${Inputvcf##*/}
echo $Out

ldhat=$VCFDir/ldhat/$anno/$window/input
ldhat_interval=$VCFDir/ldhat/$anno/$window/interval

if [ ! -d "$ldhat" ]; then
mkdir -p $ldhat
fi

if [ ! -d "$ldhat_interval" ]; then
mkdir -p $ldhat_interval
fi

for file in $InputBed/tremula_$chr.{1..9}*.bed
do
bedfile=${file##*/}
output=${bedfile%.bed}
$vcftools --gzvcf $1 --bed $file --ldhat-geno --chr Chr$chr --out $ldhat/$output
$interval -seq $ldhat/$output.ldhat.sites -loc $ldhat/$output.ldhat.locs -lk $lookup -its 1000000 -bpen 5 -samp 2000 -prefix $ldhat_interval/$output.
$stat -input $ldhat_interval/$output.rates.txt -burn 50 -loc $ldhat/$output.ldhat.locs -prefix $ldhat_interval/$output.
rm $ldhat_interval/$output.bounds.txt $ldhat_interval/$output.new_lk.txt $ldhat_interval/$output.type_table.txt $ldhat_interval/$output.rates.txt
done


