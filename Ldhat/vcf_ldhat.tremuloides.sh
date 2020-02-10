#! /bin/bash -l

# This is a small script using vcftools to remove indels

set -e
set -x

#SBATCH -A b2010014
#SBATCH -p core
#SBATCH -o Vcf_ldhat.tremuloides.out
#SBATCH -e Vcf_ldhat.tremuloides.err
#SBATCH -J Vcf_ldhat.tremuloides.job
#SBATCH -t 5-00:00:00
#SBATCH --mail-user jing.wang@emg.umu.se
#SBATCH --mail-type=ALL

set +e

module load bioinfo-tools
#module load vcftools

vcftools="/proj/b2011141/tools/vcftools/bin/vcftools"
bgzip="/proj/b2011141/tools/bgzip"

Inputvcf=$1
VCFDir=`dirname $1`
window=$2 ##100kb 500kb 1Mb
InputBed="/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremuloides/ldhat/bed/$2"
chr=$3

interval="/proj/b2011141/tools/ldhat/interval"
stat="/proj/b2011141/tools/ldhat/stat"
lkgen="/proj/b2011141/tools/ldhat/lkgen"

lookup="/proj/b2011141/tools/ldhat/lookup/3species/tremuloides.lk"

Out=${Inputvcf##*/}
echo $Out

ldhat=$VCFDir/ldhat/input/$2
ldhat_interval=$VCFDir/ldhat/interval/$2

if [ ! -d "$ldhat" ]; then
mkdir -p $ldhat
fi

if [ ! -d "$ldhat_interval" ]; then
mkdir -p $ldhat_interval
fi

for file in $InputBed/tremuloides_$chr.{1..9}*.bed
do
bedfile=${file##*/}
output=${bedfile%.bed}
$vcftools --gzvcf $1 --bed $file --ldhat-geno --chr Chr$chr --out $ldhat/$output
$interval -seq $ldhat/$output.ldhat.sites -loc $ldhat/$output.ldhat.locs -lk $lookup -its 1000000 -bpen 5 -samp 2000 -prefix $ldhat_interval/$output.
$stat -input $ldhat_interval/$output.rates.txt -burn 50 -loc $ldhat/$output.ldhat.locs -prefix $ldhat_interval/$output.
rm $ldhat_interval/$output.bounds.txt $ldhat_interval/$output.new_lk.txt $ldhat_interval/$output.type_table.txt
done


