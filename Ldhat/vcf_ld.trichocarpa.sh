#! /bin/bash -l

#SBATCH -A b2010014
#SBATCH -p core
#SBATCH -o Vcf_ld.trichcoarpa.out
#SBATCH -e Vcf_ld.trichcoarpa.err
#SBATCH -J Vcf_ld.trichcoarpa.job
#SBATCH -t 5-00:00:00
#SBATCH --mail-user jing.wang@emg.umu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
#module load vcftools

vcftools="/proj/b2011141/tools/vcftools_0.1.12/bin/vcftools"

Inputvcf=$1
VCFDir=`dirname $1`
window=$2 ##100kb 500kb 1Mb
InputBed="/proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/trichocarpa/ldhat/bed/$2"
chr=$3

Out=${Inputvcf##*/}
echo $Out

ld=$VCFDir/ld/$2

if [ ! -d "$ld" ]; then
mkdir -p $ld
fi

for file in $InputBed/trichocarpa_$chr.{1..9}*.bed
do
bedfile=${file##*/}
output=${bedfile%.bed}
$vcftools --gzvcf $1 --geno-r2 --ld-window-bp-min 1000 --ld-window-min 10 --chr Chr$chr --bed $file --out $ld/$output
done



