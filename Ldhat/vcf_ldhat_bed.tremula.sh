#! /bin/bash -l

set +e
set -x

#SBATCH -A b2010014
#SBATCH -p core
#SBATCH -o Vcf_ld.out
#SBATCH -e Vcf_ld.err
#SBATCH -J Vcf_ld.job
#SBATCH -t 5-00:00:00
#SBATCH --mail-user jing.wang@emg.umu.se
#SBATCH --mail-type=ALL


module load bioinfo-tools
#module load vcftools


##sbatch vcf_ldhat_bed.tremula.sh /proj/b2010014/GenomePaper/population_genetics/GATK/HC/snp_filter/tremula/3species.snp.rm_indel.biallelic.DP5.GQ10.bed.tremula.recode.vcf.gz 1kb
vcftools="/proj/b2011141/tools/vcftools_0.1.12/bin/vcftools"

Inputvcf=$1
VCFDir=`dirname $1`

Out=${Inputvcf##*/}
echo $Out

bed_win=$2
snp_ldhat=$VCFDir/ldhat/bed/$bed_win

if [ ! -d "$snp_ldhat" ]; then
mkdir -p $snp_ldhat
fi

Out_snp_ld=${Out%.recode.vcf.gz}

if [ $bed_win == "1kb" ] ; then
for file in {01..19}
do
pos_dir="/proj/b2010014/GenomePaper/population_genetics/pan_genome/ANGSD/tremula/tremula_$file"
awk '$14>=100' $pos_dir/tremula_$file.thetas1kbwindow.gz.pestPG | cut -f1 | cut -d ")" -f 3 | sed 's/(//g' | sed 's/,/\t/g' | awk '{print ($1-1)"\t"($2-1)}' | awk '{print "Chr'$file'",$0}' OFS="\t" > $snp_ldhat/tremula_$file.pos.bed
done
fi

if [ $bed_win == "5kb" ] ; then
for file in {01..19}
do
pos_dir="/proj/b2010014/GenomePaper/population_genetics/pan_genome/ANGSD/tremula/tremula_$file"
awk '$14>=500' $pos_dir/tremula_$file.thetas5kbwindow.gz.pestPG | cut -f1 | cut -d ")" -f 3 | sed 's/(//g' | sed 's/,/\t/g' | awk '{print ($1-1)"\t"($2-1)}' | awk '{print "Chr'$file'",$0}' OFS="\t" > $snp_ldhat/tremula_$file.pos.bed
done
fi

if [ $bed_win == "10kb" ] ; then
for file in {01..19}
do
pos_dir="/proj/b2010014/GenomePaper/population_genetics/pan_genome/ANGSD/tremula/tremula_$file"
awk '$14>=1000' $pos_dir/tremula_$file.thetas10kbwindow.gz.pestPG | cut -f1 | cut -d ")" -f 3 | sed 's/(//g' | sed 's/,/\t/g' | awk '{print ($1-1)"\t"($2-1)}' | awk '{print "Chr'$file'",$0}' OFS="\t" > $snp_ldhat/tremula_$file.pos.bed
done
fi

if [ $bed_win == "50kb" ] ; then
for file in {01..19}
do
pos_dir="/proj/b2010014/GenomePaper/population_genetics/pan_genome/ANGSD/tremula/tremula_$file"
awk '$14>=5000' $pos_dir/tremula_$file.thetas50kbwindow.gz.pestPG | cut -f1 | cut -d ")" -f 3 | sed 's/(//g' | sed 's/,/\t/g' | awk '{print ($1-1)"\t"($2-1)}' | awk '{print "Chr'$file'",$0}' OFS="\t" > $snp_ldhat/tremula_$file.pos.bed
done
fi

for file in {01..19}
do
nrow=$(cat $snp_ldhat/tremula_$file.pos.bed | wc -l)

for row in $(eval echo "{1..$nrow}")
do
echo $row
awk 'NR=='$row'' $snp_ldhat/tremula_$file.pos.bed > $snp_ldhat/tremula_$file.$row.bed
echo -e "Chrom\tChrom_Start\tChrom_End" | cat - $snp_ldhat/tremula_$file.$row.bed > $snp_ldhat/out && mv $snp_ldhat/out $snp_ldhat/tremula_$file.$row.bed 
done
done


