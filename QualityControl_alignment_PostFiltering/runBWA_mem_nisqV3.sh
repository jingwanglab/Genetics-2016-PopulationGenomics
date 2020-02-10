#! /bin/bash -l

set -e

Cores=8

module load bioinfo-tools
module load samtools/0.1.19
module load bwa


ref="/proj/b2011141/nobackup/reference/nisqV3/Ptrichocarpa_v3.0_210.fna"

fileIn=${1%.all-trimmomatic.1.fq.gz}
speciesP=${1##*QCfiltered/}
species=${speciesP%%/*}
sample=${fileIn##*/}
ind=${sample%_*}

BWA_DIR="/proj/b2011141/nobackup/private/new_nobackup_area/bam/trichocarpa/$species"

if [ ! -d "$BWA_DIR" ]; then
mkdir -p $BWA_DIR
fi

 [ -e $ref.sa ] || $bwa index $ref
 [ -e $ref.fai ] || samtools faidx $ref

#$bwa mem -t $Cores -M $ref $fileIn.all-trimmomatic.1.fq.gz $fileIn.all-trimmomatic.2.fq.gz | samtools import $ref.fai - - |samtools sort - $BWA_DIR/$sample.PE
bwa mem -t $Cores -R "@RG\tID:$sample\tPL:illumina\tPU:1\tLB:PE.$lib\tSM:$ind" -M $ref $fileIn.all-trimmomatic.1.fq.gz $fileIn.all-trimmomatic.2.fq.gz | samtools import $ref.fai - - |samtools sort - $BWA_DIR/$sample.PE

#$bwa mem -t $Cores -M $ref $fileIn.all-trimmomatic.se.fq.gz | samtools import $ref.fai - -| samtools sort - $BWA_DIR/$sample.SE
bwa mem -t $Cores -R "@RG\tID:$sample\tPL:illumina\tPU:1\tLB:SE.$lib\tSM:$ind" -M $ref $fileIn.all-trimmomatic.se.fq.gz | samtools import $ref.fai - -| samtools sort - $BWA_DIR/$sample.SE

samtools flagstat $BWA_DIR/$sample.PE.bam > $BWA_DIR/$sample.PE.stats
samtools flagstat $BWA_DIR/$sample.SE.bam > $BWA_DIR/$sample.SE.stats


