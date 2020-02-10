#! /bin/bash -l
set -e
set -x

#SBATCH -A b2010014
#SBATCH -p core
#SBATCH -n 2
#SBATCH -o python.out
#SBATCH -e python.err
#SBATCH -J python.job
#SBATCH -t 1-00:00:00
#SBATCH --mail-user jing.wang@emg.umu.se
#SBATCH --mail-type=ALL

python AnnotateRef.py /proj/b2011141/nobackup/reference/nisqV3/Ptrichocarpa_v3.0_210.fa P_trichocarpa.cds.longest.txt > P_trichocarpa.cds.longest.annotation.txt 
