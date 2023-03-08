#!/bin/bash
#SBATCH --mem=120g
#SBATCH -t 06:00:00
#SBATCH -n 8
#SBATCH -p general

module load anaconda
source activate Gilbert-work

index="/work/users/g/g/ggiri/GENOME/Human/RSEM/hg38"
bamdir="/work/users/g/g/ggiri/PCBP1/RIP_C2BBe1/STAR"
outdir="/work/users/g/g/ggiri/PCBP1/RIP_C2BBe1/RSEM"
file=$1

rsem-calculate-expression --num-threads 8 \
        --alignments \
        --seed 77 \
        ${bamdir}/${file}Aligned.toTranscriptome.out.bam \
        ${index} \
        ${outdir}/${file}
