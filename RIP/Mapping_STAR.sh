#!/bin/bash
#SBATCH --mem=200g
#SBATCH -t 24:00:00
#SBATCH -n 16
#SBATCH -p general

module load anaconda
source activate Gilbert-work

gDir="/work/users/g/g/ggiri/GENOME/Human/STAR"
gtf="/work/users/g/g/ggiri/GENOME/Human/Homo_sapiens.GRCh38.108.chr.gtf"
fastq="/work/users/g/g/ggiri/PCBP1/RIP_HEK/Data"
outDir="/work/users/g/g/ggiri/PCBP1/RIP_HEK/STAR"
file=$1

/work/users/g/g/ggiri/Tools/STAR/bin/Linux_x86_64_static/STAR \
        --runThreadN 16 \
        --sjdbGTFfile ${gtf} \
        --genomeDir ${gDir} \
        --readFilesIn ${fastq}/${file}.fastq \
        --outFilterType BySJout \
        --outFilterMismatchNmax 5 \
        --genomeLoad NoSharedMemory \
        --outFileNamePrefix ${outDir}/${file} \
        --runRNGseed 777 \
        --outSAMtype BAM Unsorted


