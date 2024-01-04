#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --mem 64000
#SBATCH -p edwards
#SBATCH -o rsem_%A.out
#SBATCH -e rsem_%A.err
#SBATCH -J rsem
#SBATCH -t 1-00:00

module purge
module load python
source activate rsem

#First prepare the reference gene files for the rsem expression calculations
rsem-prepare-reference --gtf ./ref/Calob.gtf \
                        -p 8 \
                        ./ref/Calob.fna \
                        ./ref/IndexCalob

# $1 == STAR alignments bam file
# $2 == name of STAR index
# #3 == sample name

#for i in *bam;do rsem-calculate-expression --star --num-threads ${SLURM_JOB_CPUS_PER_NODE} --alignments $i ../../IndexPmo $i;done
#rsem-calculate-expression --num-threads ${SLURM_JOB_CPUS_PER_NODE} --alignments --paired-end 305gAligned.out.bam  ../../IndexPmo 305g-sample
for i in *bam;do rsem-calculate-expression --num-threads 16 --alignments --paired-end $i ./ref/IndexCalob SAMPLE-$i;done

#prepare the reference gene files for the rsem expression calculations
rsem-prepare-reference --gtf ./ref/Calob.gtf \
                        -p 8 \
                        ./ref/Calob.fna \
                        ./ref/IndexCalob
