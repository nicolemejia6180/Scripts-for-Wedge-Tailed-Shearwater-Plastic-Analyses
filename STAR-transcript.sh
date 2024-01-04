#!/bin/bash
#SBATCH -p edwards
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --mem 50000
#SBATCH -t 10-00:00
#SBATCH -J STAR-alig
#SBATCH -o STAR-alig%j.out
#SBATCH -e STAR-alig%j.err
#SBATCH --mail-type=ALL        # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ftermignoni@fas.harvard.edu # Email to which notifications will be sent

module load GCC/7.3.0-2.30 OpenMPI/3.1.1 STAR/2.7.1a

#cd /n/holyscratch01/edwards_lab/ftermignoni/CB-project/transcriptome/STAR-align/CYUCA-RNA
#module purge
#module load python
#source activate rsem

for i in *_R1.fastq.gz; do
STAR --twopassMode Basic --quantMode TranscriptomeSAM --readFilesCommand zcat --outSAMtype BAM Unsorted --genomeDir ../IndexCmonedula --readFilesIn $i ${i%_R1.fastq.gz}_R2.fastq.gz --runThreadN 10 --outFileNamePrefix ${i%_R1.fastq.gz};done

