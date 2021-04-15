#!/bin/sh

#SBATCH -J bamtobed  # A single job name for the array
#SBATCH -p slim18                    # Partition
#SBATCH -n 6                         # 6 cores
#SBATCH -N 1                         # one node ?required
#SBATCH -t 0-4:00                    # Running time of 2 hours
#SBATCH --mem 40000                  # Memory request
#SBATCH -o out/findPeaks_%A.out          # Standard output
#SBATCH -e out/findPeaks_%A.err          # Standard error

module load ngs/samtools/1.9
module load ngs/Homer/4.9
module load ngs/UCSCutils
module load ngs/bedtools2

samtools view -h -q 10 ${SNAME}.bam -b -o ${SNAME}.f.bam
samtools sort -n -m 2G -@ 16 -o ${SNAME}.s.bam ${SNAME}.f.bam

bamToBed -i ${SNAME}.s.bam -bedpe | cut -f 1,2,6 > ${SNAME}.bed 2>/dev/null
cut -f 1,2,6 ${SNAME}.bed > ${SNAME}.s.bed


