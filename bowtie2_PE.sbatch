#! /bin/bash
#
# bowtie.sbatch
#
#SBATCH -J bowtie_array   # A single job name for the array
#SBATCH -p slim18                    # Partition
#SBATCH -n 12                        # 12 cores
#SBATCH -N 6                         # one node ?required
#SBATCH -t 0-6:00                    # Running time of 2 hours
#SBATCH --mem 20000                  # Memory request
#SBATCH -o out/bowtie_%A_%a.out          # Standard output
#SBATCH -e out/bowtie_%A_%a.err          # Standard error

# grab out filename

FILENAME=`ls *_1.txt.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
FILEBASE=`echo ${FILENAME} | sed -e "s/_1.txt.gz//g"`

module load ngs/bowtie2
module load ngs/samtools
module load ngs/bedtools2
module load ngs/Homer

BOWTIE_INDEX="/work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome"

BOWTIE_OPTS="-p 24 --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant -I 10 -X 130"
bowtie2 $BOWTIE_OPTS -x $BOWTIE_INDEX -1 ${FILEBASE}_1.txt.gz -2 ${FILEBASE}_2.txt.gz  > ${FILEBASE}_sub.sam 2> ${FILEBASE}_sub.stats
samtools view -bS -@ 12 -q 10 ${FILEBASE}_sub.sam | samtools sort -@ 12 - | tee ${FILEBASE}_sub.bam | samtools index - ${FILEBASE}_sub.bam.bai
makeTagDirectory ${FILEBASE}_sub.dir ${FILEBASE}_sub.bam -single

BOWTIE_OPTS="-p 24 --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant -I 130 -X 220"
bowtie2 $BOWTIE_OPTS -x $BOWTIE_INDEX -1 ${FILEBASE}_1.txt.gz -2 ${FILEBASE}_2.txt.gz  > ${FILEBASE}_mono.sam 2> ${FILEBASE}_mono.stats
samtools view -bS -@ 12 -q 10 ${FILEBASE}_mono.sam | samtools sort -@ 12 - | tee ${FILEBASE}_mono.bam | samtools index - ${FILEBASE}_mono.bam.bai
makeTagDirectory ${FILEBASE}_mono.dir ${FILEBASE}_mono.bam -single

BOWTIE_OPTS="-p 24 --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant -I 10 -X 220"
bowtie2 $BOWTIE_OPTS -x $BOWTIE_INDEX -1 ${FILEBASE}_1.txt.gz -2 ${FILEBASE}_2.txt.gz  > ${FILEBASE}_all.sam 2> ${FILEBASE}_all.stats
samtools view -bS -@ 12 -q 10 ${FILEBASE}_all.sam | samtools sort -@ 12 - | tee ${FILEBASE}_all.bam | samtools index - ${FILEBASE}_all.bam.bai
makeTagDirectory ${FILEBASE}_all.dir ${FILEBASE}_all.bam -single

