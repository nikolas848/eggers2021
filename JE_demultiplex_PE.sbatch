#!/bin/sh

# bowtie.sbatch

#SBATCH -J JE_demulti  # A single job name for the array
#SBATCH -p slim18                    # Partition
#SBATCH -n 8                      # x cores
#SBATCH -N 1                         # one node
#SBATCH -t 0-6:00                    # Running time of x hours
#SBATCH --mem 20000                  # Memory request
#SBATCH -o JE_%A_%a.out          # Standard output
#SBATCH -e JE_%A_%a.err          # Standard error

for f in *.fastq.gz.fastqsanger.gz; do 
    mv -- "$f" "${f%.fastq.gz.fastqsanger.gz}.fastq.gz"
done


FILEBASE=`ls *R1.fastq.gz | sed -e "s/_R1.fastq.gz//g"`

F1=${FILEBASE}_R1.fastq.gz
F2=${FILEBASE}_R4.fastq.gz
BF=${FILEBASE}_barcodes.txt
I1=${FILEBASE}_R2.fastq.gz
M=${FILEBASE}_JEstats.txt

java -jar  /opt/software/ngs/jemultiplexer/jemultiplexer_1.0.6_bundle.jar F1=$F1 F2=$F2 BF=$BF I1=$I1 M=$M FORCE=TRUE

## rename output
for f in *_[G,A,T,C]*[G,A,T,C]_1.txt.gz; do 
    mv -- "$f" "${f%_[G,A,T,C]*[G,A,T,C]_1.txt.gz}_1.txt.gz"
done

for f in *_[G,A,T,C]*[G,A,T,C]_2.txt.gz; do 
    mv -- "$f" "${f%_[G,A,T,C]*[G,A,T,C]_2.txt.gz}_2.txt.gz"
done

bash bowtie2_run_PE.sh
