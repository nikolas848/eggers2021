#!/bin/sh

#SBATCH -J findPeaksHomer   # A single job name for the array
#SBATCH -p slim16                    # Partition
#SBATCH -n 6                         # 6 cores
#SBATCH -N 1                         # one node ?required
#SBATCH -t 0-4:00                    # Running time of 2 hours
#SBATCH --mem 40000                  # Memory request
#SBATCH -o out/findPeaks_%A.out          # Standard output
#SBATCH -e out/findPeaks_%A.err          # Standard error

module load ngs/samtools/1.9
module load ngs/Homer/4.9
module load ngs/UCSCutils

#findPeaks ${FILEBASE}.dir  -i ${INPUTBASE}.dir -style histone -F 4 -L 2 -o ${FILEBASE}.h04.l02.txt


findPeaks ${FILEBASE}.dir  -i ${INPUTBASE}.dir -style factor -size 150 -F 4 -L 2 -o ${FILEBASE}.f04.l02.txt
findPeaks ${FILEBASE}.dir  -i ${INPUTBASE}.dir -style factor -size 150 -F 6 -L 2 -o ${FILEBASE}.f06.l02.txt
findPeaks ${FILEBASE}.dir  -i ${INPUTBASE}.dir -style factor -size 150 -F 8 -L 2 -o ${FILEBASE}.f08.l02.txt
for TXT in *.txt
do
		TXTBASE=`echo ${TXT} | sed -e 's/.txt//g'`
		pos2bed.pl - ${TXTBASE}.txt > ${TXTBASE}.bed
done

		makeUCSCfile ${FILEBASE}.dir  -i ${INPUTBASE}.dir  -o ${FILEBASE}.bedgraph
		gunzip ${FILEBASE}.bedgraph.gz
		bedGraphToBigWig ${FILEBASE}.bedgraph /work/project/becbec_003/drosophila_genome/dm6.chromsizes.txt ${FILEBASE}.bw

