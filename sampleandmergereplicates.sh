#!/bin/sh

module load ngs/samtools/1.9
module load ngs/Homer/4.9
module load ngs/UCSCutils
module load ngs/bedtools2

mkdir -p out

SAMPLES=`ls [0-9]*.bed | grep -v "merged" | grep -v "sample"| sed -e 's/[0-9][0-9]_\(.*\).bed/\1/'|sed -e 's/1//g'| sort| uniq `
echo ${SAMPLES}

for SAMPLE in ${SAMPLES}
	do
	A=`ls *_${SAMPLE}*.bed| grep -v "merged" | grep -v "sample"`

	for B in ${A}
		do
		wc -l ${B} >> ${SAMPLE}.txt
	done 

	for B in ${A}
		do
		MIN=`cat ${SAMPLE}.txt | sort -n | awk '{print $1}' | head -n 1` 
		bedtools sample -i ${B} -n ${MIN} > sample.${B}
		done

A=`ls sample*_${SAMPLE}*.bed`

cat ${A} > merged.${SAMPLE}.bed
makeTagDirectory merged.${SAMPLE}.dir ${A} -single -fragLength 150

done


INPUT=`ls merged*.dir -d | grep "neg"`
SAMPLES=`ls merged*.dir -d | grep -v "neg"`

INPUTBASE=`echo ${INPUT} | sed -e 's/.dir//g'`
SAMPLESBASE=`echo ${SAMPLES} | sed -e 's/.dir//g'`

echo $INPUTBASE
echo $SAMPLESBASE


for FILEBASE in ${SAMPLESBASE}
do

	sbatch --export=FILEBASE=$FILEBASE,INPUTBASE=$INPUTBASE homer.sbatch
done

