#!/bin/bash
set -e

# just random samples from 1000 genome project pilot data
SAMPLES='NA06984 NA06985 NA06989 NA06994 NA07000 NA11832 NA11840 NA11829 NA11830 NA11831'
URL='ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/data/'
regions_f='regions.txt'
range_start="6000000"
range_end="6100000"

if [ -f ./${regions_f} ]
then
 echo "regions file exist, not overwriting"
else
 for chr in `seq 1 22` X Y; do echo -e "$chr\t$range_start\t$range_end" >>${regions_f}; done
fi

echo "Getting data from $URL"

for sample in $SAMPLES
do
 echo "Getting $sample"
 samtools view -b -h  ${URL}/${sample}/alignment/${sample}.454.MOSAIK.SRP000033.2009_11.bam -L ./$regions_f >${sample}.bam
done

echo "Generation index for bam files"
parallel  samtools index ::: *.bam
