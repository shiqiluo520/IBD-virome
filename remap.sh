#!/bin/bash

prefix=$1
refs=$2
reads1=$3
reads2=$4
fid=$5

bowtie2-build $refs $prefix
bowtie2 -x $prefix -1 $reads1 -2 $reads2 -S ${fid}.sam
samtools view -bS ${fid}.sam > ${fid}_unsorted.bam
samtools sort ${fid}_unsorted.bam -o ${fid}_sorted.bam
samtools index ${fid}_sorted.bam
coverm filter --bam-files ${fid}_sorted.bam --output-bam-files ${fid}.bam -t 1 --min-read-percent-identity 0.95