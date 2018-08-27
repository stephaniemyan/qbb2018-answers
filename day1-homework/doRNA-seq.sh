#!/bin/bash

GENOME=../genomes/BDGP6
ANNOTATION=../genomes/BDGP6.Ensembl.81.gtf
FOLDER=$1

for SAMPLE in SRR072893 SRR072903 SRR072905 SRR072915
do
    mkdir $SAMPLE

    # generates quality control report for sample reads
    fastqc $FOLDER/$SAMPLE.fastq

    # maps sample reads to BDGP6, outputs as .sam
    hisat2 -x $GENOME -U $FOLDER/$SAMPLE.fastq -S $SAMPLE.sam
    
    # converts .sam file into a sorted .bam file and .bam index
    samtools sort $SAMPLE.sam -o ${SAMPLE}_sorted.bam
    samtools index ${SAMPLE}_sorted.bam ${SAMPLE}_index.bam

    # quantitates sorted .bam file
    stringtie ${SAMPLE}_sorted.bam -p 5 -e -G $ANNOTATION -o ${SAMPLE}_quantitated.gtf -B

done
