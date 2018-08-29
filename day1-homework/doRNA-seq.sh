#!/bin/bash

# doRNA-seq.sh runs FastQC, HISAT2, SAMtools, and StringTie on a and reference genome of your choice
# You must run this in a directory that includes your .fastq read files
# No other inputs needed in the command line

# UPDATE with filepath for HISAT2 index for your reference genome
GENOME=~/qbb2018-answers/genomes/BDGP6
# UPDATE with filepath for your reference genome
ANNOTATION=~/qbb2018-answers/genomes/BDGP6.Ensembl.81.gtf

# UPDATE with names of your sample files
for SAMPLE in SRR072893 SRR072903 SRR072905 SRR072915

do
	# Make new folder
    mkdir $SAMPLE

    # Generate quality control report for sample reads
	echo 'Running FastQC...'
    fastqc -o $SAMPLE $SAMPLE.fastq

    # Map sample reads to reference genome and output as .sam
	echo 'Running HISAT2...'
    hisat2 -x $GENOME -U $SAMPLE.fastq -S $SAMPLE/$SAMPLE.sam
    
    # Convert .sam file into a sorted .bam file and .bam index
	echo 'Sorting .sam file...'
    samtools sort $SAMPLE/$SAMPLE.sam -o $SAMPLE/${SAMPLE}_sorted.bam
	echo 'Making index from .sam file...'
    samtools index $SAMPLE/${SAMPLE}_sorted.bam $SAMPLE/${SAMPLE}_index.bam

    # Quantitate sorted .bam file
	echo 'Running StringTie...'
    stringtie $SAMPLE/${SAMPLE}_sorted.bam -p 5 -e -G $ANNOTATION -o $SAMPLE/${SAMPLE}_quantitated.gtf -B

done
