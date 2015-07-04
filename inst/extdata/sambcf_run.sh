#!/bin/bash -l

## Load samtools from module system (delete this line if a module system is not available)
module load samtools
module load bcftools

## Path to refenence genome needs to be set here
REF="./data/tair10.fasta"

## Remove duplicated reads (PCR artifacts)
samtools rmdup -S myfile.fastq.bam myfile.fastq.dedup.bam
samtools sort myfile.fastq.dedup.bam myfile.fastq.dedup.sorted
samtools index myfile.fastq.dedup.sorted.bam 

## Variant calling
samtools mpileup -uf $REF myfile.fastq.dedup.sorted.bam | bcftools view -bvcg -> sambcf.raw.bcf
bcftools view sambcf.raw.bcf | vcfutils.pl varFilter -D100 > sambcf.vcf

