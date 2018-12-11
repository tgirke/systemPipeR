#!/bin/bash -l

## Load samtools from module system (delete this line if a module system is not available)
module load samtools/0.1.19

## Path to refenence genome needs to be set here
REF="/rhome/tgirke/Projects/github/systemPipeR_tests/data/tair10.fasta"

## Remove duplicated reads (PCR artifacts)
samtools rmdup -S myfile.fastq.bam myfile.fastq.dedup.bam
samtools sort myfile.fastq.dedup.bam myfile.fastq.dedup.sorted
samtools index myfile.fastq.dedup.sorted.bam 

## Switch of samtools version since some utilities are not yet implemented in new version
module unload samtools/0.1.19
module load samtools/1.2
module load bcftools/1.2

## Variant calling 
## Old version <= 0.1.19
# samtools mpileup -uf $REF myfile.fastq.dedup.sorted.bam | bcftools view -bvcg -> sambcf.raw.bcf
# bcftools view sambcf.raw.bcf | vcfutils.pl varFilter -D100 > sambcf.vcf
## New version
samtools mpileup -uf $REF myfile.fastq.dedup.sorted.bam | bcftools call -cv - > sambcf.vcf
#bcftools filter -sLowQual -g3 -G10 -e'%QUAL<10 || (RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || %MAX(DV)<=3 || %MAX(DV)/%MAX(DP)<=0.3' sambcf.vcf > sambcf.vcf

