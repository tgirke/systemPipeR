#!/bin/bash -l

######################
## GATK SNP Calling ##
######################
## Please adjust the enivronment setting of this script as indicated below
## Path: run one level above data and results directories
## Note, some of the settings are adopted from this page: http://seqanswers.com/wiki/How-to/exome_analysis

## Load GATK from module system (delete or comment out these rows if a module system is not available)
module load GATK/2.4-3-g2a7af43
module load picard/1.124

## The path to your refenence genome needs to be set here
REF="./data/tair10.fasta"

## Create dictionary for FASTA reference
# java -jar picard/CreateSequenceDictionary.jar R=$REF O=./data/tair10chr.dict

## Check functionality of GATK and inputs
# java -jar GenomeAnalysisTK.jar -T CountReads -R ./data/tair10chr.fasta -I results/SRR064154.fastq.bam

## Marking PCR duplicates (if '$PICARD MarkDuplicates' fails use 'picard/MarkDuplicates.jar')
java -Xmx4g -Djava.io.tmpdir=./tmp \
-jar $PICARD MarkDuplicates \
INPUT=myfile.fastq.bam \
OUTPUT=myfile.fastq.marked.bam \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
METRICS_FILE=metrics.txt \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT

## Local realignment around indels (if '$GATK' fails use 'GenomeAnalysisTK.jar')
## (a) Create a table 'input.bam.list' of indel candidates
java -Xmx4g -jar $GATK \
-T RealignerTargetCreator \
-R $REF \
-o input.bam.list \
-I myfile.fastq.marked.bam

## (b) Realigns reads around indel candidate (if '$GATK' fails use 'GenomeAnalysisTK.jar')
java -Xmx4g -Djava.io.tmpdir=./tmp \
-jar $GATK \
-I myfile.fastq.bam \
-R $REF \
-T IndelRealigner \
-targetIntervals input.bam.list \
-o myfile.fastq.marked.realigned.bam

## (c) Fix mate information for PE data given updated alignment (if '$PICARD FixMateInformation' fails use 'picard/FixMateInformation.jar')
java -Djava.io.tmpdir=./tmp/flx-auswerter \
-jar $PICARD FixMateInformation \
INPUT=myfile.fastq.marked.realigned.bam \
OUTPUT=myfile.fastq.marked.realigned.fixed.bam \
SO=coordinate \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=true

## Phred Quality score recalibration (if '$GATK' fails use 'GenomeAnalysisTK.jar')
# skipped here since this requires known snps from SNPdb
# 1) Count covariates
# java -Xmx4g -jar $GATK \
# -l INFO \
# -R $REF \
# --DBSNP dbsnp132.txt \
# -I myfile.fastq.marked.realigned.fixed.bam \
# -T CountCovariates \
# -cov ReadGroupCovariate \
# -cov QualityScoreCovariate \
# -cov CycleCovariate \
# -cov DinucCovariate \
# -recalFile input.recal_data.csv

# 2) Table recalibration
# java -Xmx4g -jar $GATK \
# -l INFO \
# -R $REF \
# -I input.marked.realigned.fixed.bam \
# -T TableRecalibration \
# --out input.marked.realigned.fixed.recal.bam \
# -recalFile input.recal_data.csv

## Produce raw SNP calls (if '$GATK' fails use 'GenomeAnalysisTK.jar')
java -Xmx4g -jar $GATK \
-glm BOTH \
-R $REF \
-T UnifiedGenotyper \
-I myfile.fastq.marked.realigned.fixed.bam \
-o vargatk.vcf \
-metrics vargatk.metrics \
-stand_call_conf 50.0 \
-stand_emit_conf 10.0 \
-dcov 1000 \
-A AlleleBalance 

## Filter SNPs (if '$GATK' fails use 'GenomeAnalysisTK.jar')
java -Xmx4g -jar $GATK \
-R $REF \
-T VariantFiltration \
-o vargatk.recalibrated.filtered.vcf \
--variant vargatk.vcf \
--clusterWindowSize 10 \
--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
--filterName "HARD_TO_VALIDATE" \
--filterExpression "DP < 5 " \
--filterName "LowCoverage" \
--filterExpression "QUAL < 30.0 " \
--filterName "VeryLowQual" \
--filterExpression "QUAL > 30.0 && QUAL < 50.0 " \
--filterName "LowQual" \
--filterExpression "QD < 1.5 " \
--filterName "LowQD" \
--filterExpression "SB > -10.0 " \
--filterName "StrandBias" 





