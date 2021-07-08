#!/bin/bash -l
## Load trim_galore and fastqc
module load	fastqc/0.11.3
module load trim_galore/0.4.2
trim_galore --length 26 --phred33 --fastqc fastq.gz
