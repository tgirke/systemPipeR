#!/bin/bash -l

java=$1
input_vcf=$2
output_vcf=$3
ref=$4

gatk VariantFiltration \
    --java-options "${java}" \
    -V ${input_vcf} \
    -O ${output_vcf} \
    -R ${ref} \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \