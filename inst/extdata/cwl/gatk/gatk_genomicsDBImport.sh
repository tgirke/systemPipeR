#!/bin/bash -l

reference=$1
results_dir=$2
java_options=$3
DB_path=$4
shift;shift;shift;shift

## extract chromosomes from ref
chrom_list=$(awk -F ' ' '{print $1}' ${reference}.fai | \
                  sed 's/^/-L /' | \
                  paste -d " " -s | \
                  sed 's/^//')

## parse all g.vcf files names
files=$(ls ${results_dir} | grep g.vcf.gz$)

vcf_files=$(echo ${files} | \
                sed "s| | -V ${results_dir}/|g" | \
                sed "s|^|-V ${results_dir}/|" )

##echo ${chrom_list}
##echo ${vcf_files}
gatk --java-options "${java_options}" \
    GenomicsDBImport \
    --genomicsdb-workspace-path ${DB_path} \
    ${vcf_files} \
    ${chrom_list} \
    ${@}