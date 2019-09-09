################################################################
##                   Rsubread-Single_Read                     ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[Rsubread - v1.34.7](https://www.biorxiv.org/content/10.1101/377762v2): alignment and quantification of RNA sequencing reads"
label: Last updated 09/2019

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [Rscript]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: --vanilla
    valueFrom: $(inputs.rscript_rsubread_mapping)
  - prefix: --output_file
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).sam
  - prefix: --input_format
    valueFrom: "FASTQ"
  - prefix: --output_format
    valueFrom: "SAM"
  - prefix: --reference
    valueFrom: $(inputs.idx_basedir.basename)/$(inputs.idx_basename)
  - prefix: --indels
    valueFrom: "1"
  - prefix: --TH1
    valueFrom: "2"
  - prefix: --maxMismatches
    valueFrom: "3"

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  rscript_rsubread_mapping:
    type: File
  idx_basedir:
    label: "Path to the directory the index for the reference genome"
    type: Directory
  idx_basename:
    label: "Basename of the rsubread index files"
    type: string
  fq1:
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    type: File
    inputBinding:
      prefix: --readfile1
  thread:
    type: int
    inputBinding:
      prefix: --nthreads
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  rsubread_sam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).sam
  rsubread_vcf:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).sam.indel.vcf
  rsubread_summary:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).sam.summary
