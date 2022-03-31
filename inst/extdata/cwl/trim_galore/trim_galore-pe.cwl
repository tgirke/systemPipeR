################################################################
##                     Trim_Galore.cwl                        ##
################################################################

cwlVersion: v1.2.0
class: CommandLineTool
doc: "[TRIM_GALORE](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)"
label: Last updated 03/2022
hints:
  SoftwareRequirement:
    packages:
    - package: trim_galore
      version: [ 0.6.7 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [trim_galore]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: 
    valueFrom: --phred33
  - prefix: 
    valueFrom: --gzip
  - prefix: 
    valueFrom: --paired
  - prefix: -o
    valueFrom: $(inputs.results_path.path)
  - prefix: --basename
    valueFrom: $(inputs.SampleName)

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  length:
    type: int
    label: "Length of reads to be discarded"
    inputBinding:
      prefix: --length
  quality:
    type: int
    label: "Phred score for quality trimming"
    inputBinding:
      prefix: --quality
  adapter_trim_galore:
    type: string
    label: "Adapter sequence to be trimmed"
    inputBinding:
      prefix: -a
  fq1:
    type: File
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    inputBinding:
      position: 1
  fq2:
    type: File
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    inputBinding:
      position: 1
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  trim_galore_1:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_val_1.fq.gz
  trim_galore_2:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_val_2.fq.gz
  trim_galore_report:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.fq1.basename)_trimming_report.txt

