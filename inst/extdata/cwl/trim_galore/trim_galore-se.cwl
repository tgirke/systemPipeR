################################################################
##                     Trim_Galore.cwl                        ##
################################################################

cwlVersion: v1.2.0
class: CommandLineTool
doc: "[TRIM_GALORE](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)"
label: Last updated 09/2021
hints:
  SoftwareRequirement:
    packages:
    - package: trim_galore
      version: [ 0.6.5 ]

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
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory


outputs:
  trim_galore:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_trimmed.fq.gz
  trim_galore_report:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.fq1.basename)_trimming_report.txt

