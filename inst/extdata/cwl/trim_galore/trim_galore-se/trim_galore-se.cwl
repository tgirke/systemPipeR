################################################################
##                     Trim_Galore.cwl                        ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[TRIM_GALORE](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)"
hints:
  SoftwareRequirement:
    packages:
    - package: trim_galore
      version: [ 0.4.2 ]
    - package: fastqc
      version: [0.11.3]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [trim_galore]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: --phred33
    position: 1
##  - prefix: --fastqc
##    position: 1
  - prefix: --gzip
    position: 2
  - prefix: -o
    valueFrom: $(inputs.results_path.basename)
    position: 3

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  length:
    type: int
    inputBinding:
      prefix: --length
      position: 4
  quality:
    type: int
    inputBinding:
      prefix: --quality
      position: 5
  adapter:
    type: string
    inputBinding:
      prefix: -a
      position: 6
  fq1:
    type: File
    inputBinding:
      position: 7
  SampleName:
    type: string
  results_path:
    type: Directory


outputs:
  trim_galore:
    type: File
    outputBinding:
      glob:$(inputs.results_path.basename)/$(inputs.SampleName)_trimmed.fq.gz
  trimming_report:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName).fastq.gz_trimming_report.txt
  fastqc_html:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName)_trimmed_fastqc.html  
  fastqc_zip:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName)_trimmed_fastqc.zip
