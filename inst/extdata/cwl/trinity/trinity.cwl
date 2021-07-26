################################################################
##                  Trinity.cwl                               ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
label: Last updated 07/2019
hints:
  SoftwareRequirement:
    packages:
    - package: trinity-rnaseq
      version: [ 2.8.4 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [Trinity]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: --seqType
    valueFrom: 'fq'
  - prefix: --trimmomatic
  - prefix: --quality_trimming_params
    valueFrom: '''SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:35'''

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  max_memory:
    type: string
    inputBinding:
      prefix: --max_memory
  samples_file:
    type: File
    inputBinding:
      prefix: --samples_file
  thread:
    type: int
    inputBinding:
      prefix: --CPU
  results_path:
    type: Directory
    inputBinding:
      prefix: --output
  ss_lib_type:
    type: string?
  SampleName:
    type: string
    
outputs:
  trinity_results:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/Trinity.fasta
