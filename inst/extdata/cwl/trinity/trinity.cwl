################################################################
##                  Trinity.cwl                               ##
################################################################

cwlVersion: v1.2
class: CommandLineTool
label: Last updated 03/2022
hints:
  SoftwareRequirement:
    packages:
    - package: trinity-rnaseq
      version: [ 2.13.2 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [Trinity]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]


################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  seqType:
    type: string
    inputBinding:
      prefix: --seqType
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
    type: string
    inputBinding:
      prefix: --SS_lib_type
  trimmomatic:
    type: string
    inputBinding:
      prefix: --trimmomatic --quality_trimming_params

outputs:
  trinity_results:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/Trinity.fasta

