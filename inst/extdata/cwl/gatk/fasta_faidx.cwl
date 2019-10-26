################################################################
##                  Fasta-faidx-Index                         ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: ""
label: Last updated 10/1019
hints:
  SoftwareRequirement:
    packages:
    - package: samtools
      version: [ 1.9 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: ["samtools", "faidx"]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.data_path) ]
    
arguments:
  - valueFrom: $(inputs.data_path.path)/$(inputs.ref_name)

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  ref_name:
    type: string
  data_path:
    type: Directory

outputs:
  fasta_fai:
    type: File
    outputBinding:
      glob: $(inputs.data_path.path)/$(inputs.ref_name).fai
