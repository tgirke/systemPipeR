################################################################
##                  Fasta-dict-Index                          ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: ""
label: Last updated 10/1019
hints:
  SoftwareRequirement:
    packages:
    - package: gatk
      version: [ 4.1.1.0 ]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.data_path) ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [gatk, CreateSequenceDictionary]

arguments:
  - prefix: -R
    valueFrom: $(inputs.data_path.path)/$(inputs.ref_name)
    position: 1

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  ref_name:
    type: string
  data_path:
    type: Directory

outputs:
  fasta_dict:
    type: File
    outputBinding:
      glob: $(inputs.data_path.path)/$(inputs.ref_name.nameroot).dict
