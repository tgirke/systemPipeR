################################################################
##                   Kallisto-Index                           ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: 
label: Last updated 03/2022
hints:
  SoftwareRequirement:
    packages:
    - package: kallisto
      version: [ 0.46.1 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: ["kallisto", "index"]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.idx_basedir) ]
    
arguments:
  - prefix: -i
    valueFrom: $(inputs.idx_basedir.path)/$(inputs.kallisto_idx)
  - prefix: 
    valueFrom: $(inputs.idx_basedir.path)/$(inputs.idx_basename)

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  idx_basename:
    label: "Basename of the kallisto index files"
    type: string
  idx_basedir:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory
  kallisto_idx: 
    label: "Name for resulting kallisto index"
    type: string

outputs:
  kallisto_index:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/$(inputs.kallisto_idx)
