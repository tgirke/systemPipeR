################################################################
##                        Bwa-Index                           ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[BWA-MEM](http://bio-bwa.sourceforge.net/bwa.shtml): Burrows-Wheeler Alignment Tool"
label: Last updated 09/1019
hints:
  SoftwareRequirement:
    packages:
    - package: bwa
      version: [ 0.7.17 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: ["bwa", "index", "-a", "bwtsw"]

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
    label: "Basename of the bwa index files"
    type: string
  data_path:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory

outputs:
  bwa.index_files1:
    type: File
    outputBinding:
      glob: $(inputs.data_path.path)/$(inputs.ref_name).amb
  bwa.index_files2:
    type: File
    outputBinding:
      glob: $(inputs.data_path.path)/$(inputs.ref_name).ann
  bwa.index_files3:
    type: File
    outputBinding:
      glob: $(inputs.data_path.path)/$(inputs.ref_name).bwt
  bwa.index_files4:
    type: File
    outputBinding:
      glob: $(inputs.data_path.path)/$(inputs.ref_name).pac
  bwa.index_files5:
    type: File
    outputBinding:
      glob: $(inputs.data_path.path)/$(inputs.ref_name).sa
