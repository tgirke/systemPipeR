################################################################
##                        Hisat2-Index                        ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml): graph-based alignment of next generation sequencing reads to a population of genomes"
label: Last updated 03/2022
hints:
  SoftwareRequirement:
    packages:
    - package: hisat2
      version: [ 2.2.1 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [hisat2-build]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.idx_basedir) ]

arguments:
  - valueFrom: $(inputs.idx_basedir.path)/$(inputs.idx_basename)
  - valueFrom: $(inputs.idx_basedir.path)/$(inputs.idx_basename)
    
################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  idx_basedir:
    label: "Path to the directory containing the reference genome"
    type: Directory
  idx_basename:
    label: "Reference genome basename"
    type: string

outputs:
  hisat2_index1:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/$(inputs.idx_basename).1.ht2
  hisat2_index2:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/$(inputs.idx_basename).2.ht2
  hisat2_index3:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/$(inputs.idx_basename).3.ht2
  hisat2_index4:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/$(inputs.idx_basename).4.ht2
  hisat2_index5:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/$(inputs.idx_basename).5.ht2
  hisat2_index6:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/$(inputs.idx_basename).6.ht2
  hisat2_index7:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/$(inputs.idx_basename).7.ht2
  hisat2_index8:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/$(inputs.idx_basename).8.ht2
