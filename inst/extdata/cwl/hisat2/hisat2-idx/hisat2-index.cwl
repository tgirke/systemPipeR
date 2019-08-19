################################################################
##                        Hisat2-Index                        ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml): graph-based alignment of next generation sequencing reads to a population of genomes"
label: Last updated 02/2019
hints:
  SoftwareRequirement:
    packages:
    - package: hisat2
      version: [ 2.1.0 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [hisat2-build]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.hisat2_idx_basedir) ]

arguments:
  - valueFrom: $(inputs.hisat2_idx_basedir.path)/$(inputs.hisat2_idx_basename)
    position: 1

  - valueFrom: $(inputs.hisat2_idx_basedir.path)/$(inputs.hisat2_idx_basename)
    position: 2
    
################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  hisat2_idx_basedir:
    label: "Path to the directory the index for the reference genome"
    type: Directory
  hisat2_idx_basename:
    label: "Basename of the hisat2 index files"
    type: string

outputs:
  hisat2_index1:
    type: File
    outputBinding:
      glob: $(inputs.hisat2_idx_basedir.path)/$(inputs.hisat2_idx_basename).1.ht2
  hisat2_index2:
    type: File
    outputBinding:
      glob: $(inputs.hisat2_idx_basedir.path)/$(inputs.hisat2_idx_basename).2.ht2
  hisat2_index3:
    type: File
    outputBinding:
      glob: $(inputs.hisat2_idx_basedir.path)/$(inputs.hisat2_idx_basename).3.ht2
  hisat2_index4:
    type: File
    outputBinding:
      glob: $(inputs.hisat2_idx_basedir.path)/$(inputs.hisat2_idx_basename).4.ht2
  hisat2_index5:
    type: File
    outputBinding:
      glob: $(inputs.hisat2_idx_basedir.path)/$(inputs.hisat2_idx_basename).5.ht2
  hisat2_index6:
    type: File
    outputBinding:
      glob: $(inputs.hisat2_idx_basedir.path)/$(inputs.hisat2_idx_basename).6.ht2
  hisat2_index7:
    type: File
    outputBinding:
      glob: $(inputs.hisat2_idx_basedir.path)/$(inputs.hisat2_idx_basename).7.ht2
  hisat2_index8:
    type: File
    outputBinding:
      glob: $(inputs.hisat2_idx_basedir.path)/$(inputs.hisat2_idx_basename).8.ht2
