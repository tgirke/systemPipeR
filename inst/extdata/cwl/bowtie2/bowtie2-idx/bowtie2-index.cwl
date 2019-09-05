################################################################
##                        Bowtie2-Index                       ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): Fast and sensitive read alignment"
label: Last updated 09/2019
hints:
  SoftwareRequirement:
    packages:
    - package: bowtie2
      version: [ 2.3.4.1 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [ bowtie2-build ]

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
  idx_basename:
    label: "Path to the directory containing the index for the reference genome"
    type: string
  idx_basedir:
    label: "Basename of the bowtie2 index files"
    type: Directory

outputs:
  index_files1:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.1.bt2
  index_files2:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.2.bt2
  index_files3:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.3.bt2
  index_files4:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.4.bt2
  index_files5:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.rev.1.bt2
  index_files6:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.rev.2.bt2
