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
    listing: [ $(inputs.idx_basedir) ]
    
arguments:
  - valueFrom: $(inputs.idx_basedir.path)/$(inputs.idx_basename)

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  idx_basename:
    label: "Basename of the bwa index files"
    type: string
  idx_basedir:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory

outputs:
  bwa.index_files1:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.amb
  bwa.index_files2:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.ann
  bwa.index_files3:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.bwt
  bwa.index_files4:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.pac
  bwa.index_files5:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.sa
