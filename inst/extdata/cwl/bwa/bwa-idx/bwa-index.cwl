################################################################
##                       Bwa-mem-Index                        ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[BWA-MEM](http://bio-bwa.sourceforge.net/bwa.shtml): Burrows-Wheeler Alignment Tool"
label: Last updated 08/1019
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
    listing: [ $(inputs.bwa_idx_basedir) ]
    
arguments:
  - valueFrom: $(inputs.bwa_idx_basedir.path)/$(inputs.bwa_idx_basename)
    position: 1

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  bwa_idx_basename:
    type: string
  bwa_idx_basedir:
     type: Directory

outputs:
  bwa.index_files1:
    type: File
    outputBinding:
      glob: $(inputs.bwa_idx_basedir)/tair10.fasta.amb
  bwa.index_files2:
    type: File
    outputBinding:
      glob: $(inputs.bwa_idx_basedir)/tair10.fasta.ann
  bwa.index_files3:
    type: File
    outputBinding:
      glob: $(inputs.bwa_idx_basedir)/tair10.fasta.bwt
  bwa.index_files4:
    type: File
    outputBinding:
      glob: $(inputs.bwa_idx_basedir)/tair10.fasta.pac
  bwa.index_files5:
    type: File
    outputBinding:
      glob: $(inputs.bwa_idx_basedir)/tair10.fasta.sa
