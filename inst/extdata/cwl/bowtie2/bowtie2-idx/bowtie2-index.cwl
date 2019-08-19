################################################################
##                        Bowtie2-Index                        ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): Fast and sensitive read alignment"
label: Last updated 07/2019
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
    listing: [ $(inputs.bowtie2_idx_basedir) ]
    
arguments:
  - valueFrom: $(inputs.bowtie2_idx_basedir.path)/$(inputs.bowtie2_idx_basename)
    position: 1

  - valueFrom: $(inputs.bowtie2_idx_basedir.path)/$(inputs.bowtie2_idx_basename)
    position: 2
    
################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  bowtie2_idx_basename:
    type: string
  bowtie2_idx_basedir:
     type: Directory

outputs:
##  index_files:
##    type:
##      type: array
##      items: Directory
##    outputBinding:
##      glob: "*"

  index_files1:
    type: File
    outputBinding:
      glob: $(inputs.bowtie2_idx_basedir)/tair10.fasta.1.bt2
  index_files2:
    type: File
    outputBinding:
      glob: $(inputs.bowtie2_idx_basedir)/tair10.fasta.2.bt2
  index_files3:
    type: File
    outputBinding:
      glob: $(inputs.bowtie2_idx_basedir)/tair10.fasta.3.bt2
  index_files4:
    type: File
    outputBinding:
      glob: $(inputs.bowtie2_idx_basedir)/tair10.fasta.4.bt2
  index_files5:
    type: File
    outputBinding:
      glob: $(inputs.bowtie2_idx_basedir)/tair10.fasta.rev.1.bt2
  index_files6:
    type: File
    outputBinding:
      glob: $(inputs.bowtie2_idx_basedir)/tair10.fasta.rev.2.bt2
