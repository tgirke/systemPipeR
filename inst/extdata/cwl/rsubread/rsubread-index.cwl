################################################################
##                     Rsubread-index                         ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[Rsubread - v1.34.7](https://www.biorxiv.org/content/10.1101/377762v2): alignment and quantification of RNA sequencing reads"
label: Last updated 09/2019

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [Rscript]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.idx_basedir) ]

arguments:
  - prefix: --vanilla
    valueFrom: $(inputs.rscript_rsubread_idx)
  - prefix: --reference
    valueFrom: $(inputs.idx_basedir.path)/$(inputs.idx_basename)
                
################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  rscript_rsubread_idx: 
    type: File
  idx_basedir:
    label: "Path to the directory the index for the reference genome"
    type: Directory
  idx_basename:
    label: "Basename of the rsubread index files"
    type: string

outputs:
  rsubread_index1:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.00.b.array
  rsubread_index2:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.00.b.tab
  rsubread_index3:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.files
  rsubread_index4:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.log
  rsubread_index5:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.path)/tair10.fasta.reads
