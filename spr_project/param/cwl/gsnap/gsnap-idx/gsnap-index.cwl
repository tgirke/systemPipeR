################################################################
##                       gsnap-index                         ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "gmapR (1.26.0): An R interface to the GMAP/GSNAP/GSTRUCT suite - DOI: 10.18129/B9.bioc.gmapR"
label: Last updated 09/2019

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [Rscript]

arguments:
  - prefix: --vanilla
    valueFrom: $(inputs.rscript_gsnap_idx)
  - prefix: --reference
    valueFrom: $(inputs.idx_basedir.path)/$(inputs.idx_basename)
                
################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  rscript_gsnap_idx:
    type: File
  idx_basedir:
    label: "Path to the directory the index for the reference genome"
    type: Directory
  idx_basename:
    label: "Basename of the gsnap index files"
    type: string

outputs:
  gsnap_index1:
    type: File #TODO: fix SYSargs2 class to accepts
    outputBinding:
      glob: $(inputs.idx_basedir.basename)/gmap_tair10chr
  gsnap_index2:
    type: File
    outputBinding:
      glob: $(inputs.idx_basedir.basename)/gmapGenome.RData


