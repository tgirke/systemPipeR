################################################################
##                         bcftools.cwl                       ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: ""
label: Last updated 03/2022
hints:
  SoftwareRequirement:
    packages:
    - package: bcftools
      version: [ 1.15 ]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: ["bcftools", "mpileup", "-Ob"]

arguments:
  - prefix: -f
    valueFrom: $(inputs.data_path.path)/$(inputs.idx_basename)
    position: 1
  - prefix: -o
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_raw.bcf.gz
    position: 2
  - prefix:
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_sort.bam
    position: 3  

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  idx_basename:
    type: string
  data_path:
    label: "Basename of the index files"
    type: Directory
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  bcftools:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_raw.bcf.gz
