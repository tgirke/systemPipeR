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

baseCommand: ["bcftools","call", "-cv", "-Ov"]

arguments:
  - prefix: -o
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_bcf.vcf
  - prefix:
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_raw.bcf.gz

################################################################
##               Inputs and Outputs Settings                  ##                                                              
################################################################

inputs:
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  bcftools_call:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_bcf.vcf
