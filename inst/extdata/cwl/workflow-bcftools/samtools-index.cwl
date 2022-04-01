################################################################
##                       Samtools_Index                       ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[Samtools](http://www.htslib.org/doc/samtools.html): Samtools is a suite of programs for interacting with high-throughput sequencing data"
label: Last updated 03/2022
hints:
  SoftwareRequirement:
    packages:
    - package: samtools
      version: [ 1.14 ]

################################################################
##            baseCommand and arguments definitions           ##
################################################################

baseCommand: ["samtools", "index"]
      
requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: 
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_sort.bam

################################################################
##             Inputs and Outputs Settings                    ##
################################################################

inputs:
  SampleName:
    label: "Sample name"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  samtools_index:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_sort.bam.bai
