################################################################
##                       Samtools_Index                       ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[Samtools](http://www.htslib.org/doc/samtools.html): Samtools is a suite of programs for interacting with high-throughput sequencing data"
label: Last updated 02/2019
#hints:
 # SoftwareRequirement:
  #  packages:
   # - package: samtools
    #  version: [ 1.9 ]

################################################################
##            baseCommand and arguments definitions           ##
################################################################

baseCommand: ["samtools", "index", "-b"]

requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.samtools_sort_bam)
      - $(inputs.results_path)

################################################################
##             Inputs and Outputs Settings                    ##
################################################################

inputs:
  samtools_sort_bam:
    type: File
    inputBinding:
      position: 1
  SampleName:
    label: "Sample name"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  samtools-index:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName).sorted.bam.bai

