################################################################
##                        Samtools_Sort                       ##
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
##          baseCommand and arguments definitions            ##
################################################################

baseCommand: ["samtools", "sort"]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]
    
arguments:
  - prefix: -o
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_sort.bam

  - prefix:
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_mark.bam
    
################################################################
##               Inputs and Outputs Settings                   ##
################################################################

inputs:
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  samtools_sort_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_sort.bam
