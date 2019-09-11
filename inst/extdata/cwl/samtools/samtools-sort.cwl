################################################################
#                        Samtools_Sort                         #
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[Samtools](http://www.htslib.org/doc/samtools.html): Samtools is a suite of programs for interacting with high-throughput sequencing data"
hints:
  SoftwareRequirement:
    packages:
    - package: samtools
      #version: [ 1.9 ]

################################################################
#            baseCommand and arguments definitions             #
################################################################

baseCommand: ["samtools", "sort"]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: -o
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).sorted.bam

################################################################
#               Inputs and Outputs Settings                    #
################################################################

inputs:
  samtools_bam:
    type: File
    label: "BAM format input file"
    inputBinding:
      position: 1
  SampleName:
    label: "Filename to write output to"
    type: string
  thread:
    type: int
    inputBinding:
      prefix: -@
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  samtools_sort_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).sorted.bam
