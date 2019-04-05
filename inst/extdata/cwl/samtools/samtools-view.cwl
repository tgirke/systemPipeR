################################################################
##                       Samtools_View                        ##
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
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: ["samtools", "view", "-bS"]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: -o
    valueFrom: $(inputs.results_path.basename)/$(inputs.SampleName).bam

################################################################
##              Inputs and Outputs Settings                   ##
################################################################

inputs:
  sam:
    type: File
    label: "SAM format input file"
    inputBinding:
      position: 3
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  samtools_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName).bam
