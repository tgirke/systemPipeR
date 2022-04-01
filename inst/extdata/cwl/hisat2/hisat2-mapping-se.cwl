################################################################
##                   Hisat2-Single_Read                       ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml): graph-based alignment of next generation sequencing reads to a population of genomes"
label: Last updated 03/2022
hints:
  SoftwareRequirement:
    packages:
    - package: hisat2
      version: [ 2.1.0 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [hisat2]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: -S
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).sam
  - prefix: -x
    valueFrom: $(inputs.idx_basedir.path)/$(inputs.idx_basename)
  - prefix: -k
    valueFrom: '1'
  - prefix: --min-intronlen
    valueFrom: '30'
  - prefix: --max-intronlen
    valueFrom: '3000'

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  idx_basedir:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory
  idx_basename:
    label: "Basename of the hisat2 index files"
    type: string
  fq1:
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    type: File
    inputBinding:
      prefix: -U
      itemSeparator: ","
  thread:
    label: "Launch NTHREADS parallel search threads"
    type: int
    inputBinding:
      prefix: --threads
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  hisat2_sam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).sam
