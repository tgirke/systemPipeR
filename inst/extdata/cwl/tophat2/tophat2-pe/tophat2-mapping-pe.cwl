################################################################
##                  Tophat2-Paired_end.cwl                    ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[TOPHAT2](https://ccb.jhu.edu/software/tophat/index.shtml): a fast splice junction mapper for RNA-Seq reads"
label: Last updated 09/2019
hints:
  SoftwareRequirement:
    packages:
    - package: tophat2
      version: [ 2.1.1 ]
    - package: bowtie2
      version: [ 2.3.4.1 ]
      
################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [tophat2]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: -g
    valueFrom: '1'
  - prefix: --segment-length
    valueFrom: '25'
  - prefix: -i
    valueFrom: '30'
  - prefix: -I
    valueFrom: '3000'
  - prefix: -p
    valueFrom: $(inputs.thread)
  - prefix: -o
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)
  - valueFrom: $(inputs.idx_basedir.path)/$(inputs.idx_basename)
    
################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  idx_basedir:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory
  idx_basename:
    label: "Basename of the tophat2 index files"
    type: string
  fq1:
    label: "Comma-separated list of files containing mate 1s to be aligned"
    type: File
    inputBinding:
      prefix:
  fq2:
    label: "Comma-separated list of files containing mate 2s to be aligned"
    type: File
    inputBinding:
      prefix:
  SampleName:
    label: "Filename to write output to"
    type: string
  thread:
    label: "Launch NTHREADS parallel search threads"
    type: int
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  accepted_hits_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)/accepted_hits.bam
  junctions_bed:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)/junctions.bed
  insertions_bed:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)/insertions.bed
  deletions_bed:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)/deletions.bed
  align_summary:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)/align_summary.txt

