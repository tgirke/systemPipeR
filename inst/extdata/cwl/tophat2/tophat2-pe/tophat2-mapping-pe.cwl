################################################################
##                  Tophat2-Paired_end.cwl                    ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[TOPHAT2](https://ccb.jhu.edu/software/tophat/index.shtml): a fast splice junction mapper for RNA-Seq reads"
label: Last updated 07/2019
hints:
  SoftwareRequirement:
    packages:
    - package: tophat
      version: [ 2.1.1 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [tophat2]

arguments:
  - prefix: -g
    valueFrom: '1'
    position: 2
  - prefix: --segment-length
    valueFrom: '25'
    position: 3
  - prefix: -i
    valueFrom: '30'
    position: 4
  - prefix: -I
    valueFrom: '3000'
    position: 5
  - prefix: -p
    valueFrom: $(inputs.thread)
    position: 1
  - prefix: -o
    valueFrom: $(inputs.results_path.basename)/$(inputs.SampleName)
    position: 6
  - position: 7 
    valueFrom: $(inputs.reference_file)
    
################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  reference_file:
    type: File
  tophat2_idx_basedir:
    type: Directory
  tophat2_idx_basename:
    type: string
  fq1:
    type: File
    inputBinding:
      prefix:
      position: 8
  fq2:
    type: File
    inputBinding:
      prefix:
      position: 8
  SampleName:
    type: string
  thread:
    type: int
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  accepted_hits_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName)/accepted_hits.bam
  junctions_bed:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName)/junctions.bed
  insertions_bed:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName)/insertions.bed
  deletions_bed:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName)/deletions.bed
  align_summary:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName)/align_summary.txt

