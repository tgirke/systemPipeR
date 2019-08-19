#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: [tophat2]

arguments:
  - position: 7
    valueFrom: $(inputs.reference_file)
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
  - prefix: --output-dir
    valueFrom: $(inputs.output_filename)
    position: 6

inputs:
  reference_file:
    type: File
  index_basedir:
    type: Directory
  index_basename:
    type: string
  fq1:
    type: File
    inputBinding:
      position: 8
  fq2:
    type: File
    inputBinding:
      position: 8
  output_filename:
    type: string
  nthreads:
    type: int
    inputBinding:
      prefix: -p
      position: 1

outputs:
  accepted_hits_bam:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)/accepted_hits.bam
  junctions_bed:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)/junctions.bed
  insertions_bed:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)/insertions.bed
  deletions_bed:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)/deletions.bed
  align_summary:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)/align_summary.txt



