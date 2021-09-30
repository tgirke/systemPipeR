class: Workflow
cwlVersion: v1.0
inputs:
  fq1: File
  SampleName: string
  results_path: Directory
outputs:
  trim-se:
    outputSource: /trim-se
    type: File
steps:
  trim-se:
    in:
      fq1: fq1
      SampleName: SampleName
      results_path: results_path
    out: '[trim-se]'
    run: preprocessReads/trim-se.cwl
