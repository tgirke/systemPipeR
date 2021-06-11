class: Workflow
cwlVersion: v1.0
inputs:
  fq1: File
  fq2: File
  SampleName: string
  results_path: Directory
outputs:
  trim-pe_1:
    outputSource: /trim-pe_1
    type: File
  trim-pe_2:
    outputSource: /trim-pe_2
    type: File
steps:
  trim-pe:
    in:
      fq1: fq1
      fq2: fq2
      SampleName: SampleName
      results_path: results_path
    out: '[trim-pe_1, trim-pe_2]'
    run: preprocessReads/trim-pe.cwl
