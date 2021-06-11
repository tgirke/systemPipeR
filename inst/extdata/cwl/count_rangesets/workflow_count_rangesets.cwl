class: Workflow
cwlVersion: v1.0
inputs:
  fq1: File
  SampleName: string
  results_path: Directory
outputs:
  count_xls:
    outputSource: /count_xls
    type: File
steps:
  count_rangesets:
    in:
      fq1: fq1
      SampleName: SampleName
      results_path: results_path
    out: '[count_xls]'
    run: count_rangesets/count_rangesets.cwl
