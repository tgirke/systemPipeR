class: Workflow
cwlVersion: v1.0
inputs:
  length: int
  quality: int
  adapter_trim_galore: string
  fq1: File
  SampleName: string
  results_path: Directory
outputs:
  trim_galore:
    outputSource: trim_galore/trim_galore
    type: File
  trim_galore_report:
    outputSource: trim_galore/trim_galore_report
    type: File
steps:
  trim_galore:
    in:
      length: length
      quality: quality
      adapter_trim_galore: adapter_trim_galore
      fq1: fq1
      SampleName: SampleName
      results_path: results_path
    out: '[trim_galore, trim_galore_report]'
    run: trim_galore/trim_galore-se.cwl
