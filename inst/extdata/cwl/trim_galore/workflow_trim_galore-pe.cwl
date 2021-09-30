class: Workflow
cwlVersion: v1.0
inputs:
  length: int
  quality: int
  adapter_trim_galore: string
  fq1: File
  fq2: File
  SampleName: string
  results_path: Directory
outputs:
  trim_galore_1:
    outputSource: trim_galore/trim_galore_1
    type: File
  trim_galore_2:
    outputSource: trim_galore/trim_galore_2
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
      fq2: fq2
      SampleName: SampleName
      results_path: results_path
    out: '[trim_galore_1, trim_galore_2, trim_galore_report]'
    run: trim_galore/trim_galore-pe.cwl
