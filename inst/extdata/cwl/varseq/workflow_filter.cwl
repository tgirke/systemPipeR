class: Workflow
cwlVersion: v1.0
inputs:
  file1: File
  file2: File
  results_path: Directory
  SampleName: string
outputs:
  filter:
    outputSource: /filter
    type: File
steps:
  filter:
    in:
      file1: file1
      file2: file2
      results_path: results_path
      SampleName: SampleName
    out: '[filter]'
    run: varseq/filter.cwl
