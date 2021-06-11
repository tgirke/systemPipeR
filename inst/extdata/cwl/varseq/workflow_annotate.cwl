class: Workflow
cwlVersion: v1.0
inputs:
  file1: File
  file2: File
  results_path: Directory
  SampleName: string
  ref_name: string
outputs:
  annotate:
    outputSource: /annotate
    type: File
steps:
  annotate:
    in:
      file1: file1
      file2: file2
      results_path: results_path
      SampleName: SampleName
      ref_name: ref_name
    out: '[annotate]'
    run: varseq/annotate.cwl
