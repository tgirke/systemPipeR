class: Workflow
cwlVersion: v1.0
inputs:
  fq1: File
  SampleName: string
  results_path: Directory
outputs:
  annotated:
    outputSource: /annotated
    type: File
steps:
  annotate_peaks:
    in:
      fq1: fq1
      SampleName: SampleName
      results_path: results_path
    out: '[annotated]'
    run: annotate_peaks/annotate_peaks.cwl
