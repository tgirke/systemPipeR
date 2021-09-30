class: Workflow
cwlVersion: v1.0
inputs:
  fq1: File
  SampleName: string
  results_path: Directory
outputs:
  rundiff:
    outputSource: /rundiff
    type: File
steps:
  rundiff:
    in:
      fq1: fq1
      SampleName: SampleName
      results_path: results_path
    out: '[rundiff]'
    run: rundiff/rundiff.cwl
