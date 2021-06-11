class: Workflow
cwlVersion: v1.0
inputs:
  fq1: File
  SampleName: string
  results_path: Directory
outputs:
  gunzip_file:
    outputSource:
    - gunzip/gunzip_file
    - -c/gunzip_file
    type: stdout
steps:
  gunzip_-c:
    in:
      fq1: fq1
      SampleName: SampleName
      results_path: results_path
    out: '[gunzip_file]'
    run: gunzip/gunzip.cwl
