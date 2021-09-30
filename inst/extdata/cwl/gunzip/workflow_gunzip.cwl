class: Workflow
cwlVersion: v1.0
inputs:
  file: File
  SampleName: string
  ext: string
  results_path: Directory
outputs:
  gunzip_file:
    outputSource: gunzip/gunzip_file
    type: stdout
steps:
  gunzip:
    in:
      file: file
      SampleName: SampleName
      ext: ext
      results_path: results_path
    out: '[gunzip_file]'
    run: gunzip/gunzip.cwl
