class: Workflow
cwlVersion: v1.0
inputs:
  file: File
  SampleName: string
  ext: string
  results_path: Directory
outputs:
  gzip_file:
    outputSource: gzip/gzip_file
    type: stdout
steps:
  gzip:
    in:
      file: file
      SampleName: SampleName
      ext: ext
      results_path: results_path
    out: '[gzip_file]'
    run: gunzip/gzip.cwl
