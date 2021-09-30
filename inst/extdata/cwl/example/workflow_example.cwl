class: Workflow
cwlVersion: v1.0
inputs:
  message: string
  SampleName: string
  results_path: Directory
outputs:
  string:
    outputSource: echo/string
    type: stdout
steps:
  echo:
    in:
      message: message
      SampleName: SampleName
      results_path: results_path
    out: '[string]'
    run: example/example.cwl
