class: Workflow
cwlVersion: v1.0
inputs:
  rscript_test: File
  results_path: Directory
outputs:
  test1:
    outputSource: Rscript/test1
    type: File
steps:
  Rscript:
    in:
      rscript_test: rscript_test
      results_path: results_path
    out: '[test1]'
    run: test/test.cwl
