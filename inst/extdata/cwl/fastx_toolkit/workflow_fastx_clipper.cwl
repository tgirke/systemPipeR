class: Workflow
cwlVersion: v1.0
inputs:
  gunzip: File
  adaptor_fastx_toolkit: string
  SampleName: string
  results_path: Directory
outputs:
  fastx_clipper_file:
    outputSource: fastx_clipper/fastx_clipper_file
    type: File
steps:
  fastx_clipper:
    in:
      gunzip: gunzip
      adaptor_fastx_toolkit: adaptor_fastx_toolkit
      SampleName: SampleName
      results_path: results_path
    out: '[fastx_clipper_file]'
    run: fastx_toolkit/fastx_clipper.cwl
