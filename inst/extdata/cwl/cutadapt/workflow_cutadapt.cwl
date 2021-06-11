class: Workflow
cwlVersion: v1.0
inputs:
  adapter: string
  outfile: File
  minimum-length: string
  maximum-length: string
  cores: string
  fq1: File
  cutadapt: Directory
  SampleName: string
outputs:
  outfile:
    outputSource: cutadapt/outfile
    type: File
steps:
  cutadapt:
    in:
      adapter: adapter
      outfile: outfile
      minimum-length: minimum-length
      maximum-length: maximum-length
      cores: cores
      fq1: fq1
      cutadapt: cutadapt
      SampleName: SampleName
    out: '[outfile]'
    run: cutadapt/cutadapt.cwl
