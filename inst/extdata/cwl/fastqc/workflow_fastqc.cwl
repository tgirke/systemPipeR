class: Workflow
cwlVersion: v1.0
inputs:
  fq1: File
  outdir: Directory
  threads: string
outputs:
  out_1:
    outputSource: fastqc/out_1
    type: File
  out_2:
    outputSource: fastqc/out_2
    type: File
steps:
  fastqc:
    in:
      fq1: fq1
      outdir: outdir
      threads: threads
    out: '[out_1, out_2]'
    run: fastqc/fastqc.cwl
