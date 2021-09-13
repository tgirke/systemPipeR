class: Workflow
cwlVersion: v1.0
inputs:
  fq1: File
  fq2: File
  outdir: Directory
  threads: string
outputs:
  out_1:
    outputSource: fastqc/out_1
    type: File
  out_2:
    outputSource: fastqc/out_2
    type: File
  out_3:
    outputSource: fastqc/out_3
    type: File
  out_4:
    outputSource: fastqc/out_4
    type: File
steps:
  fastqc:
    in:
      fq1: fq1
      fq2: fq2
      outdir: outdir
      threads: threads
    out: '[out_1, out_2, out_3, out_4]'
    run: fastqc/fastqc.cwl
