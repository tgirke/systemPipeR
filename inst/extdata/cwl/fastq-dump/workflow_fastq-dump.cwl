class: Workflow
cwlVersion: v1.0
inputs:
  defline-seq: string
  defline-qual: string
  outdir: Directory
  maxSpotId: int
  sra: string
outputs:
  out_1:
    outputSource: fastq-dump/out_1
    type: File
  out_2:
    outputSource: fastq-dump/out_2
    type: File
steps:
  fastq-dump:
    in:
      defline-seq: defline-seq
      defline-qual: defline-qual
      outdir: outdir
      maxSpotId: maxSpotId
      sra: sra
    out: '[out_1, out_2]'
    run: fastq-dump/fastq-dump.cwl
