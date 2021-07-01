cwlVersion: v1.0
class: CommandLineTool
arguments:
- prefix: --gzip
- prefix: --split-files
baseCommand: fastq-dump
inputs:
  defline-seq:
    type: string
    inputBinding:
      prefix: --defline-seq
  defline-qual:
    type: string
    inputBinding:
      prefix: --defline-qual
  outdir:
    type: Directory
    inputBinding:
      prefix: --outdir
  maxSpotId:
    type: int
    inputBinding:
      prefix: --maxSpotId
  sra:
    type: string
    inputBinding:
      prefix: ~
outputs:
  out_1:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.sra)_1.fastq.gz
  out_2:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.sra)_2.fastq.gz
