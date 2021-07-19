cwlVersion: v1.0
class: CommandLineTool
baseCommand: fastqc
inputs:
  fq1:
    type: File
    inputBinding:
      prefix: ~
  outdir:
    type: Directory
    inputBinding:
      prefix: --outdir
  threads:
    type: string
    inputBinding:
      prefix: --threads
outputs:
  out_1:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.outdir)/$(inputs.fq1.basename)_fastqc.html
  out_2:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.outdir)/$(inputs.fq1.basename)_fastqc.zip
