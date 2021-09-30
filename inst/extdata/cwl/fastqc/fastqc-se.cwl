cwlVersion: v1.0
class: CommandLineTool
baseCommand: fastqc
inputs:
  fq1:
    type: File
    inputBinding:
      prefix: ~
  outdir:
    type: string
    inputBinding:
      prefix: --outdir
  threads:
    type: string
    inputBinding:
      prefix: --threads
  results_path:
    label: "Path to the results directory"
    type: Directory
outputs:
  out_1:
    type: File
    outputBinding:
      glob: $(inputs.outdir.path)/$(inputs.fq1.basename)_fastqc.html
  out_2:
    type: File
    outputBinding:
      glob: $(inputs.outdir.path)/$(inputs.fq1.basename)_fastqc.zip

