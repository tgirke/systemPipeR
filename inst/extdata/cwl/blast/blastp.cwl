cwlVersion: v1.0
class: CommandLineTool
baseCommand: blastp
inputs:
  query:
    type: File
    inputBinding:
      prefix: -query
  database_name:
    type: File
    inputBinding:
      prefix: -db
  evalue:
    type: int
    inputBinding:
      prefix: -evalue
  outfmt:
    type: int
    inputBinding:
      prefix: -outfmt
  out:
    type: File
    inputBinding:
      prefix: -out
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
