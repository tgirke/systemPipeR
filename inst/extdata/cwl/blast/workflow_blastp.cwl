class: Workflow
cwlVersion: v1.0
inputs:
  query: File
  database_name: File
  evalue: int
  outfmt: int
  out: File
outputs:
  out:
    outputSource: blastp/out
    type: File
steps:
  blastp:
    in:
      query: query
      database_name: database_name
      evalue: evalue
      outfmt: outfmt
      out: out
    out: '[out]'
    run: blast/blastp.cwl
