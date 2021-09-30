class: Workflow
cwlVersion: v1.0
inputs:
  ref_name: string
  data_path: Directory
outputs:
  fasta_fai:
    outputSource:
    - samtools/fasta_fai
    - faidx/fasta_fai
    type: File
steps:
  samtools_faidx:
    in:
      ref_name: ref_name
      data_path: data_path
    out: '[fasta_fai]'
    run: gatk/fasta_faidx.cwl
