class: Workflow
cwlVersion: v1.0
inputs:
  ref_name: string
  data_path: Directory
outputs:
  fasta_dict:
    outputSource:
    - gatk/fasta_dict
    - CreateSequenceDictionary/fasta_dict
    type: File
steps:
  gatk_CreateSequenceDictionary:
    in:
      ref_name: ref_name
      data_path: data_path
    out: '[fasta_dict]'
    run: gatk/fasta_dict.cwl
