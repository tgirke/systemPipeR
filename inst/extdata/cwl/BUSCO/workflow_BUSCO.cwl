class: Workflow
cwlVersion: v1.0
inputs:
  assembly: File
  lineage: Directory
  outputName: string
  cpu: int?
  results_path: Directory
outputs:
  busco_blast:
    outputSource: run_BUSCO.py/busco_blast
    type: File
  busco_table:
    outputSource: run_BUSCO.py/busco_table
    type: File
  busco_hmmer:
    outputSource: run_BUSCO.py/busco_hmmer
    type: File
  busco_missing:
    outputSource: run_BUSCO.py/busco_missing
    type: File
  busco_summary:
    outputSource: run_BUSCO.py/busco_summary
    type: File
  busco_proteins:
    outputSource: run_BUSCO.py/busco_proteins
    type: File
steps:
  run_BUSCO.py:
    in:
      assembly: assembly
      lineage: lineage
      outputName: outputName
      cpu: cpu
      results_path: results_path
    out: '[busco_blast, busco_table, busco_hmmer, busco_missing, busco_summary, busco_proteins]'
    run: BUSCO/BUSCO.cwl
