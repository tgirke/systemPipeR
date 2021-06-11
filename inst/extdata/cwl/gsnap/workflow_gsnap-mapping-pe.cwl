class: Workflow
cwlVersion: v1.0
inputs:
  rscript_gsnap_mapping: File
  fq1: File
  fq2: File
  thread: int
  molecule: string
  max_mismatches: int
  SampleName: string
  results_path: Directory
outputs:
  gsnap_bam:
    outputSource: Rscript/gsnap_bam
    type: File
  gsnap_bam_bai:
    outputSource: Rscript/gsnap_bam_bai
    type: File
steps:
  Rscript:
    in:
      rscript_gsnap_mapping: rscript_gsnap_mapping
      fq1: fq1
      fq2: fq2
      thread: thread
      molecule: molecule
      max_mismatches: max_mismatches
      SampleName: SampleName
      results_path: results_path
    out: '[gsnap_bam, gsnap_bam_bai]'
    run: gsnap/gsnap-mapping-pe.cwl
