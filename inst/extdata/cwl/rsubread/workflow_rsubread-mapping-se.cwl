class: Workflow
cwlVersion: v1.0
inputs:
  rscript_rsubread_mapping: File
  idx_basedir: Directory
  idx_basename: string
  fq1: File
  thread: int
  SampleName: string
  results_path: Directory
outputs:
  rsubread_sam:
    outputSource: Rscript/rsubread_sam
    type: File
  rsubread_vcf:
    outputSource: Rscript/rsubread_vcf
    type: File
  rsubread_summary:
    outputSource: Rscript/rsubread_summary
    type: File
steps:
  Rscript:
    in:
      rscript_rsubread_mapping: rscript_rsubread_mapping
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      fq1: fq1
      thread: thread
      SampleName: SampleName
      results_path: results_path
    out: '[rsubread_sam, rsubread_vcf, rsubread_summary]'
    run: rsubread/rsubread-mapping-se.cwl
