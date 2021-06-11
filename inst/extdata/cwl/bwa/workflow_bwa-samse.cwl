class: Workflow
cwlVersion: v1.0
inputs:
  idx_basedir: Directory
  idx_basename: string
  sai: File
  fq1: File
  SampleName: string
  results_path: Directory
outputs:
  bwa_samse_sam:
    outputSource:
    - bwa/bwa_samse_sam
    - samse/bwa_samse_sam
    type: stdout
steps:
  bwa_samse:
    in:
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      sai: sai
      fq1: fq1
      SampleName: SampleName
      results_path: results_path
    out: '[bwa_samse_sam]'
    run: bwa/bwa-samse.cwl
