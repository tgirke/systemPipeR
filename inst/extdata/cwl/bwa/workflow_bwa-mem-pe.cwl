class: Workflow
cwlVersion: v1.0
inputs:
  idx_basedir: Directory
  idx_basename: string
  fq1: File
  fq2: File
  thread: int
  SampleName: string
  results_path: Directory
outputs:
  bwa_men_sam:
    outputSource:
    - bwa/bwa_men_sam
    - mem/bwa_men_sam
    type: stdout
steps:
  bwa_mem:
    in:
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      fq1: fq1
      fq2: fq2
      thread: thread
      SampleName: SampleName
      results_path: results_path
    out: '[bwa_men_sam]'
    run: bwa/bwa-mem-pe.cwl
