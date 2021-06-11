class: Workflow
cwlVersion: v1.0
inputs:
  idx_basedir: Directory
  idx_basename: string
  fq1: File
  thread: int
  SampleName: string
  results_path: Directory
outputs:
  bwa_mem_sam:
    outputSource:
    - bwa/bwa_mem_sam
    - mem/bwa_mem_sam
    type: stdout
steps:
  bwa_mem:
    in:
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      fq1: fq1
      thread: thread
      SampleName: SampleName
      results_path: results_path
    out: '[bwa_mem_sam]'
    run: bwa/bwa-mem-se.cwl
