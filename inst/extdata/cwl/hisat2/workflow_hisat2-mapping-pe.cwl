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
  hisat2_sam:
    outputSource: hisat2/hisat2_sam
    type: File
steps:
  hisat2:
    in:
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      fq1: fq1
      fq2: fq2
      thread: thread
      SampleName: SampleName
      results_path: results_path
    out: '[hisat2_sam]'
    run: hisat2/hisat2-mapping-pe.cwl
