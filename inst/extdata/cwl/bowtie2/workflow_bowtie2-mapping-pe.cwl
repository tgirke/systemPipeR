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
  bowtie2_sam:
    outputSource: bowtie2/bowtie2_sam
    type: File
steps:
  bowtie2:
    in:
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      fq1: fq1
      fq2: fq2
      thread: thread
      SampleName: SampleName
      results_path: results_path
    out: '[bowtie2_sam]'
    run: bowtie2/bowtie2-mapping-pe.cwl
