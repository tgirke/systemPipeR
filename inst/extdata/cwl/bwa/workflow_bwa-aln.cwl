class: Workflow
cwlVersion: v1.0
inputs:
  ndis: int
  o: int
  e: int
  l: int
  k: int
  thread: int
  idx_basedir: Directory
  idx_basename: string
  fq1: File
  SampleName: string
  results_path: Directory
outputs:
  bwa_aln_sai:
    outputSource:
    - bwa/bwa_aln_sai
    - aln/bwa_aln_sai
    type: stdout
steps:
  bwa_aln:
    in:
      ndis: ndis
      o: o
      e: e
      l: l
      k: k
      thread: thread
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      fq1: fq1
      SampleName: SampleName
      results_path: results_path
    out: '[bwa_aln_sai]'
    run: bwa/bwa-aln.cwl
