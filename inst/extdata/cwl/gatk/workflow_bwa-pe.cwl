class: Workflow
cwlVersion: v1.0
inputs:
  data_path: Directory
  ref_name: string
  fq1: File
  fq2: File
  bwa_threads: int
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
      data_path: data_path
      ref_name: ref_name
      fq1: fq1
      fq2: fq2
      bwa_threads: bwa_threads
      SampleName: SampleName
      results_path: results_path
    out: '[bwa_men_sam]'
    run: gatk/bwa-pe.cwl
