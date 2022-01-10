class: Workflow
cwlVersion: v1.0
inputs:
  data_path: Directory
  ref_name: string
  fq1: File
  fq2: File
  thread: int
  SampleName: string
  results_path: Directory
outputs:
  bwa_men_sam:
    outputSource:
    - bwa/bwa_men_sam
    type: stdout
  samtools-view:
    outputSource: samtools-view/samtools_bam
    type: File
  samtools-sort:
    outputSource: samtools-sort/samtools_sort_bam
    type: File
  samtools-index:
    outputSource: samtools-index/samtools_index
    type: File
steps:
  bwa_mem:
    in:
      data_path: data_path
      ref_name: ref_name
      fq1: fq1
      fq2: fq2
      thread: thread
      SampleName: SampleName
      results_path: results_path
    out: '[bwa_men_sam]'
    run: gatk/bwa-pe.cwl
  samtools-view:
    in:
      samtools_sam: bwa_mem/bwa_men_sam
      SampleName: SampleName
      results_path: results_path
    out: [samtools_bam]
    run: samtools/samtools-view.cwl

  samtools-sort:
    in:
      samtools_bam: samtools-view/samtools_bam
      SampleName: SampleName
      thread: thread
      results_path: results_path
    out: [samtools_sort_bam]
    run: samtools/samtools-sort.cwl
   
  samtools-index:
    in:
      samtools_sort_bam: samtools-sort/samtools_sort_bam
      SampleName: SampleName
      results_path: results_path
    out: [samtools_index]
    run: samtools/samtools-index.cwl
   
