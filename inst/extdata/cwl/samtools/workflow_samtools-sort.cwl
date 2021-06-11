class: Workflow
cwlVersion: v1.0
inputs:
  samtools_bam: File
  SampleName: string
  thread: int
  results_path: Directory
outputs:
  samtools_sort_bam:
    outputSource:
    - samtools/samtools_sort_bam
    - sort/samtools_sort_bam
    type: File
steps:
  samtools_sort:
    in:
      samtools_bam: samtools_bam
      SampleName: SampleName
      thread: thread
      results_path: results_path
    out: '[samtools_sort_bam]'
    run: samtools/samtools-sort.cwl
