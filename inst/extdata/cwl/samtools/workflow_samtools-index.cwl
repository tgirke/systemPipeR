class: Workflow
cwlVersion: v1.0
inputs:
  samtools_sort_bam: File
  SampleName: string
  results_path: Directory
outputs:
  samtools_index:
    outputSource:
    - samtools/samtools_index
    - index/samtools_index
    - -b/samtools_index
    type: File
steps:
  samtools_index_-b:
    in:
      samtools_sort_bam: samtools_sort_bam
      SampleName: SampleName
      results_path: results_path
    out: '[samtools_index]'
    run: samtools/samtools-index.cwl
