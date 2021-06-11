class: Workflow
cwlVersion: v1.0
inputs:
  SampleName: string
  samtools_sort_bam: File
  results_path: Directory
outputs:
  marked_duplicates_bam:
    outputSource: picard/marked_duplicates_bam
    type: File
  marked_dup_metrics_txt:
    outputSource: picard/marked_dup_metrics_txt
    type: File
steps:
  picard:
    in:
      SampleName: SampleName
      samtools_sort_bam: samtools_sort_bam
      results_path: results_path
    out: '[marked_duplicates_bam, marked_dup_metrics_txt]'
    run: picard/picard_markduplicates.cwl
