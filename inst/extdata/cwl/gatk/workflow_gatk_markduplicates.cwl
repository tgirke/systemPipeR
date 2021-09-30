class: Workflow
cwlVersion: v1.0
inputs:
  gatk_java_options: string
  SampleName: string
  data_path: Directory
  results_path: Directory
  sort_bam: File
outputs:
  mark_bam:
    outputSource: gatk/mark_bam
    type: File
  duplicate_metrics:
    outputSource: gatk/duplicate_metrics
    type: File
steps:
  gatk:
    in:
      gatk_java_options: gatk_java_options
      SampleName: SampleName
      data_path: data_path
      results_path: results_path
      sort_bam: sort_bam
    out: '[mark_bam, duplicate_metrics]'
    run: gatk/gatk_markduplicates.cwl
