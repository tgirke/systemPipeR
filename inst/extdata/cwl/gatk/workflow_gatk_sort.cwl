class: Workflow
cwlVersion: v1.0
inputs:
  gatk_java_options: string
  SampleName: string
  data_path: Directory
  results_path: Directory
  merge_bam: File
outputs:
  sort_bam:
    outputSource: gatk/sort_bam
    type: File
steps:
  gatk:
    in:
      gatk_java_options: gatk_java_options
      SampleName: SampleName
      data_path: data_path
      results_path: results_path
      merge_bam: merge_bam
    out: '[sort_bam]'
    run: gatk/gatk_sort.cwl
