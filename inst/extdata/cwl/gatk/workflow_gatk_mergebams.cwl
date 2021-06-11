class: Workflow
cwlVersion: v1.0
inputs:
  gatk_java_options: string
  SampleName: string
  data_path: Directory
  results_path: Directory
  ref_name: string
  ubam: File
  bwa_sam: File
outputs:
  merge_bam:
    outputSource: gatk/merge_bam
    type: File
steps:
  gatk:
    in:
      gatk_java_options: gatk_java_options
      SampleName: SampleName
      data_path: data_path
      results_path: results_path
      ref_name: ref_name
      ubam: ubam
      bwa_sam: bwa_sam
    out: '[merge_bam]'
    run: gatk/gatk_mergebams.cwl
