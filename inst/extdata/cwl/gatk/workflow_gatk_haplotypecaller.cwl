class: Workflow
cwlVersion: v1.0
inputs:
  gatk_java_options: string
  SampleName: string
  data_path: Directory
  results_path: Directory
  ref_name: string
  fixed_bam: File
outputs:
  gvcf:
    outputSource: gatk/gvcf
    type: File
  gvcf_index:
    outputSource: gatk/gvcf_index
    type: File
steps:
  gatk:
    in:
      gatk_java_options: gatk_java_options
      SampleName: SampleName
      data_path: data_path
      results_path: results_path
      ref_name: ref_name
      fixed_bam: fixed_bam
    out: '[gvcf, gvcf_index]'
    run: gatk/gatk_haplotypecaller.cwl
