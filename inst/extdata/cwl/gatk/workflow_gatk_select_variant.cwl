class: Workflow
cwlVersion: v1.0
inputs:
  gatk_java_options: string
  SampleName: string
  data_path: Directory
  results_path: Directory
  ref_name: string
  cohort_filtered_vcf: File
outputs:
  vcf:
    outputSource: bcftools/vcf
    type: File
steps:
  bcftools:
    in:
      gatk_java_options: gatk_java_options
      SampleName: SampleName
      data_path: data_path
      results_path: results_path
      ref_name: ref_name
      cohort_filtered_vcf: cohort_filtered_vcf
    out: '[vcf]'
    run: gatk/gatk_select_variant.cwl
