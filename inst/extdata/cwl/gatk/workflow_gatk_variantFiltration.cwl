class: Workflow
cwlVersion: v1.0
inputs:
  gatk_java_options: string
  SampleName: string
  data_path: Directory
  results_path: Directory
  ref_name: string
  raw_vcf: File
  scripts_path: Directory
outputs:
  vcf:
    outputSource: bash/vcf
    type: File
  vcf_index:
    outputSource: bash/vcf_index
    type: File
steps:
  bash:
    in:
      gatk_java_options: gatk_java_options
      SampleName: SampleName
      data_path: data_path
      results_path: results_path
      ref_name: ref_name
      raw_vcf: raw_vcf
      scripts_path: scripts_path
    out: '[vcf, vcf_index]'
    run: gatk/gatk_variantFiltration.cwl
