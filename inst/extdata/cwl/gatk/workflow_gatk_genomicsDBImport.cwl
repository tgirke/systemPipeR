class: Workflow
cwlVersion: v1.0
inputs:
  gatk_java_options: string
  SampleName: string
  data_path: Directory
  results_path: Directory
  ref_name: string
  scripts_path: Directory
outputs:
  gvcfs:
    outputSource: bash/gvcfs
    type: Directory?
steps:
  bash:
    in:
      gatk_java_options: gatk_java_options
      SampleName: SampleName
      data_path: data_path
      results_path: results_path
      ref_name: ref_name
      scripts_path: scripts_path
    out: '[gvcfs]'
    run: gatk/gatk_genomicsDBImport.cwl
