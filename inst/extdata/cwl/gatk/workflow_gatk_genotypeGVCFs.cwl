class: Workflow
cwlVersion: v1.0
inputs:
  gatk_java_options: string
  SampleName: string
  data_path: Directory
  results_path: Directory
  ref_name: string
  gvcf_db_folder: string
outputs:
  vcfs:
    outputSource: gatk/vcfs
    type: File
  vcfs_index:
    outputSource: gatk/vcfs_index
    type: File
steps:
  gatk:
    in:
      gatk_java_options: gatk_java_options
      SampleName: SampleName
      data_path: data_path
      results_path: results_path
      ref_name: ref_name
      gvcf_db_folder: gvcf_db_folder
    out: '[vcfs, vcfs_index]'
    run: gatk/gatk_genotypeGVCFs.cwl
