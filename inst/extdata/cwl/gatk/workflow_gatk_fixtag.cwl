class: Workflow
cwlVersion: v1.0
inputs:
  gatk_java_options: string
  SampleName: string
  data_path: Directory
  results_path: Directory
  ref_name: string
  mark_bam: File
outputs:
  fixtag_bam:
    outputSource: gatk/fixtag_bam
    type: File
  fixtag_bam_index:
    outputSource: gatk/fixtag_bam_index
    type: File
steps:
  gatk:
    in:
      gatk_java_options: gatk_java_options
      SampleName: SampleName
      data_path: data_path
      results_path: results_path
      ref_name: ref_name
      mark_bam: mark_bam
    out: '[fixtag_bam, fixtag_bam_index]'
    run: gatk/gatk_fixtag.cwl
