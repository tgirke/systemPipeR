class: Workflow
cwlVersion: v1.0
inputs:
  gatk_java_options: string
  SampleName: string
  data_path: Directory
  results_path: Directory
  fq1: File
  fq2: File?
outputs:
  ubam:
    outputSource: gatk/ubam
    type: File
steps:
  gatk:
    in:
      gatk_java_options: gatk_java_options
      SampleName: SampleName
      data_path: data_path
      results_path: results_path
      fq1: fq1
      fq2: fq2
    out: '[ubam]'
    run: gatk/gatk_fastq2ubam.cwl
