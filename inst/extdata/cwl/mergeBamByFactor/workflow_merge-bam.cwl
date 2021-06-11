class: Workflow
cwlVersion: v1.0
inputs:
  bam: File
  SampleName: string
  results_path: Directory
outputs:
  merge-bam:
    outputSource: /merge-bam
    type: File
steps:
  merge-bam:
    in:
      bam: bam
      SampleName: SampleName
      results_path: results_path
    out: '[merge-bam]'
    run: mergeBamByFactor/merge-bam.cwl
