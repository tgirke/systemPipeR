class: Workflow
cwlVersion: v1.0
inputs:
  samtools_sam: File
  SampleName: string
  results_path: Directory
outputs:
  samtools_bam:
    outputSource:
    - samtools/samtools_bam
    - view/samtools_bam
    - -bS/samtools_bam
    type: File
steps:
  samtools_view_-bS:
    in:
      samtools_sam: samtools_sam
      SampleName: SampleName
      results_path: results_path
    out: '[samtools_bam]'
    run: samtools/samtools-view.cwl
