class: Workflow
cwlVersion: v1.0
inputs:
  fq1: File
  SampleName: string
  results_path: Directory
  thread: int
  minlen: string
  leading: string
  trailing: string
  slidingwindow: string
outputs:
  trimmomatic:
    outputSource:
    - trimmomatic/trimmomatic
    - SE/trimmomatic
    type: File
steps:
  trimmomatic_SE:
    in:
      fq1: fq1
      SampleName: SampleName
      results_path: results_path
      thread: thread
      minlen: minlen
      leading: leading
      trailing: trailing
      slidingwindow: slidingwindow
    out: '[trimmomatic]'
    run: trimmomatic/trimmomatic-se.cwl
