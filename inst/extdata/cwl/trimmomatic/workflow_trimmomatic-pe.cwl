class: Workflow
cwlVersion: v1.0
inputs:
  fq1: File
  fq2: File
  SampleName: string
  results_path: Directory
  thread: int
  minlen: string
  leading: string
  trailing: string
  slidingwindow: string
outputs:
  trimmomatic_1_paired:
    outputSource:
    - trimmomatic/trimmomatic_1_paired
    - PE/trimmomatic_1_paired
    type: File
  trimmomatic_1_unpaired:
    outputSource:
    - trimmomatic/trimmomatic_1_unpaired
    - PE/trimmomatic_1_unpaired
    type: File
  trimmomatic_2_paired:
    outputSource:
    - trimmomatic/trimmomatic_2_paired
    - PE/trimmomatic_2_paired
    type: File
  trimmomatic_2_unpaired:
    outputSource:
    - trimmomatic/trimmomatic_2_unpaired
    - PE/trimmomatic_2_unpaired
    type: File
steps:
  trimmomatic_PE:
    in:
      fq1: fq1
      fq2: fq2
      SampleName: SampleName
      results_path: results_path
      thread: thread
      minlen: minlen
      leading: leading
      trailing: trailing
      slidingwindow: slidingwindow
    out: '[trimmomatic_1_paired, trimmomatic_1_unpaired, trimmomatic_2_paired, trimmomatic_2_unpaired]'
    run: trimmomatic/trimmomatic-pe.cwl
