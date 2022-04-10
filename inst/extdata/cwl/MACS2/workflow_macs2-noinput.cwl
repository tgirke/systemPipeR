class: Workflow
cwlVersion: v1.0
inputs:
  SampleName: string
  results_path: Directory
  fq1: File
  format: string
  gsize: int
  qvalue: int
  
outputs:
  peaks_xls:
    outputSource:
    - macs2/peaks_xls
    - callpeak/peaks_xls
    type: File
  control_lambda:
    outputSource:
    - macs2/control_lambda
    - callpeak/control_lambda
    type: File
  narrow_peak:
    outputSource:
    - macs2/narrow_peak
    - callpeak/narrow_peak
    type: File
  summits_bed:
    outputSource:
    - macs2/summits_bed
    - callpeak/summits_bed
    type: File
  macs2_pileup:
    outputSource:
    - macs2/macs2_pileup
    - callpeak/macs2_pileup
    type: File
steps:
  macs2_callpeak:
    in:
      SampleName: SampleName
      results_path: results_path
      fq1: fq1
    out: '[peaks_xls, control_lambda, narrow_peak, summits_bed, macs2_pileup]'
    run: MACS2/macs2-noinput.cwl
