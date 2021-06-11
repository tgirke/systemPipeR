class: Workflow
cwlVersion: v1.0
inputs:
  idx_basedir: Directory
  idx_basename: string
  kallisto_idx: string
  fq1: File
  fq2: File
  thread: int
  SampleName: string
  results_path: Directory
  bootstrap_samples: int
  fragment_length: int
  std_dev: int
outputs:
  run_info:
    outputSource:
    - kallisto/run_info
    - quant/run_info
    type: File
  kallisto:
    outputSource:
    - kallisto/kallisto
    - quant/kallisto
    type: File
  kallisto_file:
    outputSource:
    - kallisto/kallisto_file
    - quant/kallisto_file
    type: File
steps:
  kallisto_quant:
    in:
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      kallisto_idx: kallisto_idx
      fq1: fq1
      fq2: fq2
      thread: thread
      SampleName: SampleName
      results_path: results_path
      bootstrap_samples: bootstrap_samples
      fragment_length: fragment_length
      std_dev: std_dev
    out: '[run_info, kallisto, kallisto_file]'
    run: kallisto/kallisto.cwl
