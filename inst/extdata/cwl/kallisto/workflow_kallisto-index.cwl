class: Workflow
cwlVersion: v1.0
inputs:
  idx_basename: string
  idx_basedir: Directory
  kallisto_idx: string
outputs:
  kallisto_index:
    outputSource:
    - kallisto/kallisto_index
    - index/kallisto_index
    type: File
steps:
  kallisto_index:
    in:
      idx_basename: idx_basename
      idx_basedir: idx_basedir
      kallisto_idx: kallisto_idx
    out: '[kallisto_index]'
    run: kallisto/kallisto-index.cwl
