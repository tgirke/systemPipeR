class: Workflow
cwlVersion: v1.0
inputs:
  rscript_rsubread_idx: File
  idx_basedir: Directory
  idx_basename: string
outputs:
  rsubread_index1:
    outputSource: Rscript/rsubread_index1
    type: File
  rsubread_index2:
    outputSource: Rscript/rsubread_index2
    type: File
  rsubread_index3:
    outputSource: Rscript/rsubread_index3
    type: File
  rsubread_index4:
    outputSource: Rscript/rsubread_index4
    type: File
  rsubread_index5:
    outputSource: Rscript/rsubread_index5
    type: File
steps:
  Rscript:
    in:
      rscript_rsubread_idx: rscript_rsubread_idx
      idx_basedir: idx_basedir
      idx_basename: idx_basename
    out: '[rsubread_index1, rsubread_index2, rsubread_index3, rsubread_index4, rsubread_index5]'
    run: rsubread/rsubread-index.cwl
