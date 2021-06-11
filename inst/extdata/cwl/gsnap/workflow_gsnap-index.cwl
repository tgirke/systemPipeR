class: Workflow
cwlVersion: v1.0
inputs:
  rscript_gsnap_idx: File
  idx_basedir: Directory
  idx_basename: string
outputs:
  gsnap_index1:
    outputSource: Rscript/gsnap_index1
    type: File
  gsnap_index2:
    outputSource: Rscript/gsnap_index2
    type: File
steps:
  Rscript:
    in:
      rscript_gsnap_idx: rscript_gsnap_idx
      idx_basedir: idx_basedir
      idx_basename: idx_basename
    out: '[gsnap_index1, gsnap_index2]'
    run: gsnap/gsnap-index.cwl
