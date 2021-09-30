class: Workflow
cwlVersion: v1.0
inputs:
  idx_basedir: Directory
  idx_basename: string
outputs:
  hisat2_index1:
    outputSource: hisat2-build/hisat2_index1
    type: File
  hisat2_index2:
    outputSource: hisat2-build/hisat2_index2
    type: File
  hisat2_index3:
    outputSource: hisat2-build/hisat2_index3
    type: File
  hisat2_index4:
    outputSource: hisat2-build/hisat2_index4
    type: File
  hisat2_index5:
    outputSource: hisat2-build/hisat2_index5
    type: File
  hisat2_index6:
    outputSource: hisat2-build/hisat2_index6
    type: File
  hisat2_index7:
    outputSource: hisat2-build/hisat2_index7
    type: File
  hisat2_index8:
    outputSource: hisat2-build/hisat2_index8
    type: File
steps:
  hisat2-build:
    in:
      idx_basedir: idx_basedir
      idx_basename: idx_basename
    out: '[hisat2_index1, hisat2_index2, hisat2_index3, hisat2_index4, hisat2_index5,
      hisat2_index6, hisat2_index7, hisat2_index8]'
    run: hisat2/hisat2-index.cwl
