class: Workflow
cwlVersion: v1.0
inputs:
  idx_basename: string
  idx_basedir: Directory
outputs:
  index_files1:
    outputSource: bowtie2-build/index_files1
    type: File
  index_files2:
    outputSource: bowtie2-build/index_files2
    type: File
  index_files3:
    outputSource: bowtie2-build/index_files3
    type: File
  index_files4:
    outputSource: bowtie2-build/index_files4
    type: File
  index_files5:
    outputSource: bowtie2-build/index_files5
    type: File
  index_files6:
    outputSource: bowtie2-build/index_files6
    type: File
steps:
  bowtie2-build:
    in:
      idx_basename: idx_basename
      idx_basedir: idx_basedir
    out: '[index_files1, index_files2, index_files3, index_files4, index_files5, index_files6]'
    run: bowtie2/bowtie2-index.cwl
