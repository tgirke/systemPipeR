class: Workflow
cwlVersion: v1.0
inputs:
  idx_basename: string
  idx_basedir: Directory
outputs:
  bwa.index_files1:
    outputSource:
    - bwa/bwa.index_files1
    - index/bwa.index_files1
    - -a/bwa.index_files1
    - bwtsw/bwa.index_files1
    type: File
  bwa.index_files2:
    outputSource:
    - bwa/bwa.index_files2
    - index/bwa.index_files2
    - -a/bwa.index_files2
    - bwtsw/bwa.index_files2
    type: File
  bwa.index_files3:
    outputSource:
    - bwa/bwa.index_files3
    - index/bwa.index_files3
    - -a/bwa.index_files3
    - bwtsw/bwa.index_files3
    type: File
  bwa.index_files4:
    outputSource:
    - bwa/bwa.index_files4
    - index/bwa.index_files4
    - -a/bwa.index_files4
    - bwtsw/bwa.index_files4
    type: File
  bwa.index_files5:
    outputSource:
    - bwa/bwa.index_files5
    - index/bwa.index_files5
    - -a/bwa.index_files5
    - bwtsw/bwa.index_files5
    type: File
steps:
  bwa_index_-a_bwtsw:
    in:
      idx_basename: idx_basename
      idx_basedir: idx_basedir
    out: '[bwa.index_files1, bwa.index_files2, bwa.index_files3, bwa.index_files4,
      bwa.index_files5]'
    run: bwa/bwa-index.cwl
