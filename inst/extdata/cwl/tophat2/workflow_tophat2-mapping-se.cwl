class: Workflow
cwlVersion: v1.0
inputs:
  idx_basedir: Directory
  idx_basename: string
  fq1: File
  SampleName: string
  thread: int
  results_path: Directory
outputs:
  accepted_hits_bam:
    outputSource: tophat2/accepted_hits_bam
    type: File
  junctions_bed:
    outputSource: tophat2/junctions_bed
    type: File
  insertions_bed:
    outputSource: tophat2/insertions_bed
    type: File
  deletions_bed:
    outputSource: tophat2/deletions_bed
    type: File
  align_summary:
    outputSource: tophat2/align_summary
    type: File
steps:
  tophat2:
    in:
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      fq1: fq1
      SampleName: SampleName
      thread: thread
      results_path: results_path
    out: '[accepted_hits_bam, junctions_bed, insertions_bed, deletions_bed, align_summary]'
    run: tophat2/tophat2-mapping-se.cwl
