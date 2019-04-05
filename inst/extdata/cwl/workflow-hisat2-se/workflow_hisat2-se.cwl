################################################################
##                Workflow_Hisat2-Single_Read                 ##
################################################################

class: Workflow
cwlVersion: v1.0

################################################################
##              Inputs and Outputs Settings                   ##
################################################################

inputs:
  fq1: File
  hisat2_idx_basedir: Directory
  hisat2_idx_basename: string
  SampleName: string
  thread: int
  results_path: Directory

outputs:
  samtools-view:
    outputSource: samtools-view/samtools_bam
    type: File
  samtools-sort:
    outputSource: samtools-sort/samtools_sort_bam
    type: File
  samtools-index:
    outputSource: samtools-index/samtools-index
    type: File

################################################################
##                Workflow Steps Definitions                  ##
################################################################

steps:
  hisat2:
    in:
      fq1: fq1
      hisat2_idx_basedir: hisat2_idx_basedir
      hisat2_idx_basename: hisat2_idx_basename
      SampleName: SampleName
      thread: thread
      results_path: results_path
    out: [hisat2_sam]
    run: hisat2-mapping-se.cwl

  samtools-view:
    in:
      SampleName: SampleName
      sam: hisat2/hisat2_sam
      results_path: results_path
    out: [samtools_bam]
    run: samtools-view.cwl

  samtools-sort:
    in:
      SampleName: SampleName
      bam: samtools-view/samtools_bam
      thread: thread
      results_path: results_path
    out: [samtools_sort_bam]
    run: samtools-sort.cwl

  samtools-index:
    in:
      samtools_sort_bam: samtools-sort/samtools_sort_bam
      SampleName: SampleName
      results_path: results_path
    out: [samtools-index]
    run: samtools-index.cwl

