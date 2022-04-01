#######################################
## Workflow_Bowtie2-Paired-end-Reads ##
#######################################

class: Workflow
cwlVersion: v1.0

################################################################
##              Inputs and Outputs Settings                   ##
################################################################

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

inputs:
  fq1: File
  fq2: File
  idx_basedir: Directory
  idx_basename: string
  SampleName: string
  thread: int
  results_path: Directory

outputs:
  Bowtie2:
    outputSource: bowtie2/bowtie2_sam
    type: File
  samtools-view:
    outputSource: samtools-view/samtools_bam
    type: File
  samtools-sort:
    outputSource: samtools-sort/samtools_sort_bam
    type: File
  samtools-index:
    outputSource: samtools-index/samtools_index
    type: File

################################################################
##                Workflow Steps Definitions                  ##
################################################################

steps:
  bowtie2:
    in:
      fq1: fq1
      fq2: fq2
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      SampleName: SampleName
      thread: thread
      results_path: results_path
    out: [bowtie2_sam]
    run: bowtie2/bowtie2-mapping-pe.cwl
    
  samtools-view:
    in:
      samtools_sam: bowtie2/bowtie2_sam
      SampleName: SampleName
      results_path: results_path
    out: [samtools_bam]
    run: samtools/samtools-view.cwl

  samtools-sort:
    in:
      samtools_bam: samtools-view/samtools_bam
      SampleName: SampleName
      thread: thread
      results_path: results_path
    out: [samtools_sort_bam]
    run: samtools/samtools-sort.cwl
   
  samtools-index:
    in:
      samtools_sort_bam: samtools-sort/samtools_sort_bam
      SampleName: SampleName
      results_path: results_path
    out: [samtools_index]
    run: samtools/samtools-index.cwl
    
