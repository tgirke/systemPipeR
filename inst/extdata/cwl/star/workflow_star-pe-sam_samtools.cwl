######################################
## Workflow_star-Paired-end-Reads ##
######################################

class: Workflow
cwlVersion: v1.0

################################################################
##              Inputs and Outputs Settings                   ##
################################################################

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

inputs:
  ## for mapping
  Genome: File
  SA: File
  SAindex: File
  fq1: File
  fq2: File
  star_idx_basedir: Directory
  SampleName: string
  thread: int
  results_path: Directory
  runMode: string
  outSAMtype: string[]
  quantMode: string[]
  twopassMode: string
  readFilesCommand: string

outputs:
  ## for mapping
  Aligned_out_sam: 
    outputSource: star/Aligned_out_sam
    type: File
  SJ_out_tab:
    outputSource: star/SJ_out_tab
    type: File
  Log_progress_out:
    outputSource: star/Log_progress_out
    type: File
  Log_final_out:
    outputSource: star/Log_final_out
    type: File
  Log_out:
    outputSource: star/Log_out
    type: File
  Aligned_toTranscriptome_out_bam:
    outputSource: star/Aligned_toTranscriptome_out_bam
    type: File
  ReadsPerGene_out_tab:
    outputSource: star/ReadsPerGene_out_tab
    type: File
  ## for samtools
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
  star:
    in:
      fq1: fq1
      fq2: fq2
      star_idx_basedir: star_idx_basedir
      SampleName: SampleName
      thread: thread
      results_path: results_path
      runMode: runMode
      outSAMtype: outSAMtype
      quantMode: quantMode
      twopassMode: twopassMode
      readFilesCommand: readFilesCommand
      Genome: Genome
      SA: SA
      SAindex: SAindex
    out: [Aligned_out_sam, SJ_out_tab, Log_progress_out, Log_final_out, Log_out, Aligned_toTranscriptome_out_bam, ReadsPerGene_out_tab]
    run: star/star-pe/star-mapping-pe.cwl 

  samtools-view:
    in:
      SampleName: SampleName
      samtools_sam: star/Aligned_out_sam
      results_path: results_path
    out: [samtools_bam]
    run: samtools/samtools-view.cwl

  samtools-sort:
    in:
      SampleName: SampleName
      samtools_bam: samtools-view/samtools_bam
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