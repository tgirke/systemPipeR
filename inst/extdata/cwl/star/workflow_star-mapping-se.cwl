class: Workflow
cwlVersion: v1.0
inputs:
  thread: int
  runMode: string
  star_idx_basedir: Directory
  Genome: File
  SA: File
  SAindex: File
  outSAMtype: string[]
  quantMode: string[]
  twopassMode: string
  readFilesCommand: string
  fq1: File
  SampleName: string
  results_path: Directory
outputs:
  Aligned_out_sam:
    outputSource: STAR/Aligned_out_sam
    type: File
  SJ_out.tab:
    outputSource: STAR/SJ_out.tab
    type: File
  Log_progress_out:
    outputSource: STAR/Log_progress_out
    type: File
  Log_final_out:
    outputSource: STAR/Log_final_out
    type: File
  Log_out:
    outputSource: STAR/Log_out
    type: File
  Aligned_toTranscriptome_out_bam:
    outputSource: STAR/Aligned_toTranscriptome_out_bam
    type: File
  ReadsPerGene_out_tab:
    outputSource: STAR/ReadsPerGene_out_tab
    type: File
steps:
  STAR:
    in:
      thread: thread
      runMode: runMode
      star_idx_basedir: star_idx_basedir
      Genome: Genome
      SA: SA
      SAindex: SAindex
      outSAMtype: outSAMtype
      quantMode: quantMode
      twopassMode: twopassMode
      readFilesCommand: readFilesCommand
      fq1: fq1
      SampleName: SampleName
      results_path: results_path
    out: '[Aligned_out_sam, SJ_out.tab, Log_progress_out, Log_final_out, Log_out,
      Aligned_toTranscriptome_out_bam, ReadsPerGene_out_tab]'
    run: star/star-mapping-se.cwl
