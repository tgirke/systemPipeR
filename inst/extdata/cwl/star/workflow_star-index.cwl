class: Workflow
cwlVersion: v1.0
inputs:
  thread: int
  runMode: string
  outFileNamePrefix: string
  genomeFastaFiles: File
  genomeSAindexNbases: int
  sjdbGTFfile: File
  sjdbGTFfeatureExon: string
  sjdbGTFtagExonParentTranscript: string
  sjdbOverhang: int
  idx_basedir: Directory
  results_path: Directory
outputs:
  chrName_txt:
    outputSource: STAR/chrName_txt
    type: File
  chrLength_txt:
    outputSource: STAR/chrLength_txt
    type: File
  chrStart_txt:
    outputSource: STAR/chrStart_txt
    type: File
  chrNameLength_txt:
    outputSource: STAR/chrNameLength_txt
    type: File
  exonGeTrInfo_tab:
    outputSource: STAR/exonGeTrInfo_tab
    type: File
  geneInfo_tab:
    outputSource: STAR/geneInfo_tab
    type: File
  transcriptInfo_tab:
    outputSource: STAR/transcriptInfo_tab
    type: File
  exonInfo_tab:
    outputSource: STAR/exonInfo_tab
    type: File
  sjdbList_fromGTF_out_tab:
    outputSource: STAR/sjdbList_fromGTF_out_tab
    type: File
  sjdbInfo_txt:
    outputSource: STAR/sjdbInfo_txt
    type: File
  sjdbList_out_tab:
    outputSource: STAR/sjdbList_out_tab
    type: File
  genomeParameters_txt:
    outputSource: STAR/genomeParameters_txt
    type: File
  Genome:
    outputSource: STAR/Genome
    type: File
  SA:
    outputSource: STAR/SA
    type: File
  SAindex:
    outputSource: STAR/SAindex
    type: File
  idx_Log_out:
    outputSource: STAR/idx_Log_out
    type: File
steps:
  STAR:
    in:
      thread: thread
      runMode: runMode
      outFileNamePrefix: outFileNamePrefix
      genomeFastaFiles: genomeFastaFiles
      genomeSAindexNbases: genomeSAindexNbases
      sjdbGTFfile: sjdbGTFfile
      sjdbGTFfeatureExon: sjdbGTFfeatureExon
      sjdbGTFtagExonParentTranscript: sjdbGTFtagExonParentTranscript
      sjdbOverhang: sjdbOverhang
      idx_basedir: idx_basedir
      results_path: results_path
    out: '[chrName_txt, chrLength_txt, chrStart_txt, chrNameLength_txt, exonGeTrInfo_tab,
      geneInfo_tab, transcriptInfo_tab, exonInfo_tab, sjdbList_fromGTF_out_tab, sjdbInfo_txt,
      sjdbList_out_tab, genomeParameters_txt, Genome, SA, SAindex, idx_Log_out]'
    run: star/star-index.cwl
