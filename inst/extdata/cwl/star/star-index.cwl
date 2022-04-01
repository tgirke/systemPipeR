################################################################
##                        STAR-Index                        ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[STAR](https://github.com/alexdobin/STAR): Spliced Transcripts Alignment to a Reference"
label: Last updated 03/2022
hints:
  SoftwareRequirement:
    packages:
    - package: STAR
      version: [ 2.7.10a ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [STAR]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.idx_basedir) ]

arguments:
  - prefix: --runThreadN
    valueFrom: $(inputs.thread)
  - prefix: --runMode
    valueFrom: $(inputs.runMode)
  - prefix: --genomeDir
    valueFrom: $(inputs.idx_basedir)
  - prefix: --outFileNamePrefix
    valueFrom: $(inputs.outFileNamePrefix)
  - prefix: --genomeFastaFiles
    valueFrom: $(inputs.genomeFastaFiles)
  - prefix: --genomeSAindexNbases
    valueFrom: $(inputs.genomeSAindexNbases)
  - prefix: --sjdbGTFfile
    valueFrom: $(inputs.sjdbGTFfile)
  - prefix: --sjdbGTFfeatureExon
    valueFrom: $(inputs.sjdbGTFfeatureExon)
  - prefix: --sjdbGTFtagExonParentTranscript
    valueFrom: $(inputs.sjdbGTFtagExonParentTranscript)
  - prefix: --sjdbOverhang
    valueFrom: $(inputs.sjdbOverhang)

    
################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  thread:
    label: NumberOfThreads
    type: int
  runMode:
    label: genomeGenerate
    type: string
  outFileNamePrefix:
    label: /path/to/output/dir/prefix
    type: string
  genomeFastaFiles:
    label: /path/to/genome/fasta1 /path/to/genome/fasta2 ...
    type: File
  genomeSAindexNbases:
    label: min(14, log2(GenomeLength)/2 - 1)
    type: int
  sjdbGTFfile:
    label: /path/to/annotations.gtf
    type: File
  sjdbGTFfeatureExon:
    label: generally exon
    type: string
  sjdbGTFtagExonParentTranscript:
    label: Parent
    type: string
  sjdbOverhang:
    label: ReadLength-1
    type: int
  idx_basedir:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory

  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  chrName_txt:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/chrName.txt

  chrLength_txt:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/chrLength.txt

  chrStart_txt:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/chrStart.txt

  chrNameLength_txt:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/chrNameLength.txt

  exonGeTrInfo_tab:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/exonGeTrInfo.tab

  geneInfo_tab:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/geneInfo.tab

  transcriptInfo_tab:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/transcriptInfo.tab

  exonInfo_tab:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/exonInfo.tab

  sjdbList_fromGTF_out_tab:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/sjdbList.fromGTF.out.tab

  sjdbInfo_txt:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/sjdbInfo.txt

  sjdbList_out_tab:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/sjdbList.out.tab

  genomeParameters_txt:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/genomeParameters.txt

  Genome:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/Genome

  SA:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/SA

  SAindex:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/SAindex

  idx_Log_out:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/Log.out

