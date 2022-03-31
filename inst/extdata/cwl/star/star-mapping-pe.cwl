################################################################
##                     star-mapping-pe.cwl                    ##
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
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: --runThreadN
    valueFrom: $(inputs.thread)
  - prefix: --runMode
    valueFrom: $(inputs.runMode)
  - prefix: --genomeDir
    valueFrom: $(inputs.star_idx_basedir)
  - prefix: --outFileNamePrefix
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).
  - prefix: --outSAMtype
    valueFrom: $(inputs.outSAMtype)
  - prefix: --quantMode
    valueFrom: $(inputs.quantMode)
  - prefix: --twopassMode
    valueFrom: $(inputs.twopassMode)
  - prefix: --readFilesCommand
    valueFrom: $(inputs.readFilesCommand)
  - prefix: --readFilesIn
    valueFrom: $(inputs.fq1)
  - prefix: 
    valueFrom: $(inputs.fq2)

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  thread:
    label: NumberOfThreads
    type: int
  runMode:
    label: "Use alignReads for mapping"
    type: string
  star_idx_basedir:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory
  Genome:
    label: "Output by STAR generate_genome_index"
    type: File
  SA:
    label: "Output by STAR generate_genome_index"
    type: File
  SAindex:
    label: "Output by STAR generate_genome_index"
    type: File
  outSAMtype:
    label: "Output SAM or BAM options Unsorted SortedByCoordinate "
    type: string[]
  quantMode:
    label: "Output in transcript coordinates AND Counting number of reads per gene"
    type: string[]
  twopassMode:
    label: "Use Basic for automatic two pass mode, otherwise use None"
    type: string
  readFilesCommand:
    label: "Use zcat, gunzip -c or bunzip -c"
    type: string
  fq1:
    label: "use comma separated list for read1"
    type: File
  fq2:
    label: "use comma separated list for read2"
    type: File
  SampleName:
    label: "Prefix output filenames with this"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  Aligned_out_sam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).Aligned.out.sam
  SJ_out_tab:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).SJ.out.tab
  Log_progress_out:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).Log.progress.out
  Log_final_out:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).Log.final.out
  Log_out:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).Log.out
  Aligned_toTranscriptome_out_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).Aligned.toTranscriptome.out.bam
  ReadsPerGene_out_tab:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).ReadsPerGene.out.tab
#  _STARgenome:
#    type: Folder
#    outputBinding:
#      glob: $(inputs.results_path.path)/$(inputs.SampleName)._STARgenome
#  _STARpass1:
#    type: Folder
#    outputBinding:
#      glob: $(inputs.results_path.path)/$(inputs.SampleName)._STARpass1

###########
## Notes ##
###########

## If the template its used in bash script with the "cwltool" 
## output files directed to results folders will be duplicated in the
## starting working directory since it is the designated output directory for CWL
