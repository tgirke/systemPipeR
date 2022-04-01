################################
## preprocessReads-SE_end.yml ##
################################

cwlVersion: v1.2
class: CommandLineTool
doc: "Run custom read preprocessing functions"
label: Last updated 11/2021

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [Rscript]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: --vanilla
    valueFrom: $(inputs.rscript_preprocessReads)
  - prefix: --outfile1
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).fastq_trim.gz
    
################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  rscript_preprocessReads:
    type: File
  fq1:
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    type: File
    inputBinding:
      prefix: --FileName1
  Fct:
    type: string
    inputBinding:
      prefix: --Fct
  batchsize:
    type: int
    inputBinding:
      prefix: --batchsize
  overwrite:
    type: boolean
    inputBinding:
      prefix: --overwrite
  compress:
    type: boolean
    inputBinding:
      prefix: --compress
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  preprocessReads_se:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).fastq_trim.gz

