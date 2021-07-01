################################################################
##                      BUSCO.cwl                             ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
label: Last updated 07/2019
hints:
  SoftwareRequirement:
    packages:
    - package: busco
      version: [ 3.0.2 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [run_BUSCO.py]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: -f
  - prefix: -m
    valueFrom: 'transcriptome'

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  assembly:
    type: File
    inputBinding:
      prefix: -i
      position: 1
  lineage:
    type: Directory
    inputBinding:
      prefix: -l
      position: 2
  outputName:
    type: string
    inputBinding:
      prefix: -o
  cpu:
    type: int?
    inputBinding:
      prefix: --cpu
  results_path:
    type: Directory

outputs:
  busco_blast:
    type: File
    outputBinding:
      glob: run_BUSCO/blast_output
  busco_table:
    type: File
    outputBinding:
      glob: run_BUSCO/full_table_BUSCO.tsv
   busco_hmmer: 
     type: File 
     outputBinding:
      glob: run_BUSCO/hmmer_output
  busco_missing:
    type: File 
    outputBinding:
      glob: run_BUSCO/missing_busco_list_BUSCO.tsv
  busco_summary:
    type: File
    outputBinding:
      glob: run_BUSCO/short_summary_Trinity_BUSCO.txt
  busco_proteins:  
    type: File  
    outputBinding:  
      glob: run_BUSCO/translated_proteins

###########
## Notes ##
###########

## this version works with sysArgs2 (type: File)
## output directory name is hard-coded as run_outputName 
## output directory should be moved into results directory using:
## file.copy("BUSCO", "results", overwrite = recursive, recursive = TRUE)
