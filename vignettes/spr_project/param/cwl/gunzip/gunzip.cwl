################################################################
##                           gunzip                           ##
################################################################

cwlVersion: v1.0
class: CommandLineTool

################################################################
##            baseCommand and arguments definitions           ##
################################################################

baseCommand: [ "gunzip", "-c" ]

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.results_path)

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  fq1:
    type: File
    inputBinding:
      position: 1
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Directory to write output to"
    type: Directory

outputs:
    gunzip_file:
        type: stdout
    
stdout: $(inputs.results_path.basename)/$(inputs.SampleName).fastq


###########
## Notes ##
###########

## If the template its used in bash script with the "cwl-runner", run: 
# "cwl-runner --outdir <path>/ gunzip.cwl gunzip.yml"

