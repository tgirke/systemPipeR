################################################################
##                           gunzip                           ##
################################################################

cwlVersion: v1.0
class: CommandLineTool

################################################################
##            baseCommand and arguments definitions           ##
################################################################

baseCommand: [ "gunzip" ]

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.results_path)

arguments:
  - prefix: -c

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  file:
    type: File
    inputBinding:
      position: 1
  SampleName:
    label: "Filename to write output to"
    type: string
  ext:
    label: "Filename extension"
    type: string
  results_path:
    label: "Directory to write output to"
    type: Directory

outputs:
    gunzip_file:
        type: stdout
    
stdout: $(inputs.results_path.basename)/$(inputs.SampleName).csv


###########
## Notes ##
###########

## If the template its used in bash script with the "cwl-runner", run: 
# "cwl-runner --outdir <path> gunzip.cwl gunzip.yml"

