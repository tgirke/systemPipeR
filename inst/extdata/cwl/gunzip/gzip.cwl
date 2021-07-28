################################################################
##                             gzip                           ##
################################################################

cwlVersion: v1.0
class: CommandLineTool

################################################################
##            baseCommand and arguments definitions           ##
################################################################

baseCommand: [ "gzip" ]

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.results_path)

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
  gzip_file:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName).$(inputs.ext)

###########
## Notes ##
###########

## If the template its used in bash script with the "cwl-runner", run: 
# "cwl-runner --outdir <path>/ gzip.cwl gzip.yml"

