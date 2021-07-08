################################################################
##                      Var-seq_combine                      ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "Parameter file for 'combineVarReports' function: R-based function that read from input files and write to output files. This is a dummy file to support within the systemPipeR workflow."
label: Last updated 10/2019

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [ ]

arguments:
  - valueFrom: $(inputs.file1)
  - valueFrom: $(inputs.file2)

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  file1:
    type: File
  file2:
    type: File
  results_path:
    type: Directory
  SampleName:
    type: string
  ref_name:
    type: string

outputs:
  annotate:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.file1.basename)
