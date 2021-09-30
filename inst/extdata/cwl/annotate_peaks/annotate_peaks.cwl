################################################################
##                     annotatepeaks.cwl                      ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "Parameter file for 'mergeBamByFactor' function: R-based function that read from input files and write to output files. This is a dummy file to support within the systemPipeChIPseq workflow."
label: Last updated 10/2019

################################################################
##           baseCommand and arguments definitions            ##
################################################################

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  fq1:
    label: "Path to Input files"
    type: File
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  annotated:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName.basename)_peaks.annotated.xls


