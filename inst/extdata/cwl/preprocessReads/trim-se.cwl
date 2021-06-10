################################################################
##                      trim-Single_Read                      ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "Parameter file for 'preprocessReads' function: R-based function that read from input files and write to output files. This is a dummy file to support within the systemPipeR workflow."
label: Last updated 09/2019

################################################################
##           baseCommand and arguments definitions            ##
################################################################

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  fq1:
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    type: File
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  trim-se:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).fastq_trim.gz

