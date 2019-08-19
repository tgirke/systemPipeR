################################################################
##                        trim-Paired_end                     ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "Parameter file for 'preprocessReads' function: R-based function that read from input files and write to output files. This is a dummy file to support within the systemPipeR workflow."
label: Last updated 08/2019

################################################################
##           baseCommand and arguments definitions            ##
################################################################

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  fq1:
    label: "Comma-separated list of files containing mate 1s to be aligned"
    type: File
  fq2:
    label: "Comma-separated list of files containing mate 2s to be aligned"
    type: File
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  trim-pe_1:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName)_1.fastq_trim.gz
  trim-pe_2:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName)_2.fastq_trim.gz

