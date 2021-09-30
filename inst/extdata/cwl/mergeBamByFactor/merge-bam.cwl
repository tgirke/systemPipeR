################################################################
##                      Merge Bam files                      ##
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
  bam:
    label: "Path to bam files"
    type: File
  SampleName:
    label: "SampleName"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  merge-bam:
    type: File
    outputBinding:
      glob: $(inputs.bam)
