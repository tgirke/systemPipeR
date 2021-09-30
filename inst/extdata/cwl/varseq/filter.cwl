################################################################
##                        Var-seq_Filter                      ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "Parameter file for 'readVcf' function: R-based function that read from input files and write to output files. This is a dummy file to support within the systemPipeR workflow."
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

outputs:
  filter:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.file1.basename).filter.vcf.bgz
