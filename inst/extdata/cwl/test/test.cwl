################################################################
##                        test-index                          ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[test - testthat]"
label: Last updated 05/20121

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [Rscript]

arguments:
  - prefix: --vanilla
    valueFrom: $(inputs.rscript_test)

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  rscript_test:
    type: File
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  test1:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/targets_test.txt
