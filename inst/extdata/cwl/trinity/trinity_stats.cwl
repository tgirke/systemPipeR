################################################################
##                  Trinity-Paired_end.cwl                     ##
################################################################

cwlVersion: v1.1
class: CommandLineTool
label: Last updated 07/2021
hints:
  SoftwareRequirement:
    packages:
    - package: trinity-rnaseq
      version: [ 2.11.0 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [$TRINITY_HOME/util/TrinityStats.pl]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]


################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  reference:
    type: File
    inputBinding:
      position: 1
  results_path:
    type: Directory

stdout: $(input.results_path)/assembly_Statistic.txt

outputs:
  trinity_stats:
    type: stdout
