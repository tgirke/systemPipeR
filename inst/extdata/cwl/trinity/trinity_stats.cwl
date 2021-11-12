################################################################
##                  Trinity-stats.cwl                         ##
################################################################

cwlVersion: v1.2
class: CommandLineTool
label: Last updated 11/2021
hints:
  SoftwareRequirement:
    packages:
    - package: trinity-rnaseq
      version: [ 2.12.0 ]

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

stdout: $(input.results_path)/N50_assembly_Statistic.txt

outputs:
  trinity_stats:
    type: stdout
