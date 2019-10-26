################################################################
##                GATK_genomicsDBImport.cwl                   ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
label: ""
doc: |
    written by Le Zhang
        11/2019
hints:
  SoftwareRequirement:
    packages:
    - package: gatk
      version: [ 4.1.1.0 ]
      
################################################################
##           baseCommand and arguments definitions            ##
################################################################

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

baseCommand: [ bash ]

arguments:
  - prefix:
    valueFrom: $(inputs.scripts_path.path)/gatk_genomicsDBImport.sh
    position: 1

  - prefix:
    valueFrom: $(inputs.data_path.path)/$(inputs.ref_name)
    position: 2

  - prefix:
    valueFrom: $(inputs.results_path.path)
    position: 3

  - prefix:
    valueFrom: $(inputs.gatk_java_options)
    position: 4

  - prefix:
    valueFrom: $(inputs.results_path.path)/gvcfs
    position: 5

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  gatk_java_options:
    type: string

  SampleName:
    type: string

  data_path:
    type: Directory

  results_path:
    type: Directory

  ref_name:
    type: string

  scripts_path:
    type: Directory

outputs:
  gvcfs:
    type: Directory?
    outputBinding:
      glob: $(inputs.results_path.path)/gvcfs
