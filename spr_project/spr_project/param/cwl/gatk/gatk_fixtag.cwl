################################################################
##                      GATK_fixtag.cwl                       ##
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

baseCommand: [ gatk ]

arguments:
  - prefix: --java-options
    valueFrom: $(inputs.gatk_java_options)
    position: 1

  - prefix:  
    valueFrom: SetNmMdAndUqTags
    position: 2

  - prefix: --INPUT
    valueFrom: $(inputs.mark_bam)
    position: 3
  
  - prefix: --OUTPUT
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_fixed.bam
    position: 4
  
  - prefix: --REFERENCE_SEQUENCE
    valueFrom: $(inputs.data_path.path)/$(inputs.ref_name)
    position: 5
  
  - prefix: --CREATE_INDEX
    valueFrom: 'true'
    position: 6

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

  mark_bam:
    type: File

outputs:
  fixtag_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_fixed.bam

  fixtag_bam_index:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_fixed.bai
