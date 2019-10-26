################################################################
##                        GATK_sort.cwl                       ##
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

baseCommand: [ gatk ]

arguments:
  - prefix: --java-options
    valueFrom: $(inputs.gatk_java_options)
    position: 1
  
  - prefix: 
    valueFrom: SortSam 
    position: 2 
  
  - prefix: --INPUT
    valueFrom: $(inputs.merge_bam)
    position: 3
  
  - prefix: --OUTPUT
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_sort.bam
    position: 4
  
  - prefix: --SORT_ORDER
    valueFrom: coordinate
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

  merge_bam:
    type: File

outputs:
  sort_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_sort.bam
