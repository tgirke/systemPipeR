################################################################
##                 GATK_markduplicates.cwl                    ##
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
    valueFrom: MarkDuplicates 
    position: 2
  
  - prefix: -I
    valueFrom: $(inputs.sort_bam)
    position: 3
  
  - prefix: -O
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_mark.bam
    position: 4
  
  - prefix: -M
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).duplicate_metrics
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
 
  sort_bam:
    type: File

outputs:
  mark_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_mark.bam

  duplicate_metrics:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).duplicate_metrics
