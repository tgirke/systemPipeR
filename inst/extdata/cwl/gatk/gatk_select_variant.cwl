################################################################
##                  GATK_select_variant.cwl                   ##
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
    - package: bcftools
      version: [ 1.15 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################
     
requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

baseCommand: [ bcftools ]

arguments:
  - prefix: view
    valueFrom: 
    position: 1
  
  - prefix: -c1
    valueFrom: 
    position: 2
  
  - prefix: -Ov 
    valueFrom: 
    position: 3

  - prefix: -f
    valueFrom: PASS
    position: 4

  - prefix: -s
    valueFrom: $(inputs.SampleName)
    position: 5
  
  - prefix: -o
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_gatk.vcf
    position: 6

  - prefix: 
    valueFrom: $(inputs.cohort_filtered_vcf)
    position: 7

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

  cohort_filtered_vcf:
    type: File

outputs:
  vcf_raw:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_gatk.vcf
