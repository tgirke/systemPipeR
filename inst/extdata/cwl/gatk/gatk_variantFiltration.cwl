################################################################
##               GATK_variantFiltration.cwl                  ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
label: ""
doc: |
    written by Le Zhang
        11/2019

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [ bash ]

arguments:
  - prefix: 
    valueFrom: $(inputs.scripts_path.path)/gatk_variantFiltration.sh
    position: 1

  - prefix:
    valueFrom: $(inputs.gatk_java_options)
    position: 2
  
  - prefix:
    valueFrom: $(inputs.raw_vcf)
    position: 3
  
  - prefix:
    valueFrom: $(inputs.results_path.path)/samples_filter.vcf.gz
    position: 4
  
  - prefix:
    valueFrom: $(inputs.data_path.path)/$(inputs.ref_name)
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

  raw_vcf:
    type: File

  scripts_path:
    type: Directory

outputs:
  vcf:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/samples_filter.vcf

  vcf_index:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/samples_filter.vcf.idx
