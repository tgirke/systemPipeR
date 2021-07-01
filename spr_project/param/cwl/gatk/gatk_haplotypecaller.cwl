################################################################
##                 GATK_haplotypecaller.cwl                   ##
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
    valueFrom: HaplotypeCaller 
    position: 2
  
  - prefix: -I
    valueFrom: $(inputs.fixed_bam)
    position: 3
  
  - prefix: -R
    valueFrom: $(inputs.data_path.path)/$(inputs.ref_name)
    position: 4

  - prefix: -O
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).g.vcf.gz
    position: 5

  - prefix: -ERC
    valueFrom: GVCF
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

  fixed_bam:
    type: File

outputs:
  gvcf:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).g.vcf.gz

  gvcf_index:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).g.vcf.gz.tbi
