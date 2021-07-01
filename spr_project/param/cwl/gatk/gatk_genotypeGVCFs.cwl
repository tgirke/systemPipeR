################################################################
##                  GATK_genotypeGVCFs.cwl                    ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
label: ""
doc: |
    written by Le Zhang
        11/2019
        IMPORTANT: MUST add this flag when running with cwltool: --relax-path-checks
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
    valueFrom: GenotypeGVCFs 
    position: 2
  
  - prefix: -V 
    valueFrom: $(inputs.gvcf_db_folder)
    position: 3
  
  - prefix: -O
    valueFrom: $(inputs.results_path.path)/samples.vcf.gz
    position: 4
  
  - prefix: -R
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

  gvcf_db_folder:
    type: string

outputs:
  vcfs:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/samples.vcf.gz

  vcfs_index:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/samples.vcf.gz.tbi
