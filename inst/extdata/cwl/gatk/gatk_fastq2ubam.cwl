################################################################
##                     GATK_fastq2ubam.cwl                    ##
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
    valueFrom: FastqToSam
    position: 2

  - prefix: --FASTQ
    valueFrom: $(inputs.fq1)
    position: 3

  - prefix: --FASTQ2
    valueFrom: $(inputs.fq2)
    position: 4

  - prefix: --SAMPLE_NAME
    valueFrom: $(inputs.SampleName)
    position: 6

  - prefix: -O
    valueFrom: $(inputs.results_path.basename)/$(inputs.SampleName).u.bam
    position: 7

  - prefix: --LIBRARY_NAME
    valueFrom: lib1
    position: 8

  - prefix: --READ_GROUP_NAME
    valueFrom: group1
    position: 9

  - prefix: --PLATFORM
    valueFrom: illumina
    position: 10

  - prefix: --PLATFORM_UNIT
    valueFrom: unit1
    position: 11

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

  fq1:
    type: File

  fq2:
    type: File?

outputs:
  ubam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).u.bam
