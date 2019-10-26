################################################################
##                      GATK_mergebams.cwl                    ##
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
    valueFrom: MergeBamAlignment 
    position: 2
  
  - prefix: --ALIGNED_BAM 
    valueFrom: $(inputs.bwa_sam)
    position: 3
  
  - prefix: --UNMAPPED_BAM
    valueFrom: $(inputs.ubam)
    position: 4
  
  - prefix: --REFERENCE_SEQUENCE
    valueFrom: $(inputs.data_path.path)/$(inputs.ref_name)
    position: 5

  - prefix: --OUTPUT
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_merge.bam
    position: 6

  - prefix: --VALIDATION_STRINGENCY
    valueFrom: SILENT
    position: 7

  - prefix: --EXPECTED_ORIENTATIONS
    valueFrom: FR
    position: 8

  - prefix: --ATTRIBUTES_TO_RETAIN
    valueFrom: X0
    position: 9

  - prefix: --SORT_ORDER
    valueFrom: queryname
    position: 10

  - prefix: --IS_BISULFITE_SEQUENCE
    valueFrom: 'false'
    position: 11

  - prefix: --ADD_MATE_CIGAR
    valueFrom: 'true'
    position: 14

  - prefix: --MAX_INSERTIONS_OR_DELETIONS
    valueFrom: '-1'
    position: 15

  - prefix: --UNMAPPED_READ_STRATEGY
    valueFrom: COPY_TO_TAG
    position: 16

  - prefix: --ALIGNER_PROPER_PAIR_FLAGS
    valueFrom: 'true'
    position: 17

  - prefix: --UNMAP_CONTAMINANT_READS
    valueFrom: 'true'
    position: 18

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

    ubam:
        type: File

    bwa_sam:
        type: File

outputs:
  merge_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_merge.bam

