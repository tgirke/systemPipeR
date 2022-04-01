################################################################
##                           samtools.cwl                    ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: ""
label: Last updated 03/2022
hints:
  SoftwareRequirement:
    packages:
      - package: samtools/1.14

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [ samtools ]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: rmdup
    valueFrom: -S
    position: 1

  - prefix: 
    valueFrom: $(inputs.bam.path) 
    position: 2
 
  - prefix: 
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_mark.bam
    position: 3

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  bam:
    type: File

  SampleName:
    label: "Prefix output filenames with this"
    type: string

  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  marked_duplicates_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_mark.bam
