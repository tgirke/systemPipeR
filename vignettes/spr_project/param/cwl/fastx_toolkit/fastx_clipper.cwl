################################################################
##                       fastx_clipper                       ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "Removing sequencing adapters / linkers"
label: Last updated 04/2019
hints:
  SoftwareRequirement:
    packages:
    - package: fastx_toolkit
      version: [ 0.0.13 ]

################################################################
##            baseCommand and arguments definitions           ##
################################################################

baseCommand: [ fastx_clipper ]

requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.results_path)

arguments:
  - prefix: -Q
    valueFrom: '33'
    position: 1
  - prefix: -l
    valueFrom: '25'
    position: 2
  - prefix: -o
    valueFrom: $(inputs.results_path.basename)/$(inputs.SampleName).clipped.fq 
    position: 5

################################################################
##             Inputs and Outputs Settings                    ##
################################################################

inputs:
  gunzip:
    type: File
    inputBinding:
      prefix: "-i"
      position: 4
  adaptor_fastx_toolkit:
    label: "Sequence of adaptor to clip from the 3'"
    type: string
    inputBinding:
      prefix: "-a"
      position: 3
  SampleName:
    label: "Sample name"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  fastx_clipper_file:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName).clipped.fq

