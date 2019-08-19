################################################################
##                     Bowtie2-Mapping                        ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): Fast and sensitive read alignment"
label: Last updated 07/2019
hints:
  SoftwareRequirement:
    packages:
    - package: bowtie2
      version: [ 2.3.4.1 ]
      
################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [ bowtie2 ]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: -S
    valueFrom: $(inputs.results_path.basename)/$(inputs.SampleName).sam
    position: 1
  - prefix: -x
    valueFrom: $(inputs.bowtie2_idx_basedir.path)/$(inputs.bowtie2_idx_basename)
    position: 2
  - prefix: -k
    valueFrom: '50'
    position: 3
  - prefix: 
    valueFrom: --non-deterministic
    position: 4

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  bowtie2_idx_basedir:
    label: "Path to the directory the index for the reference genome"
    type: Directory
  bowtie2_idx_basename:
    label: "Basename of the bowtie2 index files"
    type: string
  fq1:
    type: File
    inputBinding:
      prefix: "-1"
      position: 5
  fq2:
    type: File
    inputBinding:
      prefix: "-2"
      position: 6
  thread:
    type: int
    inputBinding:
      prefix: -p
      position: 7
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  bowtie2_sam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName).sam
