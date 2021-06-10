################################################################
##                           blastp                           ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): Fast and sensitive read alignment"
label: Last updated 09/2019
hints:
  SoftwareRequirement:
    packages:
    - package: bowtie2
      version: [ ]
      
################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: blastp

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]


################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:

outputs:

