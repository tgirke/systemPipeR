################################################################
##                     Bwa-mem-Paired_end.yml                 ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[BWA-MEM](http://bio-bwa.sourceforge.net/bwa.shtml): Burrows-Wheeler Alignment Tool"
label: Last updated 07/1019
hints:
  SoftwareRequirement:
    packages:
    - package: bwa
      version: [ 0.7.17 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: ["bwa", "mem"]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: -M
  - prefix: -R
    valueFrom: '''@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1'''

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  reference_file:
    type: File
    secondaryFiles: [ .amb, .ann, .bwt, .pac, .sa ]
    inputBinding:
      position: 2
  fq1:
    type: File
    inputBinding:
      position: 3
  fq2:
    type: File
    inputBinding:
      position: 4
  thread:
    type: int
    inputBinding:
      prefix: -t
      position: 1
  SampleName:
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

stdout: $(inputs.results_path.basename)/$(inputs.SampleName).sam

outputs:
  bwa_men_sam:
    type: stdout
    
###########
## Notes ##
###########

## If the template its used in bash script with the "cwl-runner", run: 
## Remove two single quotation: '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1'
## "cwl-runner --outdir <path> bwaPE.cwl bwaPE.yml"
