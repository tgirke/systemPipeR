################################################################
##                     Bwa-mem-Paired_end.cwl                 ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[BWA-MEM](http://bio-bwa.sourceforge.net/bwa.shtml): Burrows-Wheeler Alignment Tool"
label: Last updated 09/2019
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
  - prefix: 
    valueFrom: $(inputs.idx_basedir.path)/$(inputs.idx_basename)

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  idx_basedir:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory
  idx_basename:
    label: "Basename of the bowtie2 index files"
    type: string
  fq1: 
    label: "Comma-separated list of files containing mate 1s to be aligned"
    type: File
    inputBinding:
      position: 1
  fq2:
    label: "Comma-separated list of files containing mate 2s to be aligned"
    type: File
    inputBinding:
      position: 2
  thread:
    label: "Launch NTHREADS parallel search threads"
    type: int
    inputBinding:
      prefix: -t
  SampleName:
    label: "Filename to write output to"
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
## Adjust
## outputs:
##  bwa_men_sam:
##    type: stdout
## "cwltool --outdir <path> bwaPE.cwl bwaPE.yml"
