################################################################
##                     Bwa-mem-Paired_end.cwl                 ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[BWA-MEM](http://bio-bwa.sourceforge.net/bwa.shtml): Burrows-Wheeler Alignment Tool"
label: Last updated 10/2019
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
    position: 1

  - prefix: 
    valueFrom: $(inputs.data_path.path)/$(inputs.ref_name)
    position: 2

  - prefix: -t
    valueFrom: $(inputs.thread)
    position: 3

  - prefix:
    valueFrom: $(inputs.fq1)
    position: 4

  - prefix: 
    valueFrom: $(inputs.fq2)
    position: 5

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  data_path:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory

  ref_name:
    label: "Basename of the bowtie2 index files"
    type: string

  fq1: 
    label: "Comma-separated list of files containing mate 1s to be aligned"
    type: File

  fq2:
    label: "Comma-separated list of files containing mate 2s to be aligned"
    type: File

  thread:
    label: "Launch NTHREADS parallel search threads"
    type: int

  SampleName:
    label: "Filename to write output to"
    type: string

  results_path:
    label: "Path to the results directory"
    type: Directory

stdout:  $(inputs.results_path.basename)/$(inputs.SampleName).sam

outputs:
  bwa_men_sam:
    type: stdout
    
###########
## Notes ##
###########
## If the template its used in bash script with the "cwl-runner", run: 
## "cwltool --outdir <path> bwaPE.cwl bwaPE.yml"
