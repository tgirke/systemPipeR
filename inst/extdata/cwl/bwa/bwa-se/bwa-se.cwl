################################################################
##                 Bwa-mem-Single_Read.cwl                    ##
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

arguments:
  - prefix: 
    valueFrom: -M
    position: 2

inputs:
  reference_file:
    type: File
    secondaryFiles: [ .amb, .ann, .bwt, .pac, .sa ]
    inputBinding:
      position: 3
  fq1:
    type: File
    inputBinding:
      position: 4
  nthreads:
    type: int
    inputBinding:
      prefix: -t
      position: 1
  SampleName:
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

stdout: $(input.results_path.basename)/$(inputs.SampleName).sam

outputs:
  bwa_mem_sam:
    type: stdout
    
