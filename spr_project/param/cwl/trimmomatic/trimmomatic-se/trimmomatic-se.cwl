################################################################
##                   Trimmomatic-se.cwl                       ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[TRIMMOMATIC](http://www.usadellab.org/cms/?page=trimmomatic)"
label: Last updated 09/2019
hints:
  SoftwareRequirement:
    packages:
    - package: trimmomatic
      version: [ 0.0.36 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: ["trimmomatic", "SE"]
## baseCommand: ["java -jar trimmomatic-0.36.jar", "SE"]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - valueFrom: -phred33
  - valueFrom: $(inputs.fq1)
  - valueFrom: $(inputs.results_path.basename)/$(inputs.SampleName).trimmed.fastq.gz
    
################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  fq1:
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    type: File
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory
  thread:
    label: "Launch NTHREADS parallel search threads"
    type: int
    inputBinding:
      prefix: -threads
  minlen:
    label: "Reads of lengths shorter will be discarded"
    type: string
    inputBinding:
      position: 1
  leading:
    label: "Minimum quality required to keep a base in the beginning of the read"
    type: string
    inputBinding:
      position: 2
  trailing:
    label: "Minimum quality required to keep a base in the end of the read"
    type: string
    inputBinding:
      position: 3
  slidingwindow:
    label: "Window size and average quality required"
    type: string
    inputBinding:
      position: 4

outputs:
  trimmomatic:
    type: File
    outputBinding:
      glob: $(inputs.results_path.basename)/$(inputs.SampleName).trimmed.fastq.gz
