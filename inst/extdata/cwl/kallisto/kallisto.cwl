################################################################
##                   Kallisto.cwl                             ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
label: Last updated 09/2019
hints:
  SoftwareRequirement:
    packages:
    - package: kallisto
      version: [ 0.46.0 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: ["kallisto", "quant"]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: -o
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)
  - prefix: -i
    valueFrom: $(inputs.idx_basedir.path)/$(inputs.kallisto_idx)

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  idx_basedir:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory
  idx_basename:
    label: "Basename of the kallisto index files"
    type: string
  kallisto_idx:
    label: "Kallisto index"                                                                                                                                                                               
    type: string 
  fq1:
    label: "Comma-separated list of files containing mate 1s to be aligned"
    type: File
    inputBinding:
      position: 7
  fq2:
    label: "Comma-separated list of files containing mate 2s to be aligned"
    type: File
    inputBinding:
      position: 8
  thread:
    label: "Launch NTHREADS parallel search threads"
    type: int
    inputBinding:
      prefix: --threads 
      position: 1
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory
  bootstrap_samples:
    label: "Number of bootstrap samples"
    type: int
    inputBinding:
      prefix: -b
      position: 2
  fragment_length:
    label: "Estimated average fragment length"
    type: int
    inputBinding:
      prefix: -l
      position: 3
  std_dev:
    label: "Estimated standard deviation of fragment length"
    type: int
    inputBinding:
      prefix: -s
      position: 4

outputs:
  run_info:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)/run_info.json
  kallisto:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)/abundance.h5
  kallisto_file:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)/abundance.tsv
