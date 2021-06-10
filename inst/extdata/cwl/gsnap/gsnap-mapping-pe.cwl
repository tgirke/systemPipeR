################################################################
##                     gsnap-Paired_end.yml                   ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "gmapR (1.26.0): An R interface to the GMAP/GSNAP/GSTRUCT suite - DOI: 10.18129/B9.bioc.gmapR"
label: Last updated 09/2019

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [Rscript]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: --vanilla
    valueFrom: $(inputs.rscript_gsnap_mapping)
  - prefix: --output_file
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).gsnap

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  rscript_gsnap_mapping:
    type: File
  fq1:
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    type: File
    inputBinding:
      prefix: --readfile1
  fq2:
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    type: File
    inputBinding:
      prefix: --readfile2
  thread:
    label: "Launch NTHREADS parallel search threads"
    type: int
    inputBinding:
      prefix: --nthreads
  molecule:
    label: ""
    type: string
    inputBinding:
      prefix: --molecule
  max_mismatches:
    label: ""
    type: int
    inputBinding:
      prefix: --max_mismatches
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  gsnap_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).gsnap.sam.bam
  gsnap_bam_bai:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).gsnap.sam.bam.bai
