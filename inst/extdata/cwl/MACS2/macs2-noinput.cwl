################################################################
##                            MACS2                           ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: 
label: Last updated 10/2019
hints:
  SoftwareRequirement:
    packages:
    - package: macs2
      version: [ 2.2.4 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: ["macs2", "callpeak"]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: -t
    valueFrom: $(inputs.results_path.path)/$(inputs.fq1.basename)
  - prefix: -n
    valueFrom: $(inputs.fq1.basename)
  - prefix: --outdir
    valueFrom: $(inputs.results_path.path)
  - prefix: -f
    valueFrom: 'BAM'
  - prefix: -g
    valueFrom: '1.2e8'
  - prefix: -B
    valueFrom:
  - prefix: -q
    valueFrom: '0.01'
  - prefix: --nomodel
    valueFrom: 

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory
  fq1:
    label: "This is the only REQUIRED parameter for MACS. The file can be in any supported format specified by --format option."
    type: File

outputs:
  peaks_xls:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.fq1.basename)_peaks.xls
  control_lambda:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.fq1.basename)_control_lambda.bdg
  narrow_peak:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.fq1.basename)_peaks.narrowPeak
  summits_bed:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.fq1.basename)_summits.bed
  macs2_pileup:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.fq1.basename)_treat_pileup.bdg
