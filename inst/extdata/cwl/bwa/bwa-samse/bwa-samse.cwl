################################################################
##                 Bwa-samse-Single_Read.cwl                    ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[BWA](http://bio-bwa.sourceforge.net/bwa.shtml): Burrows-Wheeler Alignment Tool"
label: Last updated 04/2021
hints:
  SoftwareRequirement:
    packages:
    - package: bwa
      version: [ 0.7.17 ]
      
################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: ["bwa", "samse"]

arguments:
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
  sai:
    type: File
    inputBinding:
      position: 1
  fq1:
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    type: File
    inputBinding:
      position: 2
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

stdout: $(input.results_path.basename)/$(inputs.SampleName).sam

outputs:
  bwa_samse_sam:
    type: stdout
    
