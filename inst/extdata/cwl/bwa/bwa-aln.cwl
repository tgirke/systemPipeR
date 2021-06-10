################################################################
##                 Bwa-aln-Single_Read.cwl                    ##
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

baseCommand: ["bwa", "aln"]

arguments:
  - prefix: -n
    valueFrom: $(inputs.ndis)
  - prefix: -o
    valueFrom: $(inputs.o)
  - prefix: -e
    valueFrom: $(inputs.e)
  - prefix: -l
    valueFrom: $(inputs.l)
  - prefix: -k
    valueFrom: $(inputs.k)
  - prefix: -t
    valueFrom: $(inputs.thread)
  - prefix: 
    valueFrom: $(inputs.idx_basedir.path)/$(inputs.idx_basename)

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  ndis:
    label: "Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT."
    type: int
  o:
    label: "Maximum number of gap opens [1]"
    type: int
  e:
    label: "Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps) [-1]"
    type: int
  l:
    label: "Take the first INT subsequence as seed."
    type: int
  k:
    label: "Launch NTHREADS parallel search threads"
    type: int
  thread:
    label: "Launch NTHREADS parallel search threads"
    type: int
  idx_basedir:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory
  idx_basename:
    label: "Basename of the bowtie2 index files"
    type: string
  fq1:
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    type: File
    inputBinding:
      position: 1
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

stdout: $(input.results_path.basename)/$(inputs.SampleName).sai

outputs:
  bwa_aln_sai:
    type: stdout
    
