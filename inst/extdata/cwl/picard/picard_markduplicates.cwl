################################################################
##                           picard.cwl                       ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[picard](http://broadinstitute.github.io/picard/): Identifies duplicate reads"
label: Last updated 09/2019
hints:
  SoftwareRequirement:
    packages:
      - package: picard/2.18.3

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [picard]


requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]
    
arguments:
  - prefix:
    valueFrom: MarkDuplicates
  - prefix: I=
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).sorted.bam
  - prefix: O=
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).marked_duplicates.bam
  - prefix: M=
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).marked_dup_metrics.txt

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  SampleName:
    label: "Prefix output filenames with this"
    type: string
  samtools_sort_bam:
    label: "Filename of sorted sam/bam"
    type: File
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  marked_duplicates_bam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).marked_duplicates.bam
  marked_dup_metrics_txt:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).marked_dup_metrics.txt
