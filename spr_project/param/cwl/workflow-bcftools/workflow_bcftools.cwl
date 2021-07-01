################################################################
##                       Workflow_BCFtools                    ##
################################################################

class: Workflow
cwlVersion: v1.0

requirements:
  MultipleInputFeatureRequirement: {}

################################################################
##              Inputs and Outputs Settings                   ##
################################################################

inputs:
  idx_basename: string
  SampleName: string
  data_path: Directory
  results_path: Directory
  thread: int
  bam: File

outputs:
  samtools_markduplicates:
    outputSource: mark/marked_duplicates_bam
    type: File
  samtools-sort:
    outputSource: sort/samtools_sort_bam
    type: File
  samtools-index:
    outputSource: index/samtools_index
    type: File
  bcftools:                        
    outputSource: raw_call/bcftools
    type: File
  bcftools_calls:
    outputSource: vcf_call/bcftools_call
    type: File

################################################################
##                Workflow Steps Definitions                  ##
################################################################

steps:
  mark:
    in:
      bam: bam
      SampleName: SampleName
      results_path: results_path
    out: [marked_duplicates_bam]
    run: ./param/cwl/workflow-bcftools/samtools_markduplicates.cwl
  
  sort:
    in:
     # source: mark/marked_duplicates_bam
      SampleName: SampleName
      results_path: results_path
    out: [samtools_sort_bam ]
    run: ./param/cwl/workflow-bcftools/samtools-sort.cwl
   
  index:
    in:
     # source: sort/samtools_sort_bam
      SampleName: SampleName
      results_path: results_path
    out: [samtools_index]
    run: ./param/cwl/workflow-bcftools/samtools-index.cwl

  raw_call:
    in:
     # source: [sort/samtools_sort_bam, index/samtools_index]
      SampleName: SampleName
      results_path: results_path
      data_path: data_path
      idx_basename: idx_basename
    out: [bcftools]
    run: ./param/cwl/workflow-bcftools/bcftools.cwl  

  vcf_call:
    in:
     # source: raw_call/bcftools
      SampleName: SampleName
      results_path: results_path
    out: [bcftools_call]
    run: ./param/cwl/workflow-bcftools/bcftools-call.cwl
    
###########
## Notes ##
###########

## When run with 'cwltool', add following: --relax-path-checks --leave-outputs
# Also, uncomment lines with "source:" in steps 
