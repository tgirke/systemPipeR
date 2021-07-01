################################################################
##                Workflow_Hisat2-Single_Read                 ##
################################################################

class: Workflow
cwlVersion: v1.0

################################################################
##              Inputs and Outputs Settings                   ##
################################################################

inputs:
  fq1: File
  idx_basedir: Directory
  idx_basename: string
  SampleName: string
  thread: int
  results_path: Directory

outputs:
  hisat2:
    outputSource: hisat2/hisat2_sam
    type: File
  samtools-view:
    outputSource: samtools-view/samtools_bam
    type: File
  samtools-sort:
    outputSource: samtools-sort/samtools_sort_bam
    type: File
  samtools-index:
    outputSource: samtools-index/samtools_index
    type: File

################################################################
##                Workflow Steps Definitions                  ##
################################################################

steps:
  hisat2:
    in:
      fq1: fq1
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      SampleName: SampleName
      thread: thread
      results_path: results_path
    out: [hisat2_sam]
    run: ./param/cwl/hisat2/hisat2-se/hisat2-mapping-se.cwl
  
  samtools-view:
    in:
      samtools_sam: hisat2/hisat2_sam
      SampleName: SampleName
      results_path: results_path
    out: [samtools_bam]
    run: ./param/cwl/samtools/samtools-view.cwl

  samtools-sort:
    in:
      samtools_bam: samtools-view/samtools_bam
      SampleName: SampleName
      thread: thread
      results_path: results_path
    out: [samtools_sort_bam]
    run: ./param/cwl/samtools/samtools-sort.cwl
   
  samtools-index:
    in:
      samtools_sort_bam: samtools-sort/samtools_sort_bam
      SampleName: SampleName
      results_path: results_path
    out: [samtools_index]
    run: ./param/cwl/samtools/samtools-index.cwl
    
###########
## Notes ##
###########

## If the template it is used in BASH script with the "cwl-runner" or "cwltool", do: 
## Replace the relative PATH for each step for the full PATH, and run:
## "cwltool --outdir <path> workflow_hisat2-se.cwl workflow_hisat2-se.yml"
