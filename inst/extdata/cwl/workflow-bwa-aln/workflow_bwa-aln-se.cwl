################################################################
##                Workflow_bwa-aln-Single_Read                 ##
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
  ndis: int
  o: int
  e: int
  l: int
  k: int

outputs:
  bwa_aln_sai:
    outputSource: bwa_aln_sai/bwa_aln_sai
    type: File
  bwa_samse_sam:
    outputSource: bwa_samse_sam/bwa_samse_sam
    type: File


################################################################
##                Workflow Steps Definitions                  ##
################################################################

steps:
  bwa_aln_sai:
    in:
      fq1: fq1
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      SampleName: SampleName
      thread: thread
      results_path: results_path
      ndis: ndis
      o: o
      e: e
      l: l
      k: k
    out: [bwa_sai]
    run: bwa/bwa-aln.cwl
  
  bwa_samse_sam:
    in:
      sai: bwa_aln_sai/bwa_aln_sai
      fq1: fq1
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      SampleName: SampleName
      thread: thread
      results_path: results_path
    out: [bwa_sam]
    run: bwa/bwa-samse.cwl

    
###########
## Notes ##
###########

## If the template it is used in BASH script with the "cwl-runner" or "cwltool", do: 
## Replace the relative PATH for each step for the full PATH, and run:
## "cwltool --outdir <path> workflow_bwa-aln.cwl workflow_bwa-aln.yml"
