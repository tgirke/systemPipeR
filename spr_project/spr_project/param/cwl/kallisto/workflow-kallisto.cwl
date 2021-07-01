######################################
## Workflow_Kallisto-Paired-end-Reads ##
######################################

class: Workflow
cwlVersion: v1.0

################################################################
##              Inputs and Outputs Settings                   ##
################################################################

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

inputs:
  fq1: File
  fq2: File
  idx_basedir: Directory
  idx_basename: string
  kallisto_idx: string
  SampleName: string
  thread: int
  results_path: Directory
  bootstrap_samples: int                                                                                                                                                                                      
  fragment_length: int                                                                                                                                                                                         
  std_dev: int

outputs:
  index:
    outputSource: index/kallisto-index
    type: File
  kallisto:
    outputSource: kallisto/kallisto
    type: File

################################################################
##                Workflow Steps Definitions                  ##
################################################################

steps:
  index:
    in:
      idx_basedir: idx_basedir
      idx_basename: idx_basename
      kallisto_idx: kallisto_idx
    out: [kallisto-index]
    run: ./param/cwl/kallisto/kallisto-index.cwl
    
  kallisto:
    in:
      fq1: fq1                                                                                                                                                                                                                  
      fq2: fq2                                                                                                                                                                                                                  
      idx_basedir: idx_basedir                                                                                                                                                                                                  
      idx_basename: idx_basename
      kallisto_idx: kallisto_idx
      SampleName: SampleName                                                                                                                                                                                                    
      thread: thread                                                                                                                                                                                                            
      results_path: results_path
      bootstrap_samples: bootstrap_samples
      fragment_length: fragment_length
      std_dev: std_dev
      out: [kallisto_file, kallisto, run_info]
    run: ./param/cwl/kallisto/kallisto.cwl

