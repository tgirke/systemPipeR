cwlVersion: v1.0
class: CommandLineTool
baseCommand: blastdbcmd
arguments:
  - prefix: 
    valueFrom: -db  $(inputs.database_name) -entry_batch $(inputs.entry_batch) | head -n 10
inputs:
  database_name:
    type: File
  entry_batch:
    type: File

stdout: $(input.results_path.basename)/sequences_top_hits.fasta

outputs:
  blastdbcmd:
    type: stdout
      

