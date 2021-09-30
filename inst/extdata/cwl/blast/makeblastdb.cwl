cwlVersion: v1.0
class: CommandLineTool
arguments:
- prefix: -hash_index
- prefix: -parse_seqids
baseCommand: makeblastdb
inputs:
  input_file:
    type: File
    inputBinding:
      prefix: -in
  database_name:
    type: string
    inputBinding:
      prefix: -out
  molecule_type:
    type: string
    inputBinding:
      prefix: -dbtype
outputs:
  out.phd:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).phd
  out.phi:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).phd
  out.phr:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).phr
  out.pin:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).pin
  out.pog:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).pog
  out.psd:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).psd
  out.psi:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).psi
  out.psq:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).psq
