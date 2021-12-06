cwlVersion: v1.0
class: CommandLineTool
baseCommand: cat
arguments:
  - prefix: 
    valueFrom: $(inputs.out) | head -n 10 | cut -f 2
inputs:
  out:
    type: File
results_path:
  out:
    type: Directory  
stdout: $(input.results_path.basename)/top_hits

outputs:
  top_hits:
    type: stdout

