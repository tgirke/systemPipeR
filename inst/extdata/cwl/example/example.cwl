cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
inputs:
  message:
    type: string
    inputBinding: 
      position: 1
  SampleName:
    type: string
  results_path:
    type: Directory
outputs:
    string:
        type: stdout
stdout: $(inputs.results_path.basename)/$(inputs.SampleName).txt

