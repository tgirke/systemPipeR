class: Workflow
cwlVersion: v1.0
inputs:
  input_file: File
  database_name: string
  molecule_type: string
outputs:
  out.phd:
    outputSource: makeblastdb/out.phd
    type: File
  out.phi:
    outputSource: makeblastdb/out.phi
    type: File
  out.phr:
    outputSource: makeblastdb/out.phr
    type: File
  out.pin:
    outputSource: makeblastdb/out.pin
    type: File
  out.pog:
    outputSource: makeblastdb/out.pog
    type: File
  out.psd:
    outputSource: makeblastdb/out.psd
    type: File
  out.psi:
    outputSource: makeblastdb/out.psi
    type: File
  out.psq:
    outputSource: makeblastdb/out.psq
    type: File
steps:
  makeblastdb:
    in:
      input_file: input_file
      database_name: database_name
      molecule_type: molecule_type
    out: '[out.phd, out.phi, out.phr, out.pin, out.pog, out.psd, out.psi, out.psq]'
    run: blast/makeblastdb.cwl
