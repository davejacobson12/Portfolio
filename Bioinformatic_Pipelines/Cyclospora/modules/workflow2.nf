#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process runModule2{
  label 'pythonEnv_wf2'
  label 'process_high_memory'
  publishDir "${params.matrixFolder}", pattern: "*H_matrix.csv", mode: "copy"

  input:
  path latestSheet

  output:
  path '*H_matrix.csv', emit: matrixOut

  script:
  """
  python3 $projectDir/EUKARYOTYPING/CycloTyping/BayesianHeuristic.py -infile ${latestSheet} -outbh ${latestSheet.simpleName}_bh_matrix.csv -outhh ${latestSheet.simpleName}_hh_matrix.csv -outB ${latestSheet.simpleName}_B_matrix.csv -outH ${latestSheet.simpleName}_H_matrix.csv -sampmeetlocirequire samplesMeetCutoff.txt -expectlocifile ${params.wf2MarkerList} -expectlocinumber ${params.wf2MinLociNum}

  """
}

