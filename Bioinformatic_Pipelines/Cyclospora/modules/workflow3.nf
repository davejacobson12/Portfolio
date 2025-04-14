#!/usr/bin/env nextflow
nextflow.enable.dsl=2



// process runModule3 {
//   label 'clusteringEnv'
//   label 'process_medium'
//   publishDir "${clusterFolder}", pattern: "*RESULTING*.txt", mode: "copy"

//   input:
//   path latestMatrix

//   output:
//   path "*RESULTING*txt"

//   script:
//   """

//   bash $baseDir/MODULE_3_clustering_V3_nextflow.sh ${latestMatrix} ${params.stringency}

//   """ 
// }



process runModule3_v2 {
  label 'clusteringEnv'
  label 'process_medium'
  publishDir "${params.clusterFolder}", pattern: "*RESULTING*.txt", mode: "copy"

  input:
  path latestMatrix

  output:
  path "*RESULTING*txt"

  script:
  """

  bash $projectDir/bin/MODULE_3_clustering_V4.sh ${latestMatrix}

  """ 
}


