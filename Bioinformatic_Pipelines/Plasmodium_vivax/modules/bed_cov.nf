
process BED_GENOMECOV {
  tag "$meta.id"
  label 'process_high'
  label 'bedCov'

  publishDir "$projectDir/coverage_estimates/bed_coverage", pattern: "*default*txt", mode: "copy"
  // publishDir "$projectDir/coverage_estimates", pattern: "*bedGraph*txt", mode: "copy"

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path("*default*txt"), emit: defaultBED
  // tuple val(meta), path("*bedGraph*txt"), emit: bedGraph
  // path "*bam", emit: mappedBowtie2
  // path "*mapStats*", emit: mapStat
  // path "*bedtoolsCoverage.txt"
 


  script:
  """
  bedtools genomecov -ibam ${bam} -g $projectDir/REFERENCES/pvp01/pvp01_fullGenome_genomeFile.txt > ${meta.id}.default_bedtoolsCoverage.txt

  """
}


process GENOMECOV_SUMMARY {
  tag "$meta.id"
  label 'process_low'
  label 'reporting'

  publishDir "$projectDir/coverage_estimates/perChrom_coverage", pattern: "*perChrom_coverage*", mode: "copy"

  input:
  tuple val(meta), path(defaultBED)

  output:
  tuple val(meta), path("*_pvp01_perChrom_coverage.csv"), emit: chromCov

  script:
  """
  Rscript $projectDir/bin/bedCov_summary.R ${meta.id} ${defaultBED}
  """
}


process MAPPING_SUMMARY {
  tag "$meta.id"
  label 'process_low'
  

  publishDir "$projectDir/coverage_estimates/mapping_stats", pattern: "*read_summary.txt", mode: "copy"

  input:
  tuple val(meta), path(bbdukLog), path(samtoolsStat), path(bowtie2Log)

  output:
  path("*read_summary.txt"), emit: summary

  shell:
  '''
  echo "Stat:" > file.txt
  echo "Input Reads:" >> file.txt
  echo "After BBDUK Reads:" >> file.txt
  echo "Human Mapped:" >> file.txt
  echo "Human Unmapped:" >> file.txt
  echo "pvp01 Mapped:" >> file.txt
  
  echo !{meta.id} > stats.txt
  grep "Input:" !{bbdukLog} | awk '{print$2}' >> stats.txt
  grep "Result:" !{bbdukLog} | awk '{print$2}' >> stats.txt
  grep "reads mapped:" !{samtoolsStat} | awk '{print$4}' >> stats.txt
  grep "reads unmapped:" !{samtoolsStat} | awk '{print$4}' >> stats.txt
  grep "overall alignment" !{bowtie2Log} | awk '{print$1}' >> stats.txt

  paste file.txt stats.txt > !{meta.id}.read_summary.txt

  '''
}


process COMBINE_SUMMARY {
  label 'process_low'
  

  publishDir "$projectDir/coverage_estimates/mapping_stats", pattern: "*pvp01_swgs_mapStats.csv", mode: "copy"

  input:
  path(readSummary)

  output:
  path("*pvp01_swgs_mapStats.csv"), emit: chromCov

  script:
  """

  Rscript $projectDir/bin/format_map_depth.R

  """
}
