// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process HAPS_BLAST {
    tag "$meta.id"
    label 'mappingSoftware'
    label 'process_high'
  
  publishDir "${params.tmpGenotypes}", pattern: "*GENOTYPE_noJunction", mode: "copy"

  input:
  tuple val(meta), path(myHaps)
  each(fasta)

  output:
  path "*GENOTYPE_noJunction", emit: genotypeNoJunction
  // tuple val(meta), path ("${meta.id}.GENOTYPE_noJunction"), emit: genotypeNoJunction_tuple
//   path("*GENOTYPE_noJunction")
  path("endMarker"), emit: SNP_hapsDone

  script:
  """

  hapLen="\$(wc -l ${myHaps} | awk '{print\$1}')"

  if [ \$hapLen -gt 0 ]
  then

  
  makeblastdb -in ${myHaps} -input_type fasta -dbtype nucl

  blastn -db ${myHaps} -query ${fasta} -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -num_threads ${params.max_cpus} -out ${meta.id}.GENOTYPE_noJunction -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid pident mismatch gapopen gaps sseq evalue bitscore"
  else
  touch ${meta.id}.GENOTYPE_noJunction
  fi

  echo "done" >> endMarker
  """
}

  // cp ${meta.id}.GENOTYPE_noJunction $projectDir/HAPLOTYPE_CALLER_CYCLO_V2/SPECIMEN_GENOTYPES/${meta.id}

