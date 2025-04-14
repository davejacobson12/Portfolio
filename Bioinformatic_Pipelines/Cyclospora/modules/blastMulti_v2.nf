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

process BLAST_MULTI_V2 {
    tag "$meta.id"
   label 'mappingSoftware'
   label 'process_high'

  // publishDir "${newMethGenotypes}", pattern: "*haplotypes", mode: "copy"

  input:

  tuple val(meta), path(depthFile), path(multiFile), path(assignHaps_cov)
//   path(nhr)
//   path(nin)
//   path(nsq)
//   path(fasta)
  
//   each(nhr)
//   each(nin)
//   each(nsq)
  each(fasta)

  output:
  // path("*haplotypes")
  tuple val(meta), path("*finalBlast.fasta"), emit: blastFasta, optional: true
  // tuple val("${indMulti.baseName}"), path("${indMulti.baseName}.Num.fasta"), emit:testTuple

  shell:
  '''
  cp !{fasta} toBlast.fasta
  
  makeblastdb -in toBlast.fasta -input_type fasta -dbtype nucl

  echo !{multiFile} >> myFiles
  echo !{meta.id} > checkSampleName.txt

  mkdir nonJunc_mapping
  cp !{assignHaps_cov} nonJunc_mapping

  cat myFiles | tr " " "\n" | awk -F "." '{print$1}'  > test

  cat !{multiFile} >> !{meta.id}_hapsBlast.fasta

  blastn -db toBlast.fasta -query !{meta.id}_hapsBlast.fasta -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -num_threads !{params.max_cpus} -max_target_seqs 1 -word_size 7 -dust no -soft_masking false  -outfmt "6 sseqid qseqid" -out !{meta.id}.blastout

  cat !{meta.id}.blastout | sed 's/_Read/-Read/g' | awk -F "-" '{print$1"\t"$2}' |awk -F "\t" '{print$1"\t"$2"\t"$3}' | sort -k 2 > !{meta.id}.blastFormat

  cat !{depthFile} | awk -F "," '{print$1"\t"$2}' | sort -k 1 >> !{meta.id}.myDepth

  awk -F "\t" '{print$2}' !{meta.id}.blastFormat | while read line; do grep $line -w -h !{meta.id}.myDepth !{meta.id}.blastFormat | tr "\n" "\t" ; echo ; done | awk -F "\t" '{print$3"\t"$1"\t"$2"\t"$4"\t"$5}' > !{meta.id}.hapDepth  

  mySamp=!{meta.id}

  awk -v var=!{meta.id} -F "_" '{print$0, "\t"var"."$1"_"$2"_"$3"_"$4"_"$5}'  !{meta.id}.hapDepth > tester

  awk -F "\t" '{print$6"\t"$3"\t"$1"\t"$4"\t"$4"-"$5}' tester | while read search actDepth hap clusterName readName
  do
  FILE=nonJunc_mapping/$search.assignHaps_coverageCutoff
  echo "$FILE"
  if test -f "$FILE";
  then
  dynamicCutoff=`cat nonJunc_mapping/$search.assignHaps_coverageCutoff`
  echo $dynamicCutoff >> sampleMarkerDepth
  echo !{meta.id}.$clusterName.multiOut >> multiNameTesting
  else
  dynamicCutoff=!{params.assignHap_min}
  echo "no file exists"
  fi
  if [ $actDepth -gt $dynamicCutoff ] && [ $actDepth -gt !{params.assignHap_min} ]
  then 
  cat !{meta.id}.$clusterName.multiOut >> !{meta.id}.finalBlast.fasta
  echo $hap >> !{meta.id}.haplotypes
  fi
  done



  '''
}

