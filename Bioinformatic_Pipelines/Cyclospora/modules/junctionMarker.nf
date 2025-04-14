#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Give path to folders used in this module
processReads = file("$projectDir/HAPLOTYPE_CALLER_CYCLO_V2/TMP2/PROCESSED_READS")
refHapsDir = file("$projectDir/HAPLOTYPE_CALLER_CYCLO_V2/TMP2/REF_HAPS")
tmpGenotypes = file("$projectDir/HAPLOTYPE_CALLER_CYCLO_V2/TMP2/TMP_GENOTYPES")
NEW_HAPS_DIR = file("$projectDir/HAPLOTYPE_CALLER_CYCLO_V2/REF_SEQS/BLASTING/NEW_HAPS")
Junc_mapping = file("$projectDir/HAPLOTYPE_CALLER_CYCLO_V2/TMP2/JUNC_MAPPING")
NEW_HAPS_DIR_Junc = file("$projectDir/HAPLOTYPE_CALLER_CYCLO_V2/TMP2/JUNC_New_haps")

process makeJunctionRefs {
  label 'process_low'
  input:
  each(newHaps)
  
  output:
  path 'newJunctions.fasta', emit: completeFile

  script:
  """
  cat ${newHaps} > newJunctions.fasta
  """
}

process finalJunctionNewRefs {
  label 'process_low'

  publishDir "${Junc_mapping}", pattern: "fullRefs_junction.fasta", mode: "copy"

  input:
  path allNew

  output:
  path 'fullRefs_junction.fasta', emit: outJunctionRefs
  path("fullRefs_junction.fasta")

  script:
  """
  cat ${params.origNewHaps}/*Junction*.fasta ${params.fullJunction_refs} > fullRefs_junction.fasta
  
  """
}

process mappingJunction {
  tag "$meta.id"
  label 'mappingSoftware'
  label 'process_medium'

  input:
  // path(refHaps)
  // each(sample_files)
  each(refHaps)
  tuple val(meta), path(sample_files)
  output:
  path '*mapped_only_single_line2.fasta', emit: mappedSingleLine
  tuple val(meta), path('*mapped_only_single_line2.fasta'), emit: mappedSingleLine_tuple

  script:
  """
  export LC_ALL=C
  cp ${refHaps} refHaps.fasta
  bwa index refHaps.fasta
  bwa mem -t ${params.threads} refHaps.fasta ${sample_files} > ${meta.id}.alignment.sam
  samtools view -h -F 4 ${meta.id}.alignment.sam | samtools view -bS > ${meta.id}.mapped_only.sam 

  samtools view ${meta.id}.mapped_only.sam | awk '{print("@"\$1"\\n"\$10"\\n+\\n"\$11)}' > ${meta.id}.mapped_only.fastq

  miraconvert ${meta.id}.mapped_only.fastq ${meta.id}.mapped_only.fasta

  cat ${meta.id}.mapped_only.fasta | awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' | awk 'NF > 0' > ${meta.id}.mapped_only_single_line2.fasta
  """
}



process getAssembleSeqs {
  tag "$meta.id"
  label 'mappingSoftware'
  label 'process_medium'

  publishDir "${Junc_mapping}", pattern: "*READS_MAPPING.txt", mode: "copy"
  publishDir "${Junc_mapping}", pattern: "*cluster_these_seqs.fasta", mode: "copy"

  input:
  // path(mappedSingle)
  tuple val(meta), path(mappedSingle)

  output:
  path "*cluster_these_seqs.fasta", emit: forClustering
  path("*READS_MAPPING.txt")
  path("*cluster_these_seqs.fasta")
  tuple val(meta), path("*cluster_these_seqs.fasta"), emit: forClustering_tuple

  shell:
  '''
  while read line
        do
      if [[ ${line:0:1} == '>' ]]
          then
              outfile=${line#>}.fa
              echo $line > $outfile
          else
             echo $line >> $outfile
         fi
      done < !{mappedSingle}


  find -type f -name '*.fa' | wc -l > !{meta.id}.READS_MAPPING.txt

  find -type f -name  '*.fa' | while read my_fastas; do grep -l 'CCATCTACAGC.*AACAC' $my_fastas || true ;done >> !{meta.id}.list_of_sequences_to_review.txt
  find -type f -name  '*.fa' | while read my_fastas; do grep -l 'GTGTT.*GCTGTAGATGG' $my_fastas || true ;done >> !{meta.id}.list_of_sequences_to_review.txt
  find -type f -name  '*.fa' | while read my_fastas; do grep -l 'TCTACAGC.*AACACGATC' $my_fastas || true ;done >> !{meta.id}.list_of_sequences_to_review.txt
  find -type f -name  '*.fa' | while read my_fastas; do grep -l 'GATCGTGTT.*GCTGTAGA' $my_fastas || true ;done >> !{meta.id}.list_of_sequences_to_review.txt

  sort -u !{meta.id}.list_of_sequences_to_review.txt > !{meta.id}.2.list_of_sequences_to_review.txt
  touch !{meta.id}.assemble_these_seqs 
  cat !{meta.id}.2.list_of_sequences_to_review.txt | while read myReads; do cat $myReads >> !{meta.id}.assemble_these_seqs ; done

  awk '/^>/{print ">Contig_" ++i; next}{print}' < !{meta.id}.assemble_these_seqs > !{meta.id}.cluster_these_seqs.fasta
  '''
}


process clusterSeqs_1 {
  label 'mappingSoftware'
  label 'process_medium'

  //errorStrategy 'ignore'
  publishDir "${Junc_mapping}", pattern: "*map_to_this.fastq", mode: "copy"
  
  input:
 // path(toCluster)
  tuple val(meta), path(toCluster)

  output:
  // path "*map_to_this.fastq", emit: forManifest, optional: true
  path("*map_to_this.fastq")
  tuple val(meta), path("${meta.id}.map_to_this.fastq"), emit: forManifest_tuple, optional: true

  script:
  """
  toCluster=`wc -l ${toCluster} | awk '{print\$1}'`

  if [ \$toCluster -gt 0 ]
    then
      cd-hit-est -i ${toCluster} -o ${meta.id}.junction_clusters.fasta -c 1 -g 1 -d 0 -T ${params.threads} -s 1 -M ${params.RAM}
      cap3 ${meta.id}.junction_clusters.fasta -o 95 -p 99.99
    else
      touch ${meta.id}.junction_clusters.fasta.cap.contigs
      touch ${meta.id}.junction_clusters.fasta.cap.singlets
  fi 

  contigs=`wc -l ${meta.id}.junction_clusters.fasta.cap.contigs | awk '{print\$1}'`
  singlets=`wc -l ${meta.id}.junction_clusters.fasta.cap.singlets | awk '{print\$1}'`

  if [ \$contigs -gt 0 ]
   then
    cat ${meta.id}.junction_clusters.fasta.cap.contigs > ${meta.id}.map_to_this.fasta 
  else
   touch ${meta.id}.map_to_this.fasta 
  fi

  if [ \$singlets -gt 0 ] 
   then 
    cat ${meta.id}.junction_clusters.fasta.cap.singlets >> ${meta.id}.map_to_this.fasta 
  fi

  cat ${meta.id}.map_to_this.fasta | awk '/^>/{print ">new_junction_" ++i "_sequence"; next}{print}' > ${meta.id}.map_to_this_2.fasta 
  seqtk seq -F '#' ${meta.id}.map_to_this_2.fasta > ${meta.id}.map_to_this.fastq
  """
}

process miraTrials {
  tag "$meta.id"
  label 'mappingSoftware'
  label 'process_medium'

//  beforeScript 'export LC_ALL=C'

  errorStrategy 'ignore'
  publishDir "${Junc_mapping}", pattern: "*DRAFT_0_junction_sequences.fasta", mode: "copy"
  publishDir "${Junc_mapping}", pattern: "*REV_three_prime_trim.fastq", mode: "copy"
  publishDir "${Junc_mapping}", pattern: "*newly_found_junction.fastq", mode: "copy"
  publishDir "${processReads}", pattern: "${meta.id}.clean_merged.fastq", mode: "copy"
  
  input:
  // path(myClean)
  tuple val(meta), path(cleanReads), path(myClean)

  output:
  path "*DRAFT_0_junction_sequences.fasta", emit: draft0_out
  path("*DRAFT_0_junction_sequences.fasta")
  path("*REV_three_prime_trim.fastq")
  path("*newly_found_junction.fastq")
  path("${meta.id}.clean_merged.fastq")

  script:
  """
  export LC_ALL=C
  cp ${cleanReads} sample.clean_merged.fastq.gz
  gunzip sample.clean_merged.fastq.gz
  mv sample.clean_merged.fastq ${meta.id}.clean_merged.fastq
  sed "s/data = clean_merged.fastq/data = ${meta.id}.clean_merged.fastq/g" ${params.miraManifest} > ${meta.id}.miraFirst.txt
  sed  "s/data = map_to_this.fastq/data = ${meta.id}.map_to_this.fastq/g" ${meta.id}.miraFirst.txt > ${meta.id}.miraSecond.txt 
  sed "s/project = novel_junction_finder/project = ${meta.id}.novel_junction_finder/g" ${meta.id}.miraSecond.txt > ${meta.id}.miraFinal.txt 

  mira ${meta.id}.miraFinal.txt
  miraconvert ${meta.id}.novel_junction_finder_assembly/${meta.id}.novel_junction_finder_d_results/${meta.id}.novel_junction_finder_out_new_junc.unpadded.fasta ${meta.id}.newly_found_junction.fastq

  cutadapt ${meta.id}.newly_found_junction.fastq -g ACAGC --output ${meta.id}.CUTADAPT_five_prime_trim.fastq
    cutadapt ${meta.id}.CUTADAPT_five_prime_trim.fastq -a AACAC --output ${meta.id}.CUTADAPT_three_prime_trim.fastq
    cutadapt ${meta.id}.CUTADAPT_three_prime_trim.fastq -g GTGTT --output ${meta.id}.CUTADAPT_REV_five_prime_trim.fastq
    cutadapt ${meta.id}.CUTADAPT_REV_five_prime_trim.fastq -a GCTGT --output ${meta.id}.REV_three_prime_trim.fastq

    seqtk seq -a ${meta.id}.REV_three_prime_trim.fastq > ${meta.id}.REV_three_prime_trim.fasta
    cat ${meta.id}.REV_three_prime_trim.fasta > ${meta.id}.DRAFT_0_junction_sequences.fasta
  """
}
//  cp ${processReads}/${meta.id}.clean_merged.fastq .

process firstPass_Junction {
  label 'mappingSoftware'
  // label 'mappingSoftware_conda'
  label 'process_medium'
  errorStrategy 'ignore'

  publishDir "${Junc_mapping}", pattern: "*FILTER_ON_HEURISTICS.fasta", mode: "copy"
  publishDir "${Junc_mapping}", pattern: "*DRAFT_1.00001_junction_sequences.fasta", mode: "copy"
  publishDir "${Junc_mapping}", pattern: "*DRAFT_1.00001.consensus.fasta", mode: "copy"

  input:
  path(draft0)

  output:
  path("*FILTER_ON_HEURISTICS.fasta")
  path("*DRAFT_1.00001_junction_sequences.fasta")
  path("*DRAFT_1.00001.consensus.fasta")
  path "*DRAFT_1.00002.consensus.fasta", emit: draft1_00002
  
  script:
  """
  touch ${draft0.simpleName}.DRAFT_1_junction_sequences.fasta
  while read line
        do
      if [[ \${line:0:1} == '>' ]]
          then
              outfile=\${line#>}.fa
              echo \$line > \$outfile
          else
             echo \$line >> \$outfile
         fi
      done < ${draft0}

    ls *.fa | while read my_fastas
      do
       if 

      [ `grep -il 'TGCGGAAACTGTATTTTTA.*TAAAAATTTAGTACACCTAGCC' \$my_fastas | wc -l` -gt 0 ] || [ `grep -il 'GGCTAGGTGTACTAAATTTTTA.*TAAAAATACAGTTTCCGCA' \$my_fastas | wc -l` -gt 0 ] || \

      [ `grep -il 'TTATTTAATTTTACTATTTTAAAT.*ATTTAGTACACC' \$my_fastas | wc -l` -gt 0 ] || [ `grep -il 'GGTGTACTAAAT.*ATTTAAAATAGTAAAATTAAATAA' \$my_fastas | wc -l` -gt 0 ] || \

      [ `grep -il 'GAAACTGTATTTTTA.*TAAAAATTTAGTACACCT' \$my_fastas | wc -l` -gt 0 ] || [ `grep -il 'AGGTGTACTAAATTTTTA.*TAAAAATACAGTTTC' \$my_fastas | wc -l` -gt 0 ] || \

      [ `grep -il 'AATTTTACTATTTTAAAT.*ATTTAGTACA' \$my_fastas | wc -l` -gt 0 ] || [ `grep -il 'TGTACTAAAT.*ATTTAAAATAGTAAAATT' \$my_fastas | wc -l` -gt 0 ]

          then 
            cat \$my_fastas >> ${draft0.simpleName}.DRAFT_1_junction_sequences.fasta
       fi
      done

      cp ${processReads}/${draft0.simpleName}.clean_merged.fastq .

      cd-hit-est -i ${draft0.simpleName}.DRAFT_1_junction_sequences.fasta -o ${draft0.simpleName}.DRAFT_1.00001_junction_sequences.fasta -c 1 -g 1 -d 0 -T ${params.threads} -M ${params.RAM}
      bowtie2-build ${draft0.simpleName}.DRAFT_1.00001_junction_sequences.fasta ${draft0.simpleName}.DRAFT_1.00001_junction_sequences.fasta_BT_INDEX
      bowtie2 -x ${draft0.simpleName}.DRAFT_1.00001_junction_sequences.fasta_BT_INDEX -U ${draft0.simpleName}.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5, --threads ${params.threads} --local > ${draft0.simpleName}.DRAFT_1.OUTPUT.bam

      samtools sort -T -n ${draft0.simpleName}.DRAFT_1.OUTPUT.bam -o ${draft0.simpleName}.DRAFT_1.aln.sorted.bam
      samtools index ${draft0.simpleName}.DRAFT_1.aln.sorted.bam

      bcftools mpileup -Ou -f ${draft0.simpleName}.DRAFT_1_junction_sequences.fasta ${draft0.simpleName}.DRAFT_1.aln.sorted.bam | bcftools call -mv -Oz -o ${draft0.simpleName}.DRAFT_1.calls.vcf.gz
      bcftools index ${draft0.simpleName}.DRAFT_1.calls.vcf.gz

      cat ${draft0.simpleName}.DRAFT_1.00001_junction_sequences.fasta | bcftools consensus ${draft0.simpleName}.DRAFT_1.calls.vcf.gz > ${draft0.simpleName}.DRAFT_1.00001.consensus.fasta
      seqtk seq -F '#' ${draft0.simpleName}.DRAFT_1.00001.consensus.fasta > ${draft0.simpleName}.DRAFT_1.00001_junction_sequences.fastq

     cutadapt ${draft0.simpleName}.DRAFT_1.00001_junction_sequences.fastq -g ACAGC --output ${draft0.simpleName}.CUTADAPT_five_prime_trim.fastq
     cutadapt ${draft0.simpleName}.CUTADAPT_five_prime_trim.fastq -a AACAC  --output ${draft0.simpleName}.CUTADAPT_three_prime_trim.fastq
     cutadapt ${draft0.simpleName}.CUTADAPT_three_prime_trim.fastq -g GTGTT --output ${draft0.simpleName}.CUTADAPT_REV_five_prime_trim.fastq
     cutadapt ${draft0.simpleName}.CUTADAPT_REV_five_prime_trim.fastq -a GCTGT --output ${draft0.simpleName}.REV_three_prime_trim.fastq

     seqtk seq -a ${draft0.simpleName}.REV_three_prime_trim.fastq > ${draft0.simpleName}.REV_three_prime_trim.fasta

     cutadapt ${draft0.simpleName}.REV_three_prime_trim.fasta  -g ACAGC --output ${draft0.simpleName}.REV_three_prime_trim2.fasta
     cutadapt ${draft0.simpleName}.REV_three_prime_trim2.fasta -g ACAGC --output ${draft0.simpleName}.REV_three_prime_trim3.fasta
     cutadapt ${draft0.simpleName}.REV_three_prime_trim3.fasta -g ACAGC --output ${draft0.simpleName}.REV_three_prime_trim4.fasta
     cutadapt ${draft0.simpleName}.REV_three_prime_trim4.fasta -g ACAGC --output ${draft0.simpleName}.REV_three_prime_trim5.fasta
     cutadapt ${draft0.simpleName}.REV_three_prime_trim5.fasta -g ACAGC --output ${draft0.simpleName}.REV_three_prime_trim6.fasta
     cutadapt ${draft0.simpleName}.REV_three_prime_trim6.fasta -g ACAGC --output ${draft0.simpleName}.FILTER_ON_HEURISTICS.fasta

     bowtie2-build ${draft0.simpleName}.FILTER_ON_HEURISTICS.fasta ${draft0.simpleName}.FILTER_ON_HEURISTICS.fasta_BT_INDEX
     bowtie2 -x ${draft0.simpleName}.FILTER_ON_HEURISTICS.fasta_BT_INDEX  -U ${draft0.simpleName}.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads ${params.threads} --local > ${draft0.simpleName}.DRAFT_1.OUTPUT2.bam

     samtools sort -T -n ${draft0.simpleName}.DRAFT_1.OUTPUT2.bam -o ${draft0.simpleName}.DRAFT_1.aln.sorted2.bam
     samtools index ${draft0.simpleName}.DRAFT_1.aln.sorted2.bam

     bcftools mpileup -Ou -f ${draft0.simpleName}.FILTER_ON_HEURISTICS.fasta ${draft0.simpleName}.DRAFT_1.aln.sorted2.bam | bcftools call -mv -Oz -o ${draft0.simpleName}.DRAFT_1.calls2.vcf.gz
     bcftools index ${draft0.simpleName}.DRAFT_1.calls2.vcf.gz

     cat ${draft0.simpleName}.FILTER_ON_HEURISTICS.fasta | bcftools consensus ${draft0.simpleName}.DRAFT_1.calls2.vcf.gz > ${draft0.simpleName}.DRAFT_1.00002.consensus.fasta
      """
}

process firstPass_JunctionCheck1 {
  label 'mappingSoftware'
  label 'process_medium'
  // label 'mappingSoftware_conda'

  input:
  path(inConsensus)

  output:
  path "*DRAFT_1.01_junction_sequences.fasta", emit: firstPass_out, optional: true

  shell:
  '''
  cp !{Junc_mapping}/!{inConsensus.simpleName}.FILTER_ON_HEURISTICS.fasta .
  cp !{Junc_mapping}/!{inConsensus.simpleName}.DRAFT_1.00001.consensus.fasta .
  cp !{Junc_mapping}/!{inConsensus.simpleName}.DRAFT_1.00001_junction_sequences.fasta .
  cp !{Junc_mapping}/fullRefs_junction.fasta JUNCTION_REFS.fasta

  perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' !{inConsensus} > !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta
  printf "\n" >> !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta

  for CHECK_CONSENSUS_DRAFT_1 in `grep -h '>' !{inConsensus} | sed 's/>//g'`
          do
                echo $CHECK_CONSENSUS_DRAFT_1 > reference.txt
                grep -A 1 -wFf reference.txt !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta > X_DRAFT_1_CONSENSUS.fasta
                grep -A 1 -wFf reference.txt !{inConsensus.simpleName}.FILTER_ON_HEURISTICS.fasta > X_DRAFT_1_BEFORE_consensus.fasta
            rm reference.txt

                if [ `cat X_DRAFT_1_CONSENSUS.fasta | tail -1` == `cat X_DRAFT_1_BEFORE_consensus.fasta | tail -1` ];
                    then 
                      sed -i -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta
                else touch !{inConsensus.simpleName}.outCheck1
                fi
          done

  for CHECK_CONSENSUS_DRAFT_1 in `grep -h '>' !{inConsensus.simpleName}.DRAFT_1.00001.consensus.fasta | sed 's/>//g'`
      do
        echo $CHECK_CONSENSUS_DRAFT_1 > reference.txt
        grep -A 1 -wFf reference.txt !{inConsensus.simpleName}.FILTER_ON_HEURISTICS.fasta || true > Y_DRAFT_1_BEFORE_consensus.fasta 
          grep -A 1 -wFf reference.txt !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta || true >  Y_DRAFT_1_CONSENSUS.fasta   

          if [ `cat Y_DRAFT_1_CONSENSUS.fasta | tail -1 | grep "n" | wc -l` -gt 0 ];
              then 
                sed -i -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta
        fi
        
        if [ `cat Y_DRAFT_1_BEFORE_consensus.fasta | tail -1 | grep "n" | wc -l` -gt 0 ];
              then 
                sed -i -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta
        fi

        if [ `sed -r -n -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 p" !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta | tail -1 | wc -l` -eq 0 ] && \
            [ `sed -r -n -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 p" !{inConsensus.simpleName}.FILTER_ON_HEURISTICS.fasta  | tail -1 | grep "n" | wc -l` -gt 0 ];
              then
                sed -i -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" !{inConsensus.simpleName}.FILTER_ON_HEURISTICS.fasta 
        fi

        if [ `sed -r -n -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 p" !{inConsensus.simpleName}.FILTER_ON_HEURISTICS.fasta | tail -1 | grep "n" | wc -l` -gt 0 ] && \
            [ `sed -r -n -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 p" !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta | tail -1 | grep "n" | wc -l` -gt 0 ];
         then 
          sed -i -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" !{inConsensus.simpleName}.FILTER_ON_HEURISTICS.fasta
          sed -i -e "/$CHECK_CONSENSUS_DRAFT_1/,+1 d" !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta
        fi  
      done

  makeblastdb -in JUNCTION_REFS.fasta -input_type fasta -dbtype nucl

  for SCND_CHECK_CONSENSUS_DRAFT_1 in `grep -h '>' !{inConsensus.simpleName}.DRAFT_1.00001_junction_sequences.fasta | sed 's/>//g'`
    do
            echo $SCND_CHECK_CONSENSUS_DRAFT_1 > reference.txt
            grep -A 1 -wFf reference.txt !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta || true > X_DRAFT_1_CONSENSUS.fasta
            grep -A 1 -wFf reference.txt !{inConsensus.simpleName}.FILTER_ON_HEURISTICS.fasta || true > X_DRAFT_1_BEFORE_consensus.fasta

            blastn -db JUNCTION_REFS.fasta \
            -query X_DRAFT_1_CONSENSUS.fasta \
            -word_size 7 \
            -evalue 0.001 \
            -perc_identity 100 \
            -qcov_hsp_perc 100 \
            -num_threads !{params.threads} \
            -out X_DRAFT_1_CONSENSUS.blast_result \
            -max_target_seqs 1 \
            -dust no \
            -soft_masking false \
            -outfmt "6 qseqid qseq"


            blastn -db JUNCTION_REFS.fasta \
            -query X_DRAFT_1_BEFORE_consensus.fasta \
            -word_size 7 \
            -evalue 0.001 \
            -perc_identity 100 \
            -qcov_hsp_perc 100 \
            -num_threads !{params.threads} \
            -out X_DRAFT_1_BEFORE_consensus.blast_result \
            -max_target_seqs 1 \
            -dust no \
            -soft_masking false \
            -outfmt "6 qseqid qseq"

        if [ `cat X_DRAFT_1_CONSENSUS.blast_result | wc -l` == `cat X_DRAFT_1_BEFORE_consensus.blast_result | wc -l` ];
             then 
              sed -i -e  "/$SCND_CHECK_CONSENSUS_DRAFT_1/,+1 d" !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta
            fi

            if [ `cat X_DRAFT_1_CONSENSUS.blast_result | wc -l` -eq 0 ] && [ `cat X_DRAFT_1_BEFORE_consensus.blast_result | wc -l` -gt 0 ];
             then 
              sed -i -e "/$SCND_CHECK_CONSENSUS_DRAFT_1/,+1 d" !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta
            fi

            if [ `cat X_DRAFT_1_CONSENSUS.blast_result | wc -l` -eq 0 ] && [ `cat X_DRAFT_1_BEFORE_consensus.blast_result | wc -l` -eq 0 ];
             then 
              sed -i -e "/$SCND_CHECK_CONSENSUS_DRAFT_1/,+1 d" !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta
            fi

            if [ `cat X_DRAFT_1_BEFORE_consensus.blast_result | wc -l` -eq 0 ] && [ `cat X_DRAFT_1_CONSENSUS.blast_result | wc -l` -gt 0 ];
             then 
              sed -i -e "/$SCND_CHECK_CONSENSUS_DRAFT_1/,+1 d" !{inConsensus.simpleName}.FILTER_ON_HEURISTICS.fasta
            fi
         done


      cat !{inConsensus.simpleName}.ONE_LINE.DRAFT_1.00002.consensus.fasta >> !{inConsensus.simpleName}.DRAFT_1.01_junction_sequences.fasta
      cat !{inConsensus.simpleName}.FILTER_ON_HEURISTICS.fasta >> !{inConsensus.simpleName}.DRAFT_1.01_junction_sequences.fasta
  '''
}

process clusterFirstPass {
  label 'mappingSoftware'
  label 'process_medium'
  // label 'mappingSoftware_conda'

  input:
  path(draft1_1)

  output:
  path "*DRAFT_2.consensus.fasta", emit: draft2_consensus, optional: true

  shell:
  '''
  cp !{processReads}/!{draft1_1.simpleName}.clean_merged.fastq .

  cd-hit-est -i !{draft1_1} -o !{draft1_1.simpleName}.DRAFT_1.02_junction_sequences.fasta -c 1 -g 1 -d 0 -T !{params.threads} -M !{params.RAM}
  seqtk seq -F '#' !{draft1_1.simpleName}.DRAFT_1.02_junction_sequences.fasta > !{draft1_1.simpleName}.DRAFT_1.02_junction_sequences.fastq

  cutadapt !{draft1_1.simpleName}.DRAFT_1.02_junction_sequences.fastq -g ACAGC --output !{draft1_1.simpleName}.CUTADAPT_five_prime_trim.fastq
    cutadapt !{draft1_1.simpleName}.CUTADAPT_five_prime_trim.fastq -a AACAC  --output !{draft1_1.simpleName}.CUTADAPT_three_prime_trim.fastq
    cutadapt !{draft1_1.simpleName}.CUTADAPT_three_prime_trim.fastq -g GTGTT --output !{draft1_1.simpleName}.CUTADAPT_REV_five_prime_trim.fastq
    cutadapt !{draft1_1.simpleName}.CUTADAPT_REV_five_prime_trim.fastq -a GCTGT --output !{draft1_1.simpleName}.DRAFT_1.1_junction_sequences.fastq 

    seqtk seq -a !{draft1_1.simpleName}.DRAFT_1.1_junction_sequences.fastq > !{draft1_1.simpleName}.DRAFT_1.1_junction_sequences.fasta

    cd-hit-est -i !{draft1_1.simpleName}.DRAFT_1.1_junction_sequences.fasta -o !{draft1_1.simpleName}.DRAFT_1.11_junction_sequences.fasta -c 1 -g 1 -d 0 -T !{params.threads} -M !{params.RAM}
    bowtie2-build !{draft1_1.simpleName}.DRAFT_1.11_junction_sequences.fasta !{draft1_1.simpleName}.DRAFT_1.11_junction_sequences.fasta_BT_INDEX
    bowtie2 -x !{draft1_1.simpleName}.DRAFT_1.11_junction_sequences.fasta_BT_INDEX -U !{draft1_1.simpleName}.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads !{params.threads} --local > !{draft1_1.simpleName}.DRAFT_2.OUTPUT.bam

    samtools sort -T -n !{draft1_1.simpleName}.DRAFT_2.OUTPUT.bam -o !{draft1_1.simpleName}.DRAFT_2.aln.sorted.bam
    samtools index !{draft1_1.simpleName}.DRAFT_2.aln.sorted.bam
    samtools depth -aa -d 0 !{draft1_1.simpleName}.DRAFT_2.aln.sorted.bam > !{draft1_1.simpleName}.DRAFT_2.coverage 

    touch !{draft1_1.simpleName}.DRAFT_2.contigs_with_coverage_less_than.10.txt

    for JUNK_REMOVAL in `cat  !{draft1_1.simpleName}.DRAFT_2.coverage | cut -f 1 | uniq`
        do
         if [ `cat  !{draft1_1.simpleName}.DRAFT_2.coverage | sed 's/^/>/' | sed "s/$JUNK_REMOVAL\t/$JUNK_REMOVAL\t/g" | grep "$JUNK_REMOVAL" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt 10 ]
          then 
          echo $JUNK_REMOVAL >> !{draft1_1.simpleName}.DRAFT_2.contigs_with_coverage_less_than.10.txt
         fi
        done  

    cat !{draft1_1.simpleName}.DRAFT_1.11_junction_sequences.fasta > !{draft1_1.simpleName}.DRAFT_2.consensus.fasta
    sed -i -e 's/^/>/' !{draft1_1.simpleName}.DRAFT_2.contigs_with_coverage_less_than.10.txt

    for HAPS_TO_REMOVE in `cat !{draft1_1.simpleName}.DRAFT_2.contigs_with_coverage_less_than.10.txt`
      do
         sed -i -e "/$HAPS_TO_REMOVE/,+1 d" !{draft1_1.simpleName}.DRAFT_2.consensus.fasta
      done
  '''
}

process secondPass_junction {
  label 'mappingSoftware'
  label 'process_medium'
  // label 'mappingSoftware_conda'

  publishDir "${Junc_mapping}", pattern: "*DRAFT_2.200.consensus.fasta", mode: "copy"
  errorStrategy 'ignore'

  input:
  path(draft2)

  output:
  path "*DRAFT_2.200.consensus.fasta", emit: draft200_out, optional:true
  path("*DRAFT_2.200.consensus.fasta")

  shell:
  '''
  cp !{processReads}/!{draft2.simpleName}.clean_merged.fastq .
  cp !{Junc_mapping}/fullRefs_junction.fasta JUNCTION_REFS.fasta

    bowtie2-build !{draft2} !{draft2.simpleName}.DRAFT_2_consensus.fasta_BT_INDEX
    bowtie2 -x !{draft2.simpleName}.DRAFT_2_consensus.fasta_BT_INDEX -U !{draft2.simpleName}.clean_merged.fastq -q --threads !{params.threads} --local > !{draft2.simpleName}.DRAFT_2_RELAXED.OUTPUT.bam

    samtools sort -T -n !{draft2.simpleName}.DRAFT_2_RELAXED.OUTPUT.bam -o !{draft2.simpleName}.DRAFT_2_RELAXED_sorted.bam
    samtools index !{draft2.simpleName}.DRAFT_2_RELAXED_sorted.bam

    bcftools mpileup -Ou -f !{draft2} !{draft2.simpleName}.DRAFT_2_RELAXED_sorted.bam | bcftools call -mv -Oz  -o !{draft2.simpleName}.DRAFT_2.calls3.vcf.gz
    bcftools index !{draft2.simpleName}.DRAFT_2.calls3.vcf.gz
    cat !{draft2}| bcftools consensus !{draft2.simpleName}.DRAFT_2.calls3.vcf.gz > !{draft2.simpleName}.DRAFT_2.RELAXED.consensus.fasta

    perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' !{draft2.simpleName}.DRAFT_2.RELAXED.consensus.fasta > !{draft2.simpleName}.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta
  printf "\n" >> !{draft2.simpleName}.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta

  makeblastdb -in JUNCTION_REFS.fasta -input_type fasta -dbtype nucl

    for SCND_CHECK_CONSENSUS_DRAFT_2 in `grep -h '>' !{draft2.simpleName}.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta | sed 's/>//g'`
          do
              echo $SCND_CHECK_CONSENSUS_DRAFT_2 > reference.txt

              grep -A 1 -wFf reference.txt !{draft2.simpleName}.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta > Z_DRAFT_2_CONSENSUS.fasta
              grep -A 1 -wFf reference.txt !{draft2} > Z_DRAFT_2_BEFORE_consensus.fasta

              blastn -db JUNCTION_REFS.fasta \
              -query Z_DRAFT_2_CONSENSUS.fasta \
              -word_size 7 \
              -evalue 0.001 \
              -perc_identity 100 \
              -qcov_hsp_perc 100 \
              -num_threads !{params.threads} \
              -out Z_DRAFT_2_CONSENSUS.blast_result \
              -max_target_seqs 1 \
              -dust no \
              -soft_masking false \
              -outfmt "6 qseqid qseq"


              blastn -db JUNCTION_REFS.fasta \
              -query Z_DRAFT_2_BEFORE_consensus.fasta \
              -word_size 7 \
              -evalue 0.001 \
              -perc_identity 100 \
              -qcov_hsp_perc 100 \
              -num_threads !{params.threads} \
              -out Z_DRAFT_2_BEFORE_consensus.blast_result \
              -max_target_seqs 1 \
              -dust no \
              -soft_masking false \
              -outfmt "6 qseqid qseq"



              if [ `cat Z_DRAFT_2_CONSENSUS.blast_result | wc -l` == `cat Z_DRAFT_2_BEFORE_consensus.blast_result | wc -l` ];
                    then 
                        sed -i -e "/$SCND_CHECK_CONSENSUS_DRAFT_2/,+1 d" !{draft2.simpleName}.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta
              fi

              if [ `cat Z_DRAFT_2_CONSENSUS.blast_result | wc -l` -eq 0 ] && [ `cat Z_DRAFT_2_BEFORE_consensus.blast_result | wc -l` -gt 0 ];
                  then 
                      sed -i "/$SCND_CHECK_CONSENSUS_DRAFT_2/,+1 d" !{draft2.simpleName}.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta
              fi

              if [ `cat Z_DRAFT_2_CONSENSUS.blast_result | wc -l` -eq 0 ] && [ `cat Z_DRAFT_2_BEFORE_consensus.blast_result | wc -l` -eq 0 ];
                    then 
                      sed -i -e "/$SCND_CHECK_CONSENSUS_DRAFT_2/,+1 d" !{draft2.simpleName}.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta
              fi

              if [ `cat Z_DRAFT_2_BEFORE_consensus.blast_result | wc -l` -eq 0 ] && [ `cat Z_DRAFT_2_CONSENSUS.blast_result | wc -l` -gt 0 ];
                    then 
                      sed -i -e "/$SCND_CHECK_CONSENSUS_DRAFT_2/,+1 d" !{draft2}
              fi
         done

    cat !{draft2.simpleName}.DRAFT_2_RELAXED_ONE_LINE.consensus.fasta >> !{draft2.simpleName}.AFTER_RELAX.fasta
    cat !{draft2} >> !{draft2.simpleName}.AFTER_RELAX.fasta
    cat !{draft2.simpleName}.AFTER_RELAX.fasta > !{draft2.simpleName}.DRAFT_2.200.consensus.fasta
  '''
}

process thirdPass_junction{
  label 'mappingSoftware'
  label 'process_medium'
  // label 'mappingSoftware_conda'

  errorStrategy 'ignore'

 publishDir "${Junc_mapping}", pattern: "*contigs_with_coverage_less_than*", mode: "copy"
 publishDir "${Junc_mapping}", pattern: "*DRAFT_3.consensus.fasta", mode: "copy"
 publishDir "${Junc_mapping}", pattern: "*DRAFT_3.CLUSTERS_consensus.fasta", mode: "copy"
 publishDir "${Junc_mapping}", pattern: "*DRAFT_3.coverage", mode: "copy"
 publishDir "${Junc_mapping}", pattern: "*sumCov", mode: "copy"
 publishDir "${Junc_mapping}", pattern: "*totalCov", mode: "copy"
 publishDir "${Junc_mapping}", pattern: "*DRAFT_4*", mode: "copy"
  
  input:
  path(draft200)

  output:
  path "*FINAL.validated_junction_sequences.fasta", emit: validatedJunctions, optional: true
  path("*contigs_with_coverage_less_than*")
  path("*DRAFT_3.consensus.fasta")
  path("*DRAFT_3.CLUSTERS_consensus.fasta")
  path("*DRAFT_3.coverage")
  path("*sumCov")
  path("*totalSum")
  path("*DRAFT_4*")

  shell:
  '''
  cp !{processReads}/!{draft200.simpleName}.clean_merged.fastq .
  cp !{Junc_mapping}/!{draft200.simpleName}.READS_MAPPING.txt .

  bowtie2-build !{draft200} !{draft200.simpleName}.DRAFT_2.200.consensus.fasta_BT_INDEX
  bowtie2 -x !{draft200.simpleName}.DRAFT_2.200.consensus.fasta_BT_INDEX -U !{draft200.simpleName}.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads !{params.threads} --local > !{draft200.simpleName}.DRAFT_3.OUTPUT.bam

  samtools sort -T -n !{draft200.simpleName}.DRAFT_3.OUTPUT.bam -o !{draft200.simpleName}.DRAFT_3.aln.sorted.bam
  samtools index !{draft200.simpleName}.DRAFT_3.aln.sorted.bam
  samtools depth -aa -d 0 !{draft200.simpleName}.DRAFT_3.aln.sorted.bam > !{draft200.simpleName}.DRAFT_3.coverage # exports a file containing the coverage.

  awk '{print$1}' !{draft200.simpleName}.DRAFT_3.coverage | sort | uniq | while read line; do cat !{draft200.simpleName}.DRAFT_3.coverage | sed 's/^/>/' | grep "$line" | sort -u -nrk 3n | head -1 ; done > !{draft200.simpleName}.DRAFT_3.sumCov

  awk '{print$3}' !{draft200.simpleName}.DRAFT_3.sumCov | paste -sd+ | bc > !{draft200.simpleName}.totalSum

  DRAFT_3_Total=`cat !{draft200.simpleName}.totalSum | awk '{print$1}'`
  JUNCTION_rough_coverage_cutoff_test=`echo "scale=0; ($DRAFT_3_Total*0.15)/1" | bc`

  echo $JUNCTION_rough_coverage_cutoff_test > aFile_tester

  READS_MAPPING_TO_YOUR_JUNCTION=`cat !{draft200.simpleName}.READS_MAPPING.txt | awk '{print$1}'`
  JUNCTION_rough_coverage_cutoff=`echo "scale=0; ($READS_MAPPING_TO_YOUR_JUNCTION*0.010)/1" | bc`
  JUNCTION_rough_coverage_cutoff_2=`echo "scale=0; ($READS_MAPPING_TO_YOUR_JUNCTION*0.015)/1" | bc`

  touch !{draft200.simpleName}.contigs_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff.txt

  

  for JUNK_REMOVAL in `cat !{draft200.simpleName}.DRAFT_3.coverage | cut -f 1 | uniq`
      do
      if [ `cat !{draft200.simpleName}.DRAFT_3.coverage | sed 's/^/>/' | sed "s/$JUNK_REMOVAL\t/$JUNK_REMOVAL\t/g" | grep "$JUNK_REMOVAL" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt $JUNCTION_rough_coverage_cutoff ]
       then 
        echo $JUNK_REMOVAL >> !{draft200.simpleName}.contigs_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff.txt
      fi
      done

  sed -i -e 's/^/>/' !{draft200.simpleName}.contigs_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff.txt
  cat !{draft200} > !{draft200.simpleName}.DRAFT_3.consensus.fasta

  for HAPS_TO_REMOVE in `cat !{draft200.simpleName}.contigs_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff.txt`
      do
      sed -i -e "/$HAPS_TO_REMOVE/,+1 d" !{draft200.simpleName}.DRAFT_3.consensus.fasta
      done

  cd-hit-est -i !{draft200.simpleName}.DRAFT_3.consensus.fasta -o !{draft200.simpleName}.DRAFT_3.CLUSTERS_consensus.fasta -c 1 -g 1 -d 0 -T !{params.threads} -M !{params.RAM}
  bowtie2-build !{draft200.simpleName}.DRAFT_3.CLUSTERS_consensus.fasta !{draft200.simpleName}.DRAFT_3.CLUSTERS_consensus.fasta_BT_INDEX
  bowtie2 -x !{draft200.simpleName}.DRAFT_3.CLUSTERS_consensus.fasta_BT_INDEX -U !{draft200.simpleName}.clean_merged.fastq -q -D 20 -R 3 -N 0 -L 32 -i S,2,5 --threads !{params.threads} --local > !{draft200.simpleName}.DRAFT_4.OUTPUT.bam

    samtools sort -T -n !{draft200.simpleName}.DRAFT_4.OUTPUT.bam -o !{draft200.simpleName}.DRAFT_4.aln.sorted.bam
    samtools index !{draft200.simpleName}.DRAFT_4.aln.sorted.bam
    samtools depth -aa -d 0 !{draft200.simpleName}.DRAFT_4.aln.sorted.bam > !{draft200.simpleName}.DRAFT_4.coverage

  awk '{print$1}' !{draft200.simpleName}.DRAFT_4.coverage | sort | uniq | while read line; do cat !{draft200.simpleName}.DRAFT_4.coverage | sed 's/^/>/' | grep "$line" | sort -u -nrk 3n | head -1 ; done > !{draft200.simpleName}.DRAFT_4.sumCov

  awk '{print$3}' !{draft200.simpleName}.DRAFT_4.sumCov | paste -sd+ | bc > !{draft200.simpleName}.DRAFT_4.totalSum

  DRAFT_4_Total=`cat !{draft200.simpleName}.DRAFT_4.totalSum | awk '{print$1}'`
  JUNCTION_rough_coverage_cutoff_DRAFT_4=`echo "scale=0; ($DRAFT_4_Total*0.15)/1" | bc`
  
    touch !{draft200.simpleName}.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.!{params.depthNewHap_min}.txt

    touch !{draft200.simpleName}.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.!{params.depthNewHap_min}.and.$JUNCTION_rough_coverage_cutoff_DRAFT_4.txt

    for NEW_HAPS in `cat !{draft200.simpleName}.DRAFT_4.coverage | cut -f 1 | uniq` 
        do
      if [ `cat !{draft200.simpleName}.DRAFT_4.coverage | sed 's/^/>/' | sed "s/$NEW_HAPS\t/$NEW_HAPS._\t/g" | grep "$NEW_HAPS._" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt !{params.depthNewHap_min} ] \
        && [ `cat !{draft200.simpleName}.DRAFT_4.coverage | sed 's/^/>/' | sed "s/$NEW_HAPS\t/$NEW_HAPS._\t/g" | grep "$NEW_HAPS._" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt $JUNCTION_rough_coverage_cutoff_2 ];
       then 
        echo $NEW_HAPS >> !{draft200.simpleName}.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.!{params.depthNewHap_min}.and.$JUNCTION_rough_coverage_cutoff_DRAFT_4.txt
      fi
        done
   
  for NEW_HAPS in `cat !{draft200.simpleName}.DRAFT_4.coverage | cut -f 1 | uniq` 
    do 
    if  [ `cat !{draft200.simpleName}.DRAFT_4.coverage | sed 's/^/>/' | sed "s/$NEW_HAPS\t/$NEW_HAPS._\t/g" | grep "$NEW_HAPS._" | sort -u -nrk 3n | head -1 | awk '{print $3}'` -lt $JUNCTION_rough_coverage_cutoff_DRAFT_4 ];
    then
    echo $NEW_HAPS >> !{draft200.simpleName}.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.!{params.depthNewHap_min}.and.$JUNCTION_rough_coverage_cutoff_DRAFT_4.txt
    fi
    done

  cat !{draft200.simpleName}.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.!{params.depthNewHap_min}.and.$JUNCTION_rough_coverage_cutoff_DRAFT_4.txt | sort | uniq > !{draft200.simpleName}.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.!{params.depthNewHap_min}.and.$JUNCTION_rough_coverage_cutoff_DRAFT_4.final.txt

   sed -i -e 's/^/>/' !{draft200.simpleName}.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.!{params.depthNewHap_min}.and.$JUNCTION_rough_coverage_cutoff_DRAFT_4.final.txt
   cp !{draft200.simpleName}.DRAFT_3.CLUSTERS_consensus.fasta !{draft200.simpleName}.DRAFT_4.validated_junction_sequences.fasta 

   
    for HAPS_TO_REMOVE in `cat !{draft200.simpleName}.DRAFT_4.clusters_with_coverage_less_than.$JUNCTION_rough_coverage_cutoff_2.and.!{params.depthNewHap_min}.and.$JUNCTION_rough_coverage_cutoff_DRAFT_4.final.txt`
        do
        sed -i -e "/$HAPS_TO_REMOVE/,+1 d" !{draft200.simpleName}.DRAFT_4.validated_junction_sequences.fasta 
        done

    for JUNCTION_TESTER in `cat !{draft200.simpleName}.DRAFT_4.coverage | cut -f 1 | uniq` 
        do
        coverage_at_each_base=`cat !{draft200.simpleName}.DRAFT_4.coverage | sed -n "/$JUNCTION_TESTER/p" | cut -f 3`
        echo $coverage_at_each_base > $JUNCTION_TESTER.coverage_at_each_base
      number_of_bases=`echo $coverage_at_each_base | wc -w`
        sum_of_bases=`echo $coverage_at_each_base | sed -e 's/ /+/g' | bc`
          average_coverage=`echo "($sum_of_bases/$number_of_bases)" | bc`
      echo $average_coverage > $JUNCTION_TESTER.average_coverage
      standardDeviation=$(
            echo "$coverage_at_each_base" |
              awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)^2)}'
                         )
        echo $standardDeviation > $JUNCTION_TESTER.SD 
        two_SDs=`echo "(2*$standardDeviation)" | bc`
        echo $two_SDs > $JUNCTION_TESTER.2SDs
        rounded_2SDs=`echo ${two_SDs%%.*}`
        echo $rounded_2SDs > $JUNCTION_TESTER.Rounded_2SDs

        if [ $rounded_2SDs -gt $average_coverage ];
          then 
            sed -i -e  "/$JUNCTION_TESTER/,+1 d" !{draft200.simpleName}.DRAFT_4.validated_junction_sequences.fasta
        fi
         done

     cat !{draft200.simpleName}.DRAFT_4.validated_junction_sequences.fasta > !{draft200.simpleName}.FINAL.validated_junction_sequences.fasta
  '''
}

process newJunction_finder {
  label 'mappingSoftware'
  label 'process_medium'
  // label 'mappingSoftware_conda'

  input:
  path(junctionSeqs)

  output:
  path 'allNewJunctionOut.file', emit: allNewJuncHaps, optional: true

  shell:
  '''
  cp !{Junc_mapping}/fullRefs_junction.fasta JUNCTION_REFS.fasta

    while read line
         do
            if [[ ${line:0:1} == '>' ]]
              then
                 outfile=${line#>}.fa
                 echo $line > $outfile
         else
             echo $line >> $outfile
             fi
         done < !{junctionSeqs}

      makeblastdb  -in JUNCTION_REFS.fasta -input_type fasta -dbtype nucl -title junction_check

     for FOOBAR in `ls *.fa`
        do 
          junction_count=`grep 'Junction_Hap' JUNCTION_REFS.fasta | wc -l`
          new_number=`echo "$junction_count + 1" | bc`
          cat $FOOBAR | awk '/^>/{print ">Junction_Hap_'$new_number'"; next}{print}' > !{junctionSeqs.simpleName}.junction_recount_$FOOBAR

          blastn -db JUNCTION_REFS.fasta \
          -query !{junctionSeqs.simpleName}.junction_recount_$FOOBAR \
          -word_size 7 \
          -evalue 0.001 \
          -perc_identity 100 \
          -qcov_hsp_perc 100 \
          -num_threads !{params.threads} \
          -out !{junctionSeqs.simpleName}.RESULT_$FOOBAR.blast_result \
          -max_target_seqs 1 \
          -dust no \
          -soft_masking false \
          -outfmt "6 qseqid pident mismatch gapopen gaps sseqid sseq evalue bitscore"

          match_present=`cat !{junctionSeqs.simpleName}.RESULT_$FOOBAR.blast_result | wc -l`

            if [ $match_present == 0 ];
           
               then

               repeat_length=`cat !{junctionSeqs.simpleName}.junction_recount_$FOOBAR | sed "1d" | tr -cd '[:alpha:]' | wc -m`
               rep_length_2=`echo "$repeat_length + 21 + 22" | bc` 

               sed -i -e "1s/.*/>Mt_Cmt$rep_length_2.X_Junction_Hap_$new_number/" !{junctionSeqs.simpleName}.junction_recount_$FOOBAR
               cat !{junctionSeqs.simpleName}.junction_recount_$FOOBAR >> \
               !{junctionSeqs.simpleName}.Mt_Cmt$rep_length_2.X_Junction_Hap_$new_number.fasta

               theSpecimen=`echo !{junctionSeqs.simpleName}`
               theMarker=`echo Mt_Cmt$rep_length_2.X_Junction`
               theSequence=`cat !{junctionSeqs.simpleName}.junction_recount_$FOOBAR | tail -1`

               echo $theSpecimen $theMarker $theSequence > allNewJunctionOut.file
            fi
        done
  '''
}

process nameNewJuncHaps {
  //  label 'pythonEnv'
 label 'pythonEnv_wf1'
 label 'process_low'
 
  input:
  path(allNewHaps)

  output:
  path "newJuncHaps.fasta", emit: newJuncHaps, optional: true
  
  script:
  """
  cp ${Junc_mapping}/fullRefs_junction.fasta JUNCTION_REFS.fasta
  cat ${allNewHaps} > concatHaps.txt
  grep "^>" JUNCTION_REFS.fasta > existingHapNames.txt
  python3 "$projectDir/bin/newHaps_Junction_names.py" -n concatHaps.txt -e existingHapNames.txt -o "${NEW_HAPS_DIR_Junc}" 
  touch newJuncHaps.fasta
  cat ${NEW_HAPS_DIR_Junc}/*.fasta >>  newJuncHaps.fasta
  """
} 
process makeRefs_assignHaps_junc2 {
  label 'process_low'
  input:
  path(newHaps)

  output:
  path 'refHaps_junc.fasta', emit: completeFile_junc

  script:
  """
  cat ${newHaps} > refHaps_junc.fasta
  """
}


process finalNewRefs_assignHaps_junc2 {
  label 'process_low'
  publishDir "${refHapsDir}", pattern: "hapRefs_junc.fasta", mode: "copy"

  input:
  path(allNew)
  // path(originalNew)

  output:
  path('hapRefs_junc.fasta')
  path 'hapRefs_junc.fasta', emit: outFinalRefs

  script:
  """
  cat ${allNew}  ${params.origNewHaps}/*Junction*.fasta ${params.originalHaplotypes_refs} > hapRefs_junc.fasta

  """
}


process juncBlast {
  label 'mappingSoftware'
  label 'process_medium'
  // label 'mappingSoftware_conda'
  

  publishDir "${params.tmpGenotypes}", pattern: "*GENOTYPE_withJunction", mode: "copy"
 
  input:
  each(juncFile)
  path(refDB)

  output:
  path "*GENOTYPE_withJunction", emit: genotypeJunction, optional: true
  // tuple val("${juncFile.simpleName}"), path("${juncFile.simpleName}.GENOTYPE_withJunction"), emit: genotypeJunction_tuple, optional: true
  path("*GENOTYPE_withJunction")
  path("endMarker"), emit: junction_hapsDone

  script:
  """
  juncLen="\$(wc -l ${juncFile} | awk '{print\$1}')"

  if [ \$juncLen -gt 0 ]
  then 
  makeblastdb -in ${juncFile} -input_type fasta -dbtype nucl

  blastn -db ${juncFile} -query ${refDB} -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -num_threads ${params.threads} -out ${juncFile.simpleName}.GENOTYPE_withJunction -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid pident mismatch gapopen gaps sseq evalue bitscore"

  else
  touch ${juncFile.simpleName}.GENOTYPE_withJunction
  fi

  echo "done" >> endMarker
  """
}
