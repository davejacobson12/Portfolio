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

process BEDEXT_NEWHAPS {
    tag "$meta.id"
    label 'process_medium'
    label 'mappingSoftware'
    // shell '/bin/bash', '-euo', 'pipefail'
    // errorStrategy 'ignore'

    input:
    tuple val(meta), path(sorted), path(index)
    each(marker)
    each(refDB)

    output:
    // tuple val(meta), path("*mapped_only.fastq")      , emit: renamed     , optional:false
    tuple val(meta), path("*assignHaps_coverageCutoff"), emit: markerCov, optional: true
    // path("*.f*")
    // path "*multi-seq",  emit: multiOut, optional: true
    // path("*mapped_only.fastq"), emit: finalClipped_fastq
    // tuple val(meta), path("*clusterCounted.csv"), path("*multi-seq"), emit:bedTuple, optional: true
    tuple val(meta), path("*mapped_only.fastq"), emit: tryBedFastq, optional:true
    path "allNewHapsOut.file", emit: allNewHaps, optional: true

    when:
    task.ext.when == null || task.ext.when


    script:
    // def args = task.ext.args ?: ''
    // def args2 = task.ext.args2 ?: ''
    // def prefix = task.ext.prefix ?: '!{meta.id}'
    
    """
    myName="\$(echo $marker)"

    shortName="\$(echo \$myName | sed 's/_PART_*.//g')"
    grep \$myName ${params.bedFile} > \$myName.myMatch.bed
    left_trim="\$(cat \$myName.myMatch.bed | cut -f2)"
    right_trim="\$(cat \$myName.myMatch.bed | cut -f3)"

    echo \$left_trim \$right_trim > \$myName.myMatch.trim
    echo \$shortName > \$shortName.markerNames


    samtools view -bu -F 4 ${sorted} \$shortName | java -jar ${params.samjs} -e "record.alignmentStart <= \$left_trim && record.alignmentEnd >= \$right_trim" --samoutputformat BAM -o ${meta.id}.\$myName.bam

    echo "after first step" 

    java -jar ${params.pcrclipreads} --bed \$myName.myMatch.bed ${meta.id}.\$myName.bam --samoutputformat BAM -o ${meta.id}.\$myName.softClipped.bam

    java -jar ${params.biostar84452}  ${meta.id}.\$myName.softClipped.bam > ${meta.id}.\$myName.finalClipped.map_to_marker_only.bam

    samtools view -h -F 4 ${meta.id}.\$myName.finalClipped.map_to_marker_only.bam| samtools view -bS > ${meta.id}.\$myName.mapped_only.sam

    samtools fastq ${meta.id}.\$myName.mapped_only.sam > ${meta.id}.\$myName.mapped_only.fastq


    samtools sort -T -n ${meta.id}.\$myName.finalClipped.map_to_marker_only.bam -o ${meta.id}.\$myName.finalClipped.map_to_marker_only_aln_sorted.bam

    samtools index ${meta.id}.\$myName.finalClipped.map_to_marker_only_aln_sorted.bam

    samtools depth -aa -d 0 ${meta.id}.\$myName.finalClipped.map_to_marker_only_aln_sorted.bam | awk "/^\$shortName/" > firstTry1.txt
    
    head -\$right_trim firstTry1.txt > firstTry2.txt
    awk "NR > \$left_trim" firstTry2.txt > ${meta.id}.\$myName.coverage

    avgCov="\$(awk '{ total += \$3; count++ } END { print total/count }' ${meta.id}.\$myName.coverage)"
    echo "after avgCov"
    echo "scale=0; (\$avgCov * 0.10)/1" | bc > ${meta.id}.\$myName.assignHaps_coverageCutoff
    echo "scale=0; (\$avgCov * 0.25)/1" | bc > ${meta.id}.\$myName.newHaps_coverageCutoff
    newHaps_cutoff=`echo "scale=0; (\$avgCov * 0.25)/1" | bc`

    newHaps2_cutoff='\$(cat ${meta.id}.\$myName.newHaps_coverageCutoff)'

    echo \$newHaps_cutoff > ${meta.id}.\$myName.newHapsTest_coverageCutoff

    cd-hit-est -i ${meta.id}.\$myName.mapped_only.fastq -o ${meta.id}.\$myName.clean_merged_Clusters.fq -c 1 -g 1 -d 0 -T ${params.max_cpus} -s 1


    clus="\$(wc -l ${meta.id}.\$myName.clean_merged_Clusters.fq | awk '{print\$1}')"

    if [ \$clus -gt 0 ]
    then
    seqret -sequence ${meta.id}.\$myName.clean_merged_Clusters.fq -outseq ${meta.id}.\$myName.clean_merged_Clusters.fasta

    seqret -sequence ${meta.id}.\$myName.mapped_only.fastq -outseq ${meta.id}.\$myName.mapped_only.fasta

    perl ${params.multiSeq_pl} ${meta.id}.\$myName.clean_merged_Clusters.fasta ${meta.id}.\$myName.clean_merged_Clusters.fq.clstr ${meta.id}.\$myName.newHaps_multi-seq \$newHaps_cutoff

    cat ${meta.id}.\$myName.clean_merged_Clusters.fq.clstr  | sed 's/Cluster /Cluster/g' > ${meta.id}.\$myName.clean_merged_Clusters_toCount.fq.clstr
    python3 $projectDir/bin/cluster_counter_v2.py -x ${meta.id}.\$myName.clean_merged_Clusters_toCount.fq.clstr -y ${meta.id}.\$myName.clusterCounted.csv
    
    fi


    ls ${meta.id}.\$myName.newHaps_multi-seq  | grep -v ".fasta" > ${meta.id}.\$myName.newHaps_multi-seq.list


    cat ${meta.id}.\$myName.newHaps_multi-seq.list | while read i
    do

    clusDepth="\$(grep "Cluster\$i"  ${meta.id}.\$myName.clusterCounted.csv | awk -F ',' '{print\$2}')"  

    echo \$clusDepth > myDepth_\$i


    cat ${meta.id}.\$myName.newHaps_multi-seq/\$i > ${meta.id}.\$myName.newHaps_multi-seq/\$i.fasta


    makeblastdb -in ${meta.id}.\$myName.newHaps_multi-seq/\$i.fasta -input_type fasta -dbtype nucl

    blastn -db ${meta.id}.\$myName.newHaps_multi-seq/\$i.fasta -query ${refDB} -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -num_threads 4 -out ${meta.id}.\$myName.\$i.singleMatch_blastOut -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid"

    match_present="\$(wc -l ${meta.id}.\$myName.\$i.singleMatch_blastOut | awk '{print\$1}')"


    covCutoff=`cat "${meta.id}.\$myName.newHaps_coverageCutoff" | awk '{print\$1}'`

    if [ \$match_present == 0  ] && [ \$clusDepth -gt ${params.depthNewHap_min} ] && [ \$clusDepth -gt \$covCutoff ]
    
    then 

    blastn -db ${meta.id}.\$myName.newHaps_multi-seq/\$i.fasta -query ${refDB} -evalue 0.001 -perc_identity 75 -qcov_hsp_perc 85 -num_threads 4 -out ${meta.id}.\$myName.\$i.newHaps_blast_RED_SIM -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid"
    
    theMarker=`head -1 ${meta.id}.\$myName.\$i.newHaps_blast_RED_SIM | sed 's/_Hap_.*/_Hap_/'`
    theSequence=`cat ${meta.id}.\$myName.newHaps_multi-seq/\$i | tail -3 | tr -d '\n' | sed 's/^>Read_[0-9]*//g'`
    theSequence_FDA=`cat ${meta.id}.\$myName.newHaps_multi-seq/\$i | tail -1`
    theSpecimen=`echo ${meta.id}`

    echo \$theSpecimen \$theMarker \$theSequence > allNewHapsOut.file

    fi
    done 
    """
}

