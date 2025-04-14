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

process FIRST_CLUSTERING_V2 {
    tag "$meta.id"
    label 'mappingSoftware'
    label 'process_medium'

    // publishDir "${renamedReads}", pattern: "RENAMED*fasta", mode: "copy"

    input:
    // path(clippedReads)
    tuple val(meta), path(clippedReads)

    output:
    path("RENAMED*fasta"), emit: sampleFasta ,  optional: true
    // path '*multi-seq', emit: multiSeq_out, optional: true
    // path '${clippedReads.simpleName}.*.multiOut', emit: multiSeq_Ind, optional: true
    // path '*toCount.fq.clstr', emit:clusterDepth, optional: true
    // path '${clippedReads.simpleName}.*.depth', emit:clusterDepthSingle, optional: true
    // tuple val("$meta.id"), path("*.depth"), emit: depthTuple, optional: true
    // tuple val("$meta.id"), path("*.multiOut"), emit: multiTuple, optional: true
    tuple val(meta),path("*.depth"), path("*.multiOut"), emit: combinedTuple, optional: true


    shell:
    '''
    cat !{clippedReads} | awk '{if(NR%4==1) $0=sprintf("@Read_%d",(1+i++)); print;}' > RENAMED.!{meta.id}.clipped_for_hap_calling.fastq
    
    clipLen="\$(wc -l RENAMED.!{meta.id}.clipped_for_hap_calling.fastq | awk '{print$1}')"
    if [ $clipLen -gt 0 ]
    then
    seqret -sequence RENAMED.!{meta.id}.clipped_for_hap_calling.fastq -outseq RENAMED.!{meta.id}.clipped_for_hap_calling.fasta

    cd-hit-est -i RENAMED.!{meta.id}.clipped_for_hap_calling.fastq -o !{meta.id}.clean_merged_Clusters.fq -c 1 -g 1 -d 0 -T !{params.max_cpus} -s 1 -M 10000

    seqret -sequence !{meta.id}.clean_merged_Clusters.fq -outseq !{meta.id}.clean_merged_Clusters.fasta

    perl !{params.multiSeq_pl} !{meta.id}.clean_merged_Clusters.fasta !{meta.id}.clean_merged_Clusters.fq.clstr !{meta.id}.multi-seq !{params.assignHap_min}
    
    else
    touch RENAMED.!{meta.id}.clipped_for_hap_calling.fasta
    mkdir !{meta.id}.multi-seq
    touch !{meta.id}.clean_merged_Clusters.fq.clstr 
    fi

    ls !{meta.id}.multi-seq > !{meta.id}.list

    cat !{meta.id}.clean_merged_Clusters.fq.clstr | sed 's/Cluster /Cluster/g' > !{meta.id}.toCount.fq.clstr

    clusLen="\$(wc -l !{meta.id}.toCount.fq.clstr | awk '{print$1}')"

    if [ $clusLen -gt 0 ]
    then
    python3 /scicomp/groups-pure/Projects/CycloDeploy/cyclone/cyclone_nf_core/cyclo-eightmarker/bin/cluster_counter_v2.py -x !{meta.id}.toCount.fq.clstr -y !{meta.id}.clusterCounted.csv

    cat !{meta.id}.list | while read j
    do
    cp !{meta.id}.multi-seq/$j !{meta.id}.Cluster$j.trial
    cat !{meta.id}.Cluster$j.trial | sed "s/>Read/>Cluster${j}-Read/g" > !{meta.id}.Cluster$j.multiOut
    grep "Cluster$j," !{meta.id}.clusterCounted.csv > !{meta.id}.Cluster$j.depth
    done
    fi

    '''
}

