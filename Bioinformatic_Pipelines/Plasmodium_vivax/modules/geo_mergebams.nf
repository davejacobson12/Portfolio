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

process MERGE_BAMS {
    tag "$meta.id"
    label 'bioinformaticProcessing'
    label 'process_high'
    // shell '/bin/bash', '-euo', 'pipefail'
    // errorStrategy 'ignore'

    input:
    tuple val(meta), path(bams)
    // each(marker)
    // each(refDB)

    output:
    tuple val(meta), path("*mapped_sorted_merged.bam"), emit: bam, optional: true

    when:
    task.ext.when == null || task.ext.when


    script:
    // def args = task.ext.args ?: ''
    // def args2 = task.ext.args2 ?: ''
    // def prefix = task.ext.prefix ?: '!{meta.id}'
    
    """

    mkdir bam_folder
    cp ${bams} bam_folder
    cd bam_folder
    ls *bam > myList.txt
    var=`cat myList.txt | tr "\\n" " "`

    samtools merge ${meta.id}_mapped_merged.bam \$var

    samtools sort ${meta.id}_mapped_merged.bam > ../${meta.id}.mapped_sorted_merged.bam



    """
}

