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

process READ_RENAME{
    tag "$meta.id"
    label 'mappingSoftware'
    label 'process_medium'
    // shell '/bin/bash', '-euo', 'pipefail'

    input:
    tuple val(meta), path(samMapped)

    output:
    tuple val(meta), path("*mapped_only.fastq")      , emit: renamed     , optional:false

    when:
    task.ext.when == null || task.ext.when

    shell:
    // def args = task.ext.args ?: ''
    // def args2 = task.ext.args2 ?: ''
    // def prefix = task.ext.prefix ?: '!{meta.id}'

    '''
    samtools view !{samMapped} | awk '{print("@"\$1"\\n"\$10"\\n+\\n"\$11)}' | awk '{if(NR%4==1) \$0=sprintf("@Read_%d",(1+i++)); print;}' > !{meta.id}.mapped_only.fastq
    
    '''
}

process EMPTY_GENOTYPE {
    tag "$meta.id"
    label 'process_low'
    // label 'mappingSoftware'
    // shell '/bin/bash', '-euo', 'pipefail'

    input:
    tuple val(meta), path(reads)

    output:
    // tuple val(meta), path("*mapped_only.fastq")      , emit: renamed     , optional:false

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    // def args2 = task.ext.args2 ?: ''
    // def prefix = task.ext.prefix ?: '!{meta.id}'

    """
    sampleName=`echo ${meta.id} | sed 's/_T1//g' | sed 's/-/_/g'`

    echo \$sampleName > testFile.txt
    touch $projectDir/HAPLOTYPE_CALLER_CYCLO_V2/SPECIMEN_GENOTYPES/\$sampleName
    """
}

