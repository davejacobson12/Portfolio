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

process MAPREFS_BOWTIE2 {
    tag "$meta.id"
    label 'bioinformaticProcessing'
    label 'process_low'
    // shell '/bin/bash', '-euo', 'pipefail'
    // errorStrategy 'ignore'

    input:
    tuple val(meta), path(reads)
    each(marker)
    // each(refDB)

    output:
    tuple val(meta), path("*mapped.bam"), emit: bam, optional: true

    when:
    task.ext.when == null || task.ext.when


    script:
    // def args = task.ext.args ?: ''
    // def args2 = task.ext.args2 ?: ''
    // def prefix = task.ext.prefix ?: '!{meta.id}'
    """


    myName="\$(echo $marker)"
    var="\$(echo \$myName | sed 's/\\n//g')"

    #Use reads[0] reads[1] command when you have unmerged reads from WGS data
    #bowtie2 -x $projectDir/REFERENCES/pv_geo_refs/\$var\\_bt2 -1 ${reads[0]} -2 ${reads[1]} --local  --rg-id ${meta.id} --rg SM:${meta.id} --rg LB:${meta.id} --rg PU:${meta.id} --rg PL:ILLUMINA ${meta.id}_\$var.sam


    bowtie2 -x $projectDir/REFERENCES/pv_geo_refs/\$var\\_bt2 -U ${reads} --local  --rg-id ${meta.id} --rg SM:${meta.id} --rg LB:${meta.id} --rg PU:${meta.id} --rg PL:ILLUMINA ${meta.id}_\$var.sam

    samtools view -b -F 4 ${meta.id}_\$var.sam > ${meta.id}_\$var.mapped.bam
    """

}

    // bowtie2 -x $projectDir/REFERENCES/pv_geo_refs/\$var\\_bt2 -1 ${reads[0]} -2 ${reads[1]} --local  --rg-id ${meta.id} --rg SM:${meta.id} --rg LB:${meta.id} --rg PU:${meta.id} --rg PL:ILLUMINA ${meta.id}_\$var.sam

    // var=`sed 's/\n//g' ${marker}`
    //     echo !{marker} > file1.txt
    // var=`cat !{marker} | tr "\n" ""`
    // echo $var > file2.txt