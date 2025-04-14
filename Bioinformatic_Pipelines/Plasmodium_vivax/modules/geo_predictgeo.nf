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

process BALK_PREDICT {
    label 'balkClassifier'
    label 'process_low'
    // shell '/bin/bash', '-euo', 'pipefail'

  input:
  path sampleBarcode

  output:
  path "samples_gatk_region_predicted.out", emit: gatkRegion
  path "samples_gatk_country_predicted.out", emit: gatkCountry

  script:
  """

  ${params.balkClassifier} predict ${sampleBarcode} -o samples_gatk_region_predicted.out --jlibfile  ${params.balkRegion_JLIB}

  ${params.balkClassifier} predict ${sampleBarcode} -o samples_gatk_country_predicted.out --jlibfile  ${params.balkCountry_JLIB}
  """
}



