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

process HAPLOTYPE_SHEET {
    // tag "$meta.id"
    label 'pythonEnv_wf1'
    label 'process_medium'

    publishDir "$projectDir/haplotype_sheets", pattern: "*haplotype_shee*.txt", mode: "copy"

    input:
    path(haps)

    output:
    path "*haplotype_sheet.txt", emit: finalHapSheet

    script:
    """
    python3 $projectDir/bin/haplotypeSheet_gen_forNextflow.py  -g ${params.hapSheet_specificGenotypes} -r ${params.hapSheet_specificReferences} -t $projectDir/HAPLOTYPE_CALLER_CYCLO_V2/tempDir -o haplotype_sheet.txt -f $projectDir/HAPLOTYPE_CALLER_CYCLO_V2/FAILED_SPECIMENS -H $projectDir/haplotype_sheets -Q $projectDir/HAPLOTYPE_CALLER_CYCLO_V2/QClogs 
    """


}

process combineGenotypes {
  label 'process_low'
//   publishDir "${tmpGenotypes}", pattern: "*GENOTYPE_final", mode: "copy"
  
  input:
//   tuple val(meta), path(fastq)
  each (inFASTQ)
  path (genotypes1)
  path (genotypes2)
//   tuple val(meta), path(genotypes1), path(genotypes2)

  output:
  path "*GENOTYPE_final", emit: finalGenotypes
  path "*GENOTYPE_final"
  path "endMarker", emit: complete


  script:
  """

  cat ${params.tmpGenotypes}/${inFASTQ.simpleName}.GENOTYPE_noJunction ${params.tmpGenotypes}/${inFASTQ.simpleName}.GENOTYPE_withJunction >> ${inFASTQ.simpleName}.GENOTYPE_final

  sampleName=`echo ${inFASTQ.simpleName} | sed 's/_T1//g' | sed 's/-/_/g'`
  echo \$sampleName > testFile.txt
  cp ${inFASTQ.simpleName}.GENOTYPE_final $projectDir/HAPLOTYPE_CALLER_CYCLO_V2/SPECIMEN_GENOTYPES/\$sampleName

  echo "done" >> endMarker
  """
}

// process copyGenotypes{
//   label 'process_low'
//   input:
//   path endGenotype

//   output:
//   path 'endMarker', emit: endFile

//   script:
//   """
//   cp "${tmpGenotypes}/${endGenotype.simpleName}.GENOTYPE_final" "${finalGenotypes}/${endGenotype.simpleName}"
//   echo "done" >> endMarker
//   """
// }