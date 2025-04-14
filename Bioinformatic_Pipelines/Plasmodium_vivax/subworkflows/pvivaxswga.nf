/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowPvivaxclustering.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/test/fastqc/main'
include { MULTIQC                     } from '../modules/test/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/test/custom/dumpsoftwareversions/main'
include { BBMAP_BBDUK }                 from '../modules/test/bbmap/bbduk/main' 
include { BBMAP_BBMERGE }               from '../modules/test/bbmap/bbmerge/main' 
include { BOWTIE2_ALIGN }               from '../modules/test/bowtie2/align/main' 
include { SNIPPY_RUN }               from '../modules/test/snippy/run/main' 
include { SNIPPY_CORE }               from '../modules/test/snippy/core/main' 
include { SNPSITES }               from '../modules/test/snpsites/main' 
include { IQTREE }               from '../modules/test/iqtree/main' 

// MODULE: Local Modules Also for Haplotype Calling Workflow
include { BBMAP_COMBINE         }       from '../modules/local/combine_bbmerge'

// SUBWORKFLOW: Installed directly from test/subworkflows
include { SHORTREAD_HOSTREMOVAL }       from '../subworkflows/nf-core/shortread_hostremoval/'

//Local modules for checking coverage
include { BED_GENOMECOV  }              from '../modules/local/bed_cov'
include { GENOMECOV_SUMMARY  }          from '../modules/local/bed_cov'
include { MAPPING_SUMMARY  }            from '../modules/local/bed_cov'
include { COMBINE_SUMMARY }             from '../modules/local/bed_cov'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow VIVAX_SWGA {

    ch_versions = Channel.empty()
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPvivaxclustering.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowPvivaxclustering.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()


    //  take in files from the quality input quality check
    BBMAP_BBDUK (
        INPUT_CHECK.out.reads,
        params.adapters
    )

    //  include ability to merge reads and comibine into single end file (BBMAP_COMBINE)
    //  We want to merge reads for amplicon data
    // BBMAP_BBMERGE(
    //     BBMAP_BBDUK.out.reads,
    //     false //interleave
    // )

    // //  Local process for combining merged and unmerged into single file
    // BBMAP_COMBINE(
    //     BBMAP_BBMERGE.out.merged,
    //     BBMAP_BBMERGE.out.unmerged
    // )

    // //  add 'single_end:true' to meta. Necessary for proper file handling in nf-core modules
    // ch_combined_reads = BBMAP_COMBINE.out.clean_merged.map { meta, file ->
    //     [ meta + [ single_end: true], file ]
    // }

    //  SHORTREAD_HOSTREMOVAL from TAXPROFILER NF-Core workflow
    //  Map to human genome and keep the unmapped reads
    //  using merged reads
    SHORTREAD_HOSTREMOVAL (
        BBMAP_BBDUK.out.reads,
        params.hostremoval_reference,
        [[],params.shortread_hostremoval_index] // need to set a tuple for the index
    )


   //align to the pvp01 genome
    BOWTIE2_ALIGN (
        SHORTREAD_HOSTREMOVAL.out.reads,
        [[],"$projectDir/REFERENCES/pvp01/"],
        [[],"$projectDir/REFERENCES/pvp01/pvp01_genomic.fna"],
        false, //save unaligned
        true //sort bam

    )

    //processes to check coverage of PvP01 genome in a couple of different ways.
    //Results are printed to the $projectDir/coverage_estimates folder
    BED_GENOMECOV (
        BOWTIE2_ALIGN.out.bam
    )

    GENOMECOV_SUMMARY (
        BED_GENOMECOV.out.defaultBED
    )

    ch_first_merge = BBMAP_BBDUK.out.log.combine(SHORTREAD_HOSTREMOVAL.out.stats, by:0)
    ch_map_stats = ch_first_merge.combine(BOWTIE2_ALIGN.out.log, by:0)

    MAPPING_SUMMARY (
        ch_map_stats
    )

    COMBINE_SUMMARY (
        MAPPING_SUMMARY.out.summary.collect()
    )
}


workflow VIVAX_SWGA_SNIPPY {

    ch_versions = Channel.empty()
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPvivaxclustering.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowPvivaxclustering.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()


    //  take in files from the quality input quality check
    BBMAP_BBDUK (
        INPUT_CHECK.out.reads,
        params.adapters
    )

    //  include ability to merge reads and comibine into single end file (BBMAP_COMBINE)
    //  We want to merge reads for amplicon data
    // BBMAP_BBMERGE(
    //     BBMAP_BBDUK.out.reads,
    //     false //interleave
    // )

    // //  Local process for combining merged and unmerged into single file
    // BBMAP_COMBINE(
    //     BBMAP_BBMERGE.out.merged,
    //     BBMAP_BBMERGE.out.unmerged
    // )

    // //  add 'single_end:true' to meta. Necessary for proper file handling in nf-core modules
    // ch_combined_reads = BBMAP_COMBINE.out.clean_merged.map { meta, file ->
    //     [ meta + [ single_end: true], file ]
    // }

    //  SHORTREAD_HOSTREMOVAL from TAXPROFILER NF-Core workflow
    //  Map to human genome and keep the unmapped reads
    //  using merged reads
    SHORTREAD_HOSTREMOVAL (
        BBMAP_BBDUK.out.reads,
        params.hostremoval_reference,
        [[],params.shortread_hostremoval_index] // need to set a tuple for the index
    )

    SNIPPY_RUN (
        SHORTREAD_HOSTREMOVAL.out.reads,
        "$projectDir/REFERENCES/pvp01/pvp01_genomic.fna"
    )

    //align to the pvp01 genome to get mapping stats
    BOWTIE2_ALIGN (
        SHORTREAD_HOSTREMOVAL.out.reads,
        [[],"$projectDir/REFERENCES/pvp01/"],
        [[],"$projectDir/REFERENCES/pvp01/pvp01_genomic.fna"],
        false, //save unaligned
        true //sort bam

    )

    //processes to check coverage of PvP01 genome in a couple of different ways.
    //Results are printed to the $projectDir/coverage_estimates folder
    BED_GENOMECOV (
        BOWTIE2_ALIGN.out.bam
    )

    GENOMECOV_SUMMARY (
        BED_GENOMECOV.out.defaultBED
    )

    ch_first_merge = BBMAP_BBDUK.out.log.combine(SHORTREAD_HOSTREMOVAL.out.stats, by:0)
    ch_map_stats = ch_first_merge.combine(BOWTIE2_ALIGN.out.log, by:0)

    MAPPING_SUMMARY (
        ch_map_stats
    )

    COMBINE_SUMMARY (
        MAPPING_SUMMARY.out.summary.collect()
    )

    // don't really need to run snippy core at this point, its only really needed when running all samples you want to have in the tree
    /*
    ch_vcf = SNIPPY_RUN.out.vcf.map { meta, vcf -> 
        [vcf]
        }.collect()
    ch_aligned = SNIPPY_RUN.out.aligned_fa.map{ meta, aligned ->
        [aligned]
        }.collect()

    ch_aligned_full =ch_aligned.map {file ->
        meta = [id:"vivax_snippy"]
        [meta,[file]]
    }

    ch_vcf_full =ch_vcf.map {file ->
        meta = [id:"vivax_snippy"]
        [meta,[file]]
    }

    SNIPPY_CORE (
        ch_vcf,
        ch_aligned,
         "$projectDir/REFERENCES/pvp01/pvp01_genomic.fna"

    )
    */
    
}

workflow VIVAX_SWGA_SNIPPY_EXISTING {


    ch_vcf = Channel.fromPath("$projectDir/snippy/samples_vcfs_forTree/*vcf")
    ch_aligned = Channel.fromPath("$projectDir/snippy/samples_aligned_forTree/*aligned.fa")
   
    SNIPPY_CORE (
        ch_vcf.collect(),
        ch_aligned.collect(),
         "$projectDir/REFERENCES/pvp01/pvp01_genomic.fna"

    )

    SNPSITES (
        SNIPPY_CORE.out.full_aln
    )

    ch_phylo_aln = SNPSITES.out.phylo_aln.map {file ->
        meta = [id:"${params.tree_timestamp}"]
        [meta,[file]]
    }

   

    IQTREE (
        ch_phylo_aln, //tuple val(meta), path(alignment)
        [],// path(tree)
        [],// path(tree_te)
        [],// path(lmclust)
        [],// path(mdef)
        [],// path(partitions_equal)
        [],// path(partitions_proportional)
        [],// path(partitions_unlinked)
        [],// path(guide_tree)
        [],// path(sitefreq_in)
        [],// path(constraint_tree)
        [],// path(trees_z)
        [],// path(suptree)
        [],// path(trees_rf)

    )


}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow QC_ONLY {

    ch_versions = Channel.empty()
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPvivaxclustering.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowPvivaxclustering.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
