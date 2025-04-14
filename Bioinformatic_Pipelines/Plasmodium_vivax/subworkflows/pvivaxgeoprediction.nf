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
include { SAMTOOLS_INDEX }              from '../modules/test/samtools/index/main'
include { BOWTIE2_ALIGN }               from '../modules/test/bowtie2/align/main' 

// MODULE: Local Modules Also for Haplotype Calling Workflow
include { BBMAP_COMBINE         }   from '../modules/local/combine_bbmerge'

//MODULE: Local Modules for Geo Prediction
include { MAPREFS_BOWTIE2 }         from '../modules/local/geo_mapbowtie2'
include { MERGE_BAMS }              from '../modules/local/geo_mergebams'
include { GATHER_GVCFS }            from '../modules/local/geo_combinegvcfs'
include { GENOTYPE_GATK }           from '../modules/local/geo_genotypegatk'
include { VARIANTSTOTABLE_GATK }    from '../modules/local/geo_variantstotable'
include { BARCODE_GATK }            from '../modules/local/geo_gatkbarcode'
include { BALK_PREDICT }            from '../modules/local/geo_predictgeo'
include { FORMAT_PREDICTION_DB }    from '../modules/local/geo_formatR'

//Module install from git stopped working for some reason, so I now am using my own installs of tools
include { PICARD_CLEANSAM }             from '../modules/local/picard/cleansam/main'
include { GATK4_HAPLOTYPECALLER }       from '../modules/local/gatk4/haplotypecaller/main'

// SUBWORKFLOW: Installed directly from nf-core/subworkflows
include { SHORTREAD_HOSTREMOVAL }       from '../subworkflows/nf-core/shortread_hostremoval/'

//MODULE: Local Modules for Geo Prediction - old versions of workflow
// include { mapToRefs_bowtie2; cleanSAM_bowtie2; indexBAM_bowtie2; gatkHapCaller_bowtie2; gatherGVCFs; genotypeGATK; variantsToTable_gatk; gatkBarcode; predictGeo_gatk }        from '../modules/local/geo_hapcalling'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow VIVAX_GEO {

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
    BBMAP_BBMERGE(
        BBMAP_BBDUK.out.reads,
        false //interleave
    )

    //  Local process for combining merged and unmerged into single file
    BBMAP_COMBINE(
        BBMAP_BBMERGE.out.merged,
        BBMAP_BBMERGE.out.unmerged
    )

    //  add 'single_end:true' to meta. Necessary for proper file handling in nf-core modules
    ch_combined_reads = BBMAP_COMBINE.out.clean_merged.map { meta, file ->
        [ meta + [ single_end: true], file ]
    }

    //  SHORTREAD_HOSTREMOVAL from TAXPROFILER NF-Core workflow
    //  Map to human genome and keep the unmapped reads
    //  using merged reads
    SHORTREAD_HOSTREMOVAL (
        ch_combined_reads,
        params.hostremoval_reference,
        [[],params.shortread_hostremoval_index] // need to set a tuple for the index
    )

    //  using unmerged reads - will need to change the MAPREFS_BOWTIE2 process to use -1 -2 instead of -U
    //  Probably should include this as a flag but we're unlikely to ever use unmerged reads for ampliseq data.
    //  May be worth including for WGS data
    //     SHORTREAD_HOSTREMOVAL (
    //         BBMAP_BBDUK.out.reads,
    //         params.hostremoval_reference,
    //         [[],params.shortread_hostremoval_index] // need to set a tuple for the index
    //     )

    // New Adapation of workflow to split out each amplicon reference into a different process, rather than a large loop for each sample
    Channel
        .fromPath("$projectDir/REFERENCES/pv_geo_refs/eachAmpliconGeo.txt")
        .splitText()
        .set{ loci_ch }
   
    /* Abbreviated marker list for testing
    Channel
        .fromPath("$projectDir/REFERENCES/pv_geo_refs/eachAmpliconGeo_short.txt")
        .splitText()
        .set{ loci_ch }
    */
   
    //  Own script for bowtie2 mapping to each marker listed in the loci_ch. Wanted more flexibility
    MAPREFS_BOWTIE2 (
        SHORTREAD_HOSTREMOVAL.out.reads,
        loci_ch
    )

    //  Merge all marker bams for sample into a single bam.
    //  This is a locally created module because otherwise I could not get the bams into the same order for different samples
    MERGE_BAMS (
        MAPREFS_BOWTIE2.out.bam.groupTuple()
    )

    //  Nf-core modules to clean and index bam
    PICARD_CLEANSAM (
        MERGE_BAMS.out.bam
    )

    SAMTOOLS_INDEX (
        PICARD_CLEANSAM.out.bam
    )

    //  Combine bam and bai into 1 tuple. Empty lists are required to match cardinality needed for optional GATK4_HAPLOTYPECALLER module
    ch_bam_for_gatk = PICARD_CLEANSAM.out.bam.combine(SAMTOOLS_INDEX.out.bai, by:0).map { meta, file1, file2 ->
        [meta, file1, file2, [], []]
    }

    //  NF-Core module for haplotypecaller. Empty meta tags for reference files
    //  Added -ERC GVCF to modules.conf
    GATK4_HAPLOTYPECALLER (
        ch_bam_for_gatk,
        [[],"$projectDir/REFERENCES/pv_geo_refs/P_VIVAX_ALL_LONG_HAPS_READ_RECOVERY_ALL_CHROMOSOMES_JULY_26_2023.fasta"],
        [[],"$projectDir/REFERENCES/pv_geo_refs/P_VIVAX_ALL_LONG_HAPS_READ_RECOVERY_ALL_CHROMOSOMES_JULY_26_2023.fasta.fai"],
        [[],"$projectDir/REFERENCES/pv_geo_refs/P_VIVAX_ALL_LONG_HAPS_READ_RECOVERY_ALL_CHROMOSOMES_JULY_26_2023.dict"],
        [[],[]],
        [[],[]]
    )

    //  a bit taken from mycosnp nf-core made by Hunter Seabolt - he used a local gatk combingvcfs (https://github.com/CDCgov/mycosnp-nf/blob/master/workflows/mycosnp.nf)
    //  This drops the meta from the output tuple. 
    ch_vcf = GATK4_HAPLOTYPECALLER.out.vcf.map{meta, vcf ->[ vcf ]  }.collect()
    ch_vcf_idx = GATK4_HAPLOTYPECALLER.out.tbi.map{meta, idx ->[ idx ]  }.collect()

    /////////////// From here down, these are all modules that are adaptations of modules used in original workflow. Slightly adapated for this nf-core-like workflow  ///////////////////////////////////

    //  Combine each samples gVCF into
    GATHER_GVCFS(
        ch_vcf.collect(),
        ch_vcf_idx.collect()
    )

    //  Genotype
    GENOTYPE_GATK (
        GATHER_GVCFS.out.combinedGVCFs,
        GATHER_GVCFS.out.combinedIndex
    )

    //  Convert to table
    VARIANTSTOTABLE_GATK (
        GENOTYPE_GATK.out.genotypedVCFs
    )

    //  Format for balk classifier
    BARCODE_GATK(
        VARIANTSTOTABLE_GATK.out.vcfTable, ch_vcf.collect()
    )

    //  Run balk classifier
    BALK_PREDICT (
        BARCODE_GATK.out.barcodeOut
    )

    // format geographic prediction for down stream
    updatedDB = Channel.fromPath("$projectDir/GeoPrediction_output/predictedOut/*geoPrediction_simple.txt").toSortedList().flatten().last()
    FORMAT_PREDICTION_DB (
        updatedDB,
        BALK_PREDICT.out.gatkCountry,
        BALK_PREDICT.out.gatkRegion,
        BARCODE_GATK.out.barcodeOut
    )
}


workflow VIVAX_GEO_existingVCF {

    ch_vcf = Channel.fromPath("$projectDir/GeoPrediction_output/variantCalls/*vcf*")
            

    //  Combine each samples gVCF into
    GATHER_GVCFS(
        ch_vcf.collect()
        // ch_vcf_idx.collect()
    )

    //  Genotype
    GENOTYPE_GATK (
        GATHER_GVCFS.out.combinedGVCFs,
        GATHER_GVCFS.out.combinedIndex
    )

    //  Convert to table
    VARIANTSTOTABLE_GATK (
        GENOTYPE_GATK.out.genotypedVCFs
    )

    //  Format for balk classifier
    BARCODE_GATK(
        VARIANTSTOTABLE_GATK.out.vcfTable, ch_vcf.collect()
    )

    //  Run balk classifier
    BALK_PREDICT (
        BARCODE_GATK.out.barcodeOut
    )

    // format geographic prediction for down stream
    updatedDB = Channel.fromPath("$projectDir/GeoPrediction_output/predictedOut/*geoPrediction_simple.txt").toSortedList().flatten().last()
    FORMAT_PREDICTION_DB (
        updatedDB,
        BALK_PREDICT.out.gatkCountry,
        BALK_PREDICT.out.gatkRegion,
        BARCODE_GATK.out.barcodeOut
    )
}

workflow VIVAX_GEO_SWGA {

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
    // SHORTREAD_HOSTREMOVAL (
    //     ch_combined_reads,
    //     params.hostremoval_reference,
    //     [[],params.shortread_hostremoval_index] // need to set a tuple for the index
    // )

    //  using unmerged reads - will need to change the MAPREFS_BOWTIE2 process to use -1 -2 instead of -U
    //  Probably should include this as a flag but we're unlikely to ever use unmerged reads for ampliseq data.
    //  May be worth including for WGS data
        SHORTREAD_HOSTREMOVAL (
            BBMAP_BBDUK.out.reads,
            params.hostremoval_reference,
            [[],params.shortread_hostremoval_index] // need to set a tuple for the index
        )

    // New Adapation of workflow to split out each amplicon reference into a different process, rather than a large loop for each sample
    Channel
        .fromPath("$projectDir/REFERENCES/pv_geo_refs/eachAmpliconGeo.txt")
        .splitText()
        .set{ loci_ch }
   
    /* Abbreviated marker list for testing
    Channel
        .fromPath("$projectDir/REFERENCES/pv_geo_refs/eachAmpliconGeo_short.txt")
        .splitText()
        .set{ loci_ch }
    */
   
    //  Own script for bowtie2 mapping to each marker listed in the loci_ch. Wanted more flexibility
    MAPREFS_BOWTIE2_WGS (
        SHORTREAD_HOSTREMOVAL.out.reads,
        loci_ch
    )

    //  Merge all marker bams for sample into a single bam.
    //  This is a locally created module because otherwise I could not get the bams into the same order for different samples
    MERGE_BAMS (
        MAPREFS_BOWTIE2.out.bam.groupTuple()
    )

    //  Nf-core modules to clean and index bam
    PICARD_CLEANSAM (
        MERGE_BAMS.out.bam
    )

    SAMTOOLS_INDEX (
        PICARD_CLEANSAM.out.bam
    )

    //  Combine bam and bai into 1 tuple. Empty lists are required to match cardinality needed for optional GATK4_HAPLOTYPECALLER module
    ch_bam_for_gatk = PICARD_CLEANSAM.out.bam.combine(SAMTOOLS_INDEX.out.bai, by:0).map { meta, file1, file2 ->
        [meta, file1, file2, [], []]
    }

    //  NF-Core module for haplotypecaller. Empty meta tags for reference files
    //  Added -ERC GVCF to modules.conf
    GATK4_HAPLOTYPECALLER (
        ch_bam_for_gatk,
        [[],"$projectDir/REFERENCES/pv_geo_refs/P_VIVAX_ALL_LONG_HAPS_READ_RECOVERY_ALL_CHROMOSOMES_JULY_26_2023.fasta"],
        [[],"$projectDir/REFERENCES/pv_geo_refs/P_VIVAX_ALL_LONG_HAPS_READ_RECOVERY_ALL_CHROMOSOMES_JULY_26_2023.fasta.fai"],
        [[],"$projectDir/REFERENCES/pv_geo_refs/P_VIVAX_ALL_LONG_HAPS_READ_RECOVERY_ALL_CHROMOSOMES_JULY_26_2023.dict"],
        [[],[]],
        [[],[]]
    )

    //  a bit taken from mycosnp nf-core made by Hunter Seabolt - he used a local gatk combingvcfs (https://github.com/CDCgov/mycosnp-nf/blob/master/workflows/mycosnp.nf)
    //  This drops the meta from the output tuple. 
    ch_vcf = GATK4_HAPLOTYPECALLER.out.vcf.map{meta, vcf ->[ vcf ]  }.collect()
    ch_vcf_idx = GATK4_HAPLOTYPECALLER.out.tbi.map{meta, idx ->[ idx ]  }.collect()

    /////////////// From here down, these are all modules that are adaptations of modules used in original workflow. Slightly adapated for this nf-core-like workflow  ///////////////////////////////////

    //  Combine each samples gVCF into
    GATHER_GVCFS(
        ch_vcf.collect(),
        ch_vcf_idx.collect()
    )

    //  Genotype
    GENOTYPE_GATK (
        GATHER_GVCFS.out.combinedGVCFs,
        GATHER_GVCFS.out.combinedIndex
    )

    //  Convert to table
    VARIANTSTOTABLE_GATK (
        GENOTYPE_GATK.out.genotypedVCFs
    )

    //  Format for balk classifier
    BARCODE_GATK(
        VARIANTSTOTABLE_GATK.out.vcfTable, ch_vcf.collect()
    )

    //  Run balk classifier
    BALK_PREDICT (
        BARCODE_GATK.out.barcodeOut
    )

    // format geographic prediction for down stream
    updatedDB = Channel.fromPath("$projectDir/GeoPrediction_output/predictedOut/*geoPrediction_simple.txt").toSortedList().flatten().last()
    FORMAT_PREDICTION_DB (
        updatedDB,
        BALK_PREDICT.out.gatkCountry,
        BALK_PREDICT.out.gatkRegion,
        BARCODE_GATK.out.barcodeOut
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



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
