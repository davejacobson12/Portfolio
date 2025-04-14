#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pvivax/pvivaxclustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/pvivax/pvivaxclustering
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
// params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//workflows to create reports
include { QC_ONLY }                     from './workflows/pvivaxreports'
include { STATE_REPORTS }               from './workflows/pvivaxreports'

//workflows for haplotype calling/clustering
include { PVIVAXCLUSTERING }            from './workflows/pvivaxclustering'
include { PVIVAXCLUSTERING_SWGA }       from './workflows/pvivaxclustering'
include {HAP_SHEET_ONLY }               from './workflows/pvivaxclustering'
include { DISTANCE_MATRIX_CLUSTERING }  from './workflows/pvivaxclustering'
include { CLUSTER_ONLY }                from './workflows/pvivaxclustering'

//workflows for geo prediction
include { VIVAX_GEO }                   from './workflows/pvivaxgeoprediction'
include { VIVAX_GEO_existingVCF }       from './workflows/pvivaxgeoprediction'
include { VIVAX_GEO_SWGA }              from './workflows/pvivaxgeoprediction'


//workflows for drug resistance screening
include { VIVAX_MaRS }                  from './workflows/pvivaxmars'
include { VIVAX_MaRS_WGS }              from './workflows/pvivaxmars'


//workflows for whole genome sequence data
include { VIVAX_SWGA }                  from './workflows/pvivaxswga'
include { VIVAX_SWGA_SNIPPY }           from './workflows/pvivaxswga'
include { VIVAX_SWGA_SNIPPY_EXISTING }  from './workflows/pvivaxswga'





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//workflows to create reports
workflow RUN_MULTIQC {
    QC_ONLY ()
}

workflow PLASMODIUM_REPORTS {
    STATE_REPORTS ()
}

//workflows for haplotype calling/clustering
workflow PVIVAX_PVIVAXCLUSTERING {
    PVIVAXCLUSTERING ()
}

workflow PVIVAX_PVIVAXCLUSTERING_SWGA {
    PVIVAXCLUSTERING_SWGA ()
}


workflow HAP_SHEET {
    HAP_SHEET_ONLY ()
}

workflow Mods2_3 {
    DISTANCE_MATRIX_CLUSTERING ()
}

workflow Mod3 {
    CLUSTER_ONLY ()
}

//workflows for geo prediction
workflow PVIVAX_PVIVAXGEO {
    VIVAX_GEO ()
}

workflow PVIVAX_PVIVAXGEO_SWGA {
    VIVAX_GEO_SWGA ()
}


workflow PVIVAXGEO_EXISTINGVCFS {
    VIVAX_GEO_existingVCF ()
}

//workflows for drug resistance screening
workflow PVIVAX_DRUG_RES {
    VIVAX_MaRS ()
}

workflow PVIVAX_DRUG_RES_WGS {
    VIVAX_MaRS_WGS ()
}

//workflows for whole genome sequence data
workflow PVIVAX_WGS {
    VIVAX_SWGA ()
}

workflow PVIVAX_WGS_SNIPPY {
    VIVAX_SWGA_SNIPPY ()
}

workflow PVIVAX_WGS_SNIPPY_EXISTING {
    VIVAX_SWGA_SNIPPY_EXISTING ()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
