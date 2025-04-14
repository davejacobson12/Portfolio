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

WorkflowEightmarker.initialise(params, log)

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



// MODULE: Local Modules for Haplotype Calling
include { FINALNEWREFSV2                }   from '../modules/local/finalnewrefsv2'
include { BBMAP_COMBINE                 }   from '../modules/local/combine_bbmerge'
include { READ_RENAME; EMPTY_GENOTYPE   }   from '../modules/local/rename_bam2fastq'
include { BEDEXT_NEWHAPS                }   from '../modules/local/bedExtraction_plusNewHaps' 
include { MAKE_FINAL_READS_V2           }   from '../modules/local/makeFinalReads_v2'
include { NAME_NEW_HAPS                 }   from '../modules/local/nameNewHaps'
// include { MAKE_REFS_ASSIGN_HAPS_V2      }   from '../modules/local/makeRefs_assignHaps2_v2'
include { FINAL_NEW_REFS_ASSIGN_HAPS_V2 }   from '../modules/local/finalNewRefs_assignHaps2' 
// include { BLAST_REF_FINAL }                 from '../modules/local/blastDB'
include { FIRST_CLUSTERING_V2 }             from '../modules/local/firstClustering_v2'
include { BLAST_MULTI_V2 }                  from '../modules/local/blastMulti_v2'
include { HAPS_BLAST }                      from '../modules/local/hapsBlast'
include { HAPLOTYPE_SHEET; combineGenotypes}                  from '../modules/local/haplotype_sheet'

include { runModule2 } from '../modules/local/workflow2'
include { runModule3_v2 } from '../modules/local/workflow3'
include { MERGE_MASTERS; TGC_DETECTION; HTML_GENERATION } from '../modules/local/report_generation'

//for the junction marker

include { makeJunctionRefs; finalJunctionNewRefs; mappingJunction; getAssembleSeqs; clusterSeqs_1; miraTrials; firstPass_Junction; firstPass_JunctionCheck1; clusterFirstPass; secondPass_junction; thirdPass_junction; newJunction_finder; nameNewJuncHaps; makeRefs_assignHaps_junc2; finalNewRefs_assignHaps_junc2; juncBlast } from '../modules/local/junctionMarker'


//  MODULE: Local modules for distance matrix and clustering

// include { mod2_importData } from '../modules/local/distance_matrix_workflow'
// include { mod2_mammelCode_v2 } from '../modules/local/distance_matrix_workflow'

// include { runModule3_Pvivax } from '../modules/local/clustering_workflow'
// include { runModule3_parnas } from '../modules/local/clustering_workflow'

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
include { FASTQC                      } from '../modules/test/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/test/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/test/nf-core/custom/dumpsoftwareversions/main'

include { BBMAP_BBDUK }                 from '/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/pvivax_nf_core/pvivax-pvivaxclustering/modules/test/bbmap/bbduk/main' 
include { BBMAP_BBMERGE }               from '/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/pvivax_nf_core/pvivax-pvivaxclustering/modules/test/bbmap/bbmerge/main' 
include { BWA_INDEX }                   from '/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/pvivax_nf_core/pvivax-pvivaxclustering/modules/test/bwa/index/main'  
include { BWA_MEM }                     from '/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/pvivax_nf_core/pvivax-pvivaxclustering/modules/test/bwa/mem/main' 
include { BOWTIE2_ALIGN }               from '/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/pvivax_nf_core/pvivax-pvivaxclustering/modules/test/bowtie2/align/main' 
include { SAMTOOLS_INDEX }              from '/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/pvivax_nf_core/pvivax-pvivaxclustering/modules/test/samtools/index/main'

include { BOWTIE2_BUILD }               from '/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/pvivax_nf_core/pvivax-pvivaxclustering/modules/local/bowtie2/build/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

tmpDir = file("$projectDir/HAPLOTYPE_CALLER_CYCLO_V2/TMP2")
tmpDir.deleteDir()
tmpDir.mkdir()

// Info required for completion email and summary
def multiqc_report = []

workflow EIGHTMARKER {

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
    workflow_summary    = WorkflowEightmarker.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowEightmarker.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
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

        // Create an empty genotype for each sample analyzed
    EMPTY_GENOTYPE {
        INPUT_CHECK.out.reads
    }

    // Take in reads from input check
    BBMAP_BBDUK (
        INPUT_CHECK.out.reads,
        params.adapters
    )



    // merge and combine amplicon reads
    BBMAP_BBMERGE(
        BBMAP_BBDUK.out.reads,
        false //interleave
    )

    BBMAP_COMBINE(
        BBMAP_BBMERGE.out.merged,
        BBMAP_BBMERGE.out.unmerged
    )
    /*
        concat all reference haplotypes into two files
            First file has the original full length references plus any new haps
            Second file has just the short marker part references plus any new haps
    */

    FINALNEWREFSV2 (
        MULTIQC.out.report
    )

    //  Convert those files into a channel with a map (map is necessary for nf-core bwa and bowtie2 indexing)
    ch_full_refs = FINALNEWREFSV2.out.outFinalRefs.map {file ->
        meta = [id:"fullRefs"]
        [meta,[file]]
    }

    ch_short_refs = FINALNEWREFSV2.out.outShortRefs.map {file ->
        meta = [id:"shortRefs"]
        [meta,[file]]
    }

    BWA_INDEX (
        ch_full_refs
    )

    // using the unmerged reads
    // BWA_MEM (
    //     BBMAP_BBDUK.out.reads,
    //     BWA_INDEX.out.index,
    //     ch_full_refs,
    //     false
    // )

    // use merged reads to do initial mapping to ampliseq refs
    BWA_MEM (
        BBMAP_COMBINE.out.clean_merged,
        BWA_INDEX.out.index,
        ch_full_refs,
        false
    )

    //equivalent of mapping2fastq - rename to have reads named Read1, Read2, etc
    READ_RENAME (
        BWA_MEM.out.bam
    )

    //set the single_end meta to true for bowtie alignment
    ch_rename_reads = READ_RENAME.out.renamed.map { meta, file ->
        [ meta + [ single_end: true], file ]
    }

    //now set up channel to bowtie2 align to full length markers
    ch_read_rec = [
        [id:'readRecoveryRefs'],[params.fullNonJunction_refs]
    ]

    // build bowtie database of refs and previously bwa aligned reads to the bowtie2 ref
    BOWTIE2_BUILD (
        ch_read_rec
    )

    BOWTIE2_ALIGN (
        ch_rename_reads,
        BOWTIE2_BUILD.out.index,
        ch_read_rec,
        false, //save unaligned
        true //sort_bam
    )

    // index the alignment for each sample and combine into single tuple with the alignment
    SAMTOOLS_INDEX(
        BOWTIE2_ALIGN.out.bam
    )

    ch_mapped_reads = BOWTIE2_ALIGN.out.bam.combine(SAMTOOLS_INDEX.out.bai, by:0)

    // create a channel with each marker and create jobs for each marker. Make sure the full marker list is used
    Channel
        .fromPath("${params.loci}")
        .splitText()
        .set{ loci_ch }

    // adapt previously existing workflow to run with nf-core. BEDEXT_NEWHAPS is a complex module that has many parts, not split into different modules for ease of use
    BEDEXT_NEWHAPS (
        ch_mapped_reads, 
        loci_ch, 
        FINALNEWREFSV2.out.outShortRefs
    )

    //  get final set of reads for haplotype calling for each sample
    MAKE_FINAL_READS_V2(
        BEDEXT_NEWHAPS.out.tryBedFastq.groupTuple()
    )
  
    // assign names to any new haplotypes 
    NAME_NEW_HAPS (
        BEDEXT_NEWHAPS.out.allNewHaps.collectFile(), 
        FINALNEWREFSV2.out.outShortRefs
    )

    // make updated ref db, including newly identified haplotypes 
    FINAL_NEW_REFS_ASSIGN_HAPS_V2(
        // makeRefs_v2.out.completeFile.collectFile(), //don't have this yet, need to check on it
        NAME_NEW_HAPS.out.newHaps_singleFile.collectFile().ifEmpty("$projectDir/assets/empty.txt")
    )

    //  cluster identical reads together for each sample and find the depth of each of those reads
    FIRST_CLUSTERING_V2 (
        MAKE_FINAL_READS_V2.out.clipped
    )

    //  now combine each clustered read and depth with the cutoff required to call a haplotype at each marker into a single tuple
    ch_hap_cov = BEDEXT_NEWHAPS.out.markerCov.groupTuple()
    ch_cluster_hap_cov = FIRST_CLUSTERING_V2.out.combinedTuple.combine(ch_hap_cov, by: 0)

    // BLAST for each clustered read, using hap coverage to filter out markers with low coverage
    BLAST_MULTI_V2(
        ch_cluster_hap_cov,
        FINAL_NEW_REFS_ASSIGN_HAPS_V2.out.outFinalRefs
    )

    // final blast to assign haplotypes for each sample
    HAPS_BLAST (
        BLAST_MULTI_V2.out.blastFasta, 
        FINAL_NEW_REFS_ASSIGN_HAPS_V2.out.outFinalRefs
    )

    newJunctions = Channel.fromPath("${params.origNewHaps}/*Junction*.fasta")

    //makeJunctionRefs(newJunctions.ifEmpty(""))
    finalJunctionNewRefs(MULTIQC.out.report)

    //This is if we use santizeme for read cleaning
    // mappingJunction(finalJunctionNewRefs.out.outJunctionRefs, SanitizeMe2.out.cleaned)

    //This is if we use a read count and quality score cutoff to decide what goes into the workflow
    // mappingJunction(finalJunctionNewRefs.out.outJunctionRefs, fastqc_afterMerge.out.qcFASTQs )

    //This is if we don't use any initial filtering step
    // ch_junctionReads = BBMAP_COMBINE.out.clean_merged.map{meta, file ->
    //     [file]
    // }
    mappingJunction(finalJunctionNewRefs.out.outJunctionRefs, BBMAP_COMBINE.out.clean_merged)
    // getAssembleSeqs(mappingJunction.out.mappedSingleLine)
    getAssembleSeqs(mappingJunction.out.mappedSingleLine_tuple)

    // cluster the sequences (i don't think we need the bbduk ready reads here but i will test)
    // clusterSeqs_1(getAssembleSeqs.out.forClustering)
    clusterSeqs_1(getAssembleSeqs.out.forClustering_tuple)

    ch_for_mira = BBMAP_COMBINE.out.clean_merged.combine(clusterSeqs_1.out.forManifest_tuple, by:0)
    // ch_for_mira.view()

    // use mira to assemble the reads and the remaining passes are repeated mapping/trimming steps that follow logical routines to remove artifact sequences
    // miraTrials(clusterSeqs_1.out.forManifest)
    miraTrials(ch_for_mira)

    firstPass_Junction(miraTrials.out.draft0_out)
    firstPass_JunctionCheck1(firstPass_Junction.out.draft1_00002)
    clusterFirstPass(firstPass_JunctionCheck1.out.firstPass_out)
    secondPass_junction(clusterFirstPass.out.draft2_consensus)
    thirdPass_junction(secondPass_junction.out.draft200_out)

    // map to existing junction haplotypes to see if a new junction haplotype was found
    newJunction_finder(thirdPass_junction.out.validatedJunctions)

    // collect all new junction haplotypes and use python script to ensure that no duplicate haplotypes are named
    nameNewJuncHaps(newJunction_finder.out.allNewJuncHaps.collectFile())
    // makeRefs_assignHaps_junc2(newJunctions.ifEmpty(""))
    // finalNewRefs_assignHaps_junc2(nameNewJuncHaps.out.newJuncHaps.collectFile().ifEmpty("$projectDir/assets/empty.txt"), makeRefs_assignHaps_junc2.out.completeFile_junc.collectFile())
    finalNewRefs_assignHaps_junc2(nameNewJuncHaps.out.newJuncHaps.collectFile().ifEmpty("$projectDir/assets/empty.txt"))
    juncBlast(thirdPass_junction.out.validatedJunctions, finalNewRefs_assignHaps_junc2.out.outFinalRefs)

    ch_sample = BBMAP_COMBINE.out.clean_merged.map{meta, fastq ->
        fastq
    }
    // ch_sample.view()
    // because this workflow runs the non junction and junction processes in parallel, we need to combine the junction and non junction genotypes for each specimen
    // combineGenotypes(HAPS_BLAST.out.genotypeNoJunction.ifEmpty("$projectDir/assets/emptySNP_genotype").collectFile(), juncBlast.out.genotypeJunction.ifEmpty("$projectDir/assets/emptyJunc_genotype").collectFile())
    //   combineGenotypes(ch_sample, HAPS_BLAST.out.genotypeNoJunction.ifEmpty("$projectDir/assets/emptySNP_genotype"), juncBlast.out.genotypeJunction.ifEmpty("$projectDir/assets/emptySNP_genotype"))
      combineGenotypes(ch_sample, HAPS_BLAST.out.genotypeNoJunction.collect(), juncBlast.out.genotypeJunction.collect())

    // place the genotypes in the correct folder with the correct name
    // copyGenotypes(combineGenotypes.out.finalGenotypes.collectFile())

    //
    HAPLOTYPE_SHEET(
        combineGenotypes.out.finalGenotypes.collect()
    ) 
}


workflow HAP_SHEET_ONLY {
    HAPLOTYPE_SHEET (
        file(params.input)
    )
}


workflow Workflow2 {
  hapSheets = Channel.fromPath("$projectDir/haplotype_sheets/*.txt").toSortedList().flatten().last()
  paramInput_wf2 = Channel.fromPath(params.wf2HapSheet)

  if( params.wf2HapSheet == "$projectDir/assets/empty.txt" ){
    runModule2(hapSheets)
  }
  else{
    runModule2(paramInput_wf2)
  }

}

workflow Workflow3_v2 {
  myMat = Channel.fromPath("$projectDir/ensemble_matrices/*H_matrix.csv").toSortedList().flatten().last()
  paramInput_wf3 = Channel.fromPath(params.wf3Matrix)

  if( params.wf3Matrix == "$projectDir/assets/empty.txt"){
    runModule3_v2(myMat)
  }
  else{
    runModule3_v2(paramInput_wf3)
  }
}

workflow Workflow_2and3_v2 {
  hapSheets = Channel.fromPath("$projectDir/haplotype_sheets/*.txt").toSortedList().flatten().last()
  paramInput_wf2 = Channel.fromPath(params.wf2HapSheet)

  if( params.wf2HapSheet == "$projectDir/assets/empty.txt" ){
    runModule2(hapSheets)
    runModule3_v2(runModule2.out.matrixOut)
  }
  else{
    runModule2(paramInput_wf2)
    runModule3_v2(runModule2.out.matrixOut)
  }
}

workflow MODULE4 {
    MERGE_MASTERS ()
}
workflow MODULE5 {
    TGC_DETECTION ()
}

workflow MODULE6 {
    HTML_GENERATION ()
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
