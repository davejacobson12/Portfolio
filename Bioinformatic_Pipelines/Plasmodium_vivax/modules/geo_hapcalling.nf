process mapToRefs_bowtie2 {
  // errorStrategy 'ignore'
  tag "$meta.id"
  label 'bioinformaticProcessing'

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path("*mapped_sorted_merged.bam"), emit: bowtieMapped

  shell:
  '''
  cat !{params.REFERENCES}/eachAmpliconGeo.txt | while read line; do bowtie2 -x !{params.REFERENCES}/$line\\_bt2 -1 !{reads[0]}  -2 !{reads[1]} --local  --rg-id !{meta.id} --rg SM:!{meta.id} --rg LB:!{meta.id} --rg PU:!{meta.id} --rg PL:ILLUMINA !{meta.id}_$line.sam ;done
  
  cat !{params.REFERENCES}/eachAmpliconGeo.txt | while read line; do samtools view -b -F 4 !{meta.id}_$line.sam > !{meta.id}_$line\\_mapped.bam ; done

  
  samtools merge !{meta.id}_mapped_merged.bam !{meta.id}_*mapped.bam 

  samtools sort !{meta.id}_mapped_merged.bam > !{meta.id}.mapped_sorted_merged.bam

  rm -f -- *_mapped_merged.bam
  rm -f -- *mapped.bam
  rm -f -- *.sam

  '''
}

process cleanSAM_bowtie2 {
  tag "$meta.id"
  label 'bioinformaticProcessing'

  input:
  tuple val(meta), path(mapBowtieSAM)
//   path mapBowtieSAM

  output:
  tuple val(meta), path("*mapped_sorted_merged_clean.bam"), emit: cleanedSAM

  script:
  """
  gatk  CleanSam -R ${params.ampliseqRef} -I ${meta.id}.mapped_sorted_merged.bam -O ${meta.id}.mapped_sorted_merged_clean.bam
  """

}

process indexBAM_bowtie2 {
    tag "$meta.id"
  // label 'variantCalling'
  label 'bioinformaticProcessing'
  
  input:
  tuple val(meta), path(markedBAM)
  output:
  tuple val(meta), path("*.mapped*bai"), emit: indexOut
  script:
  """
  samtools index ${markedBAM}
  """
  
}
//   cp ${markedBAM} ${markedBAM.simpleName}.mapped_sorted_merged_clean_index.bam




process gatkHapCaller_bowtie2{
    tag "$meta.id"
  label 'bioinformaticProcessing'
  
//   publishDir "${params.variantFolder}", pattern: "*gatk_bowtie2.g.vcf", mode: "copy"

  input:
  tuple val(meta), path(bam), path(index)

  output:
  tuple val(meta), path("*gatk_bowtie2.g.vcf"), emit: gatkOut
//   path "*gatk_bowtie2.g.vcf", emit: gatkOut
  script:
  """
  gatk  HaplotypeCaller -R ${params.ampliseqRef} -I ${bam} --native-pair-hmm-threads 8  -ERC GVCF  -O ${meta.id}.ampliseq_gatk_bowtie2.g.vcf 
  """

}

process gatherGVCFs{
 label 'bioinformaticProcessing'

  input:
  path vcfList

  output:
  tuple val("gatk"), path("gatk.g.vc*"), emit: combinedGVCFs

  script:
  """
  echo "${vcfList.join('\n')}" > vcf.list

  gatk CombineGVCFs -R ${params.ampliseqRef} --variant vcf.list -O gatk.g.vcf
  """
}

process genotypeGATK {
 label 'bioinformaticProcessing'

  input:
  tuple val("gatk"), path(gvcf)

  output:
  path "gatk_allLoci.vcf", emit: gatkVCF

  script:
  """
  gatk GenotypeGVCFs -R ${params.ampliseqRef} -V gatk.g.vcf --include-non-variant-sites -O gatk_allLoci.vcf

  """

}

process variantsToTable_gatk {
  label 'bioinformaticProcessing'
  // publishDir "${params.variantFolder}", pattern: "*gatk_variantTable.table", mode: "copy"

  input:
  path gatkVCF

  output:
  path "gatk_variantTable.table", emit: gatkTable

  script:
  """
  gatk VariantsToTable -V ${gatkVCF} -F CHROM -F POS -GF GT -GF DP --show-filtered -O gatk_variantTable.table
  """

}

//Fix this to have the ref files better set up
process gatkBarcode {
  label 'bioinformaticProcessing'

  publishDir "${params.variantFolder}", pattern: "*gatk_pvivax_ampliseq_barcode.csv", mode: "copy"
  publishDir "${params.variantFolder}", pattern: "*.table", mode: "copy"

  input:
  path gatkTable
  path vcfList

  output:
  path "*gatk_pvivax_ampliseq_barcode.csv", emit:barcodeOut
  path "*.table"

  script:
  """
  echo "${vcfList.join('\n')}" > vcf.list
  sed 's/.ampliseq_gatk_bowtie2.g.vcf//g' vcf.list > sampleList.txt

  python3 $projectDir/bin/snpTable_toBarcode.py --balkKey $projectDir/REFERENCES/pv_geo_refs/balkKey2.txt --snpFile ${gatkTable} --refBarcode $projectDir/REFERENCES/pv_geo_refs/geoCombined_barcode_allSamples.csv --sampleList sampleList.txt --outfile gatk_pvivax_ampliseq_barcode.csv 

  """
}

//Could also include steps for creating the classifier as a separte workflow
//Probably also add a timestamp or something uniq to prediction name - maybe same thing with barcodes above
process predictGeo_gatk {
  label 'balkClassifier'

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

