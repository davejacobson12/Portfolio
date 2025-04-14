

process FORMAT_PREDICTION_DB {
    label 'reporting'
    label 'process_low'
    // shell '/bin/bash', '-euo', 'pipefail'

    publishDir "$projectDir/GeoPrediction_output/predictedOut/", pattern: "*.txt", mode: "copy"

 input:
  path(existingSimpleDB)
  path(gatkCountry)
  path(gatkRegion)
  path(gatkBarcode)

  output:
  path "*prediction_withMetadata.txt", emit:barcodeMeta
  path "*pvivax_geoPrediction_simple.txt", emit:barcodeSimple
    //simple
    //updated formatted predictions

  script:
  """
    Rscript $projectDir/bin/geoprediction_formatOut.R ${existingSimpleDB} ${gatkCountry} ${gatkRegion} ${gatkBarcode}

  """
}
