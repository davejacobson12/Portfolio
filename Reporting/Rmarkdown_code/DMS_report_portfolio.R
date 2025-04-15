library(knitr)
library(kableExtra)
library(tibble)
library(dplyr)
library(stringr)
# library(readxl)


############Arguments for Pfal and Pvivax
args = commandArgs(trailingOnly = T)


##### Info for making Rmarkdown html file
document_number = "Not Cleared"
version_number = " - Draft Document"

doc_control <- paste(document_number, version_number, sep="")

time_now <- format(Sys.time(), "%H%M")
date_now <- Sys.Date()
printDate <- paste("Date of Report Generation (yyyy/mm/dd): ", date_now)
fullPath <- "/Users/davejacobson/Desktop/DataScience_Trainings/Portfolio/Reporting/Rmarkdown_code/"

####################### universal preprocessing for pvivax reports ##################################

Abbreviations <- c("AF", "EAS", "ESEA", "LAM", "MSEA", "OCE", "WAS", "WSEA")
RegionNames <- c("Africa", 
"East Asia", 
"East Southeast Asia",
"Central / South America", 
"Malaysia Region Southeast Asia", 
"Oceania", 
"Western Asia",
"Western Southeast Asia")
BranchColor <- rep("",8)
geographicRegions <- data.frame(Abbreviations, RegionNames)
geographicRegionsLegend <- data.frame(Abbreviations, RegionNames, BranchColor)
myColVec <- c("#FF3030", "#0000FF", "#008B00", "#A020F0", "#FFA500", "#FFFF00", "#FFC0CB", "#A9A9A9")


#pvivax geo prediction regions
geoKable <- knitr::kable(geographicRegions, "html", padding = 40, line_sep = 2, align = "cc") %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center", font_size = 12)

geoKableLegend <- knitr::kable(geographicRegionsLegend, "html", padding = 40, line_sep = 2, align = "ccc", caption = "Legend: Branch Colors on Tree Correspond to Predicted Geographic Region") %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center", font_size = 12) %>%
  column_spec(c(3),background = myColVec)

#Some data cleaning
# colnames(parnasClusters) <- c("LSDB_Sequence_ID", "Cluster")
# parnasClusters$Cluster <- parnasClusters$Cluster +1
# fullMeta$Travel <- gsub(" ", "", fullMeta$Travel)
# colnames(fullMeta)[1] <- "Seq_ID"
# fullMeta$Seq_ID <- gsub("-","_", fullMeta$Seq_ID)
# fullMeta$Date_Collected <- as.Date(as.numeric(fullMeta$Date_Collected), origin = "1899-12-30")

#Extract the columns of interest from the metadata sheet
# fullMeta <- fullMeta[,c(1,2,3,4,7,12,13)]

# merge1 <- merge(fullMeta, parnasClusters, by = "Seq_ID")
# merge2 <- merge(merge1, geographic, by = "Seq_ID")

#more data cleaning so that the rows are in order and grouped together by cluster. Remove any duplicate rows (There are some floating around in different sheets)
# merge2 <- merge2[order(merge2$Cluster),]
# merge2 <- merge2[!duplicated(merge2),]

###Finished with p vivax data cleaning ##############

#Make a table that list species detetected by PET / analyses completed in Plasmodium fal/vivax for DMS
petTable <- as.data.frame(matrix(ncol = 3, nrow = 4))

colnames(petTable) <- c("Plasmodium Species", "PET-PCR Detection", "Additional Analyses")
petTable[1,1] <- "P. falciparum"
petTable[2,1] <- "P. vivax"
# petTable[3,1] <- "P. knowlesi"
petTable[3,1] <- "P. ovale"
petTable[4,1] <- "P. malariae"
petTable[c(3,4),3] <- "None"
petTable[c(1,2,3,4),2] <- "Yes"
petTable[1,3] <- "Drug resistance screening, geographic prediction"
petTable[2,3] <- "Genetic clustering, geographic prediction"

#kable for malaria analyses
petTableKable <- knitr::kable(petTable, align = "c",caption = "Table 1: Malaria Genotyping Analyses") %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12) %>%
      column_spec(1, italic = T)


#Make a table that lists the continents analyzed for P fal geopredictions
AbbreviationsPfal <- c("AF","AF-W", "AF-R", "AS-SE", "AS-S", "OC-NG", "SA")
RegionNamesPfal <- c("Africa", "West Africa", "Non-west Africa", "Southeast Asia", "South Asia", "Oceania", "Central / South America")
AnalysisType <- c("Two-Gene Panel", "211-SNP Panel", "211-SNP Panel", "211-SNP Panel", "211-SNP Panel", "Both Approaches", "Both Approaches")
geographicRegionsPfal <- data.frame(AbbreviationsPfal, RegionNamesPfal, AnalysisType)
colnames(geographicRegionsPfal) <- c("Abbreviations", "Geographic Region", "Analysis Approach")

#kable for pfal continents
geoKablePfal <- knitr::kable(geographicRegionsPfal, "html", padding = 40, line_sep = 2, align = "cc") %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center", font_size = 12)

#This stateList is just for testing, eventually delete this to include all states
stateList <- c('NY','TX','CA','FL')
metaIn <- read.csv(paste(fullPath, "input_data/malaria_initial_meta.csv",sep=""), header = T)

#Loop through each state and make reports
for(stateID in stateList){
  # setwd("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/Plasmodium_reporting_folder/pythonOutFiles_tmp")
  print(stateID)
  
  #read in the results file from Dhruvi's python code (this has all metadata and pfal results)
  # setwd("pythonOutFiles_tmp")
  # resultFile <- list.files(pattern =paste(stateID,"_results.csv", sep =""))
  # resultsIn = read.csv(resultFile, header = T, na.strings =c(""))
  resultsIn <- metaIn[grepl(stateID, metaIn$Submitter_State),]
  resultsIn <- resultsIn[order(resultsIn$LSDB_Specimen_ID),]
  print(resultsIn)
  
  #put not tested for hrp Status for now
  #resultsIn['hrpStatus'] <- 'Not_Tested'
  
  rownames(resultsIn) <- NULL
  print("reading in by state")
  resultsIn$Collection_date[resultsIn$Collection_date == "not provided "] <- NA
  resultsIn$Collection_date[resultsIn$Collection_date == "not provided"] <- NA
  resultsIn$Collection_date[resultsIn$Collection_date == "unknown"] <- NA

  #format the columns that have dates. This code looks complicated, but just makes sure that the dates that are read in are converted into mm/dd/yyyy format
  dateColumns <- c(7,16)
  for(i in dateColumns){
    resultsIn[, i] <- sapply(resultsIn[, i], as.character)
  }
    for(i in 1:length(rownames(resultsIn))){
      for(j in dateColumns){
        if(is.na(resultsIn[i,j])){
          resultsIn[i,j] <- resultsIn[i,j]
        }else if(nchar(resultsIn[i,j])==5){
          resultsIn[i,j] <- format(as.Date(as.numeric(resultsIn[i,j]),origin = "1899-12-30"), "%m/%d/%Y")
        }else if(substr(resultsIn[i,j],1,2) == "20"){
          test1 <- as.character(resultsIn[i,j])
          test2 <- substr(test1,1,10)
          test3 <- format(as.Date(test2,"%Y-%m-%d"),"%m/%d/%Y")
          resultsIn[i,j] <- test3
        }
      }
    }

  # resultsForTable1 <- resultsIn[,c(9,2,6,5,4,8,1,12,7,10)]

  #Select metadata from the results file and put into order that data flows
    #LSDB LAB ID | CSID | State Lab ID | Collection Date | Travel History | Specimen Status | ReportStatusDate | Species Detected | Geo Prediction
  names(resultsIn)[names(resultsIn) == 'Specimen_Status'] <- 'Sample_Status'
  names(resultsIn)[names(resultsIn) == 'Patient_Travel_History_Interpretation'] <- 'Reported_Travel_History'
  names(resultsIn)[names(resultsIn) == 'GeoPrediction1'] <- 'GeoPrediction'
  names(resultsIn)[names(resultsIn) == 'Malaria.Species'] <- 'Malaria_Species'
  names(resultsIn)[names(resultsIn) == 'hrpStatus'] <- 'HRP2/3 Deletions'
  names(resultsIn)[names(resultsIn) == 'Final_COI_predction..Mono.multi.strain.infection.'] <- 'Multiplicity of Infection'
  resultsIn$Malaria_Species[is.na(resultsIn$Malaria_Species)] <- "N/A"
  resultsIn$Malaria_Species[resultsIn$Malaria_Species == ""] <- "N/A"
  resultsIn$Sample_Status[resultsIn$Sample_Status == ""] <- "N/A"
  resultsIn$Sample_Status[is.na(resultsIn$Sample_Status)] <- "N/A"

  resultsIn$Reported_Travel_History[is.na(resultsIn$Reported_Travel_History)] <- "Not Provided"
  resultsIn$Reported_Travel_History[resultsIn$Reported_Travel_History == ""] <- "Not Provided"
  resultsIn$Reported_Travel_History[resultsIn$Reported_Travel_History == "NOT PROVIDED"] <- "Not Provided"
  resultsIn$Reported_Travel_History[resultsIn$Reported_Travel_History == "not provided"] <- "Not Provided"
  print("check in colnames")
  print(colnames(resultsIn))

#    [1] "Sample_Status"                                    
#  [2] "CSID"                                             
#  [3] "LSDB_Specimen_ID"                                 
#  [4] "LSDB_Sequence_ID"                                 
#  [5] "State_Lab_ID"                                     
#  [6] "Case_ID"                                          
#  [7] "Collection_date"                                  
#  [8] "Submitter_State"                                  
#  [9] "Reported_Travel_History"                          
# [10] "Malaria_Species"                                  
# [11] "HRP2/3 Deletions"                                 
# [12] "GeoPrediction"                                    
# [13] "GeneticClusterNumber"                             
# [14] "Final_COI_predction..Mono.multi.strain.infection."
# [15] "State_name"                                       
# [16] "ReportStatusDate"      
  #stop("after results in colnames")
  # stop("results in colnames")
#   Sample_Status"           "CSID"                   
#  [3] "LSDB_Specimen_ID"        "LSDB_Sequence_ID"       
#  [5] "State_Lab_ID"            "Collection_date"        
#  [7] "Submitter_State"         "Reported_Travel_History"
#  [9] "Malaria_Species"         "hrpStatus"              
# [11] "GeoPrediction"           "State_name"             
# [13] "ReportStatusDate"        "Pfcrt.C72S

  resultsForTable1 <- resultsIn[,c(3,2,5,6,7,9,1,16,10,11,12,14)]
  # resultsForTable1 <- resultsIn[,c(,2,6,5,4,1,12,7,10)]

  resultsForTable1_noGeo <- resultsIn[,c(3,5,6,7,9,1,16,10)]

  #make a name for this file when it is eventually printed (necessary for including as download link in html file)
  specimenProcessingCSV <- paste(fullPath,date_now, "_", time_now,"_", stateID, "_malaria_sample_processing_summary.csv", sep ="")

  print("before table 1")
  #Make a table summarizing everything we've received at CDC / status / species/ geo prediction
  print(colnames(resultsForTable1_noGeo))
  table1 <- knitr::kable(resultsForTable1_noGeo, align = "c",caption = "Table 2: Samples Processed for Domestic Malaria Surveillance") %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12) %>%
      add_header_above(c("CDC IDs" = 1,"Sample Metadata" = 4, "CDC Processing" = 3),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", color = "black", background = "#C6E2FF", extra_css = "line-height: 20pt;") %>%
      column_spec(c(1,5), border_right = "2px solid gray")

  #grab all genotypes plust limited metadata (used for p fal reporting)
  # completeSNPTable <- resultsIn[c(3,5,6,9,11,12,14,17:83)] 
  completeSNPTable <- resultsIn[c(3,5,6,9,11,12,14,17:67)] 
  print(colnames(completeSNPTable)) 
  completeMaRSName <- paste(fullPath,date_now, "_", time_now,"_", stateID, "_MaRS_genotyping_results.csv", sep ="")

  #now read in the CRT results table from Dhruvi's code. This is first used to see if the state has any specimens processed in MaRS
  crttable_file <- paste(fullPath,"input_data/",stateID,"_crtTable.csv", sep ="")
  # test <-paste(fullPath,"input_data/",stateID,"_crtTable.csv", sep ="")
  PFCRT_copy = read.csv(crttable_file, header = T, stringsAsFactors=FALSE, colClasses = c("character"))
  

  PFCRT_copy <- PFCRT_copy %>% rename_with(~ str_replace(., 'Pfcrt.', ''))
  print("check 1")
#Table of specimens analyzed in MaRS
  if(length(rownames(PFCRT_copy))>0){
    analyzedMars <- PFCRT_copy$CSID
    #get samples anlyzed in mars from the larger table of all specimens submitted by the state
    marsTableStart <- resultsForTable1[resultsForTable1$CSID %in% analyzedMars,]
    rownames(marsTableStart) <- NULL
    marsTableStart <- marsTableStart[,c(1,3:12)]
    marsTable <- knitr::kable(marsTableStart, align = "c", caption = "Table 3: Samples Analyzed in the MaRS Bioinformatic Workflow") %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12) %>%
      add_header_above(c("CDC IDs" = 1,"Sample Metadata" = 4, "CDC Processing" = 6),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", color = "black", background = "#C6E2FF", extra_css = "line-height: 20pt;") %>%
      column_spec(c(1,5), border_right = "2px solid gray")

   }else{
    #If no specimens processed in MaRS, make an empty table
    marsTableStart <- data.frame(matrix(ncol = 11, nrow = 1))
    marsTableStart[1,1] <- "No Samples Tested"
    marsTableStart[1,c(2:11)] <- "N/A"
    resultsForTable1 <-resultsForTable1[,c(1,3:12)]
    colnames(marsTableStart) <- colnames(resultsForTable1)
    marsTable <- knitr::kable(marsTableStart, align = "c", caption = "Table 3: Samples Analyzed in the MaRS Bioinformatic Workflow") %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12) %>%
      add_header_above(c("No Samples Processed for P. falciparum Drug Resistance Screening" = 11), italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", background = "#FF7F00", color = "black", extra_css = "line-height: 20pt;")
     
  }

  PFCRT_copy <- PFCRT_copy[,c(3,2,6,5,7,9,11,14,17,18,19,20,21)]
  if(length(rownames(PFCRT_copy)) == 0){
    PFCRT_copy[1,c(1:14)] <- "N/A"
  }
  
#Loop through each specimen's results for crt and see if has wildtype or alternative allele. This is used to color the cells
#Grab/order metadata and the  mars results
print(colnames(PFCRT_copy))

  # PFCRT_copy <- PFCRT_copy[,c(2,1,4,5,3,6,7,9,10:14)]
  

  names(PFCRT_copy)[names(PFCRT_copy) == 'Patient_Travel_History_Interpretation'] <- 'Reported_Travel_History'
  names(PFCRT_copy)[names(PFCRT_copy) == 'GeoPrediction1'] <- 'GeoPrediction'
  names(PFCRT_copy)[names(PFCRT_copy) == 'hrpStatus'] <- 'HRP2/3 Deletions'
  names(PFCRT_copy)[names(PFCRT_copy) == 'Final_COI_predction..Mono.multi.strain.infection.'] <- 'Multiplicity of Infection'
  PFCRT_copy  <- PFCRT_copy[order(PFCRT_copy$LSDB_Specimen_ID),]
  rownames(PFCRT_copy) <- NULL
  colVec6 <- c()
  colVec7 <- c()
  colVec8 <- c()
  colVec9 <- c()
  colVec10 <- c()
  print("check 2")
  #Repeat this loop for all markers (in other parts of this code). First if statement is the wildtype, if is the wild type, don't change color. 
  #Second statement checks for x (no amplification), color this blue
  #Third statement check for N/A (no samples analyzed), leave this color as is
  #Finally, if it is none of the above, it must be alternate allele. Color it red
  addCol <- NULL
  print(PFCRT_copy)
  for(k in 1:length(rownames(PFCRT_copy))){
    if(PFCRT_copy[k,9] == "C"){
      addCol <- ""
    }else if(PFCRT_copy[k,9] == "x"){
      addCol <- "#A4D3EE"
    }else if(PFCRT_copy[k,9] == "N/A"){
      addCol <- ""
    }
    else{
      addCol <- "#FFA07A"
    }
    colVec6 <- c(colVec6, addCol)
  }

  addCol <- NULL
  for(k in 1:length(rownames(PFCRT_copy))){
    if(PFCRT_copy[k,10] == "V"){
      addCol <- ""
    }else if(PFCRT_copy[k,10] == "x"){
      addCol <- "#A4D3EE"
    } else if(PFCRT_copy[k,10] == "N/A"){
      addCol <- ""
    }
    else{
      addCol <- "#FFA07A"
    }
    colVec7 <- c(colVec7, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(PFCRT_copy))){
    if(PFCRT_copy[k,11] == "M"){
      addCol <- ""
    }else if(PFCRT_copy[k,11] == "x"){
      print("is blue")
      addCol <- "#A4D3EE"
    } else if(PFCRT_copy[k,11] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec8 <- c(colVec8, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(PFCRT_copy))){
    if(PFCRT_copy[k,12] == "N"){
      addCol <- ""
    }else if(PFCRT_copy[k,12] == "x"){
      print("is blue")
      addCol <- "#A4D3EE"
    } else if(PFCRT_copy[k,12] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec9 <- c(colVec9, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(PFCRT_copy))){
    if(PFCRT_copy[k,13] == "K"){
      addCol <- ""
    }else if(PFCRT_copy[k,13] == "x"){
      print("is blue")
      addCol <- "#A4D3EE"
    } else if(PFCRT_copy[k,13] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec10 <- c(colVec10, addCol)
  }
  print(colVec10)
#Now make the crt table. Use the color column vectors for each genotype column
  PFCRT_copy <- PFCRT_copy[,c(1,3:13)]
  
  print(PFCRT_copy)
  crtTable <- knitr::kable(PFCRT_copy, align = "c",caption = "Table 4: Genetic markers in <i>Pfcrt</i> associated with resistance to the antimalarial drug chloroquine", escape= FALSE) %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12, protect_latex = TRUE) %>%
      add_header_above(c("Sample Metadata" = 7, "Pfcrt SNPs" = 5),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", color = "black", background = "#CDCDC1", extra_css = "line-height: 20pt;")  %>%
      column_spec(c(7), border_right = "2px solid gray") %>%
      column_spec(c(8), background = colVec6) %>%
      column_spec(c(9), background = colVec7) %>%
      column_spec(c(10), background = colVec8) %>%
      column_spec(c(11), background = colVec9) %>%
      column_spec(c(12), background = colVec10)


############# Now take the same process that we used for Crt and apply it to K13, mdr1, cytb, and dhps/dhfr
  k13table_file <- paste(fullPath,"input_data/",stateID,"_k13Table.csv", sep ="")
  K13 = read.csv(k13table_file, header = T, , stringsAsFactors=FALSE, colClasses = c("character"))
  K13 <- K13 %>% rename_with(~ str_replace(., 'Pfk13.', ''))
  print(colnames(K13))
  K13 <- K13[,c(3,2,6,5,7,9,11,14,18,20,21,22,38,27,28,29,30,32,34,36,37)]
  if(length(rownames(K13)) == 0){
    K13[1,c(1:22)] <- "N/A"
  }
  names(K13)[names(K13) == 'Patient_Travel_History_Interpretation'] <- 'Reported_Travel_History'
  names(K13)[names(K13) == 'GeoPrediction1'] <- 'GeoPrediction'
  names(K13)[names(K13) == 'hrpStatus'] <- 'HRP2/3 Deletions'
  names(K13)[names(K13) == 'Final_COI_predction..Mono.multi.strain.infection.'] <- 'Multiplicity of Infection'
  K13 <- K13[order(K13$LSDB_Specimen_ID),]
  rownames(K13) <- NULL
  # K13<- K13[,c(2,1,4,5,3,6,7,9, 10:22)]
  
  colVec6 <- c()
  colVec7 <- c()
  colVec8 <- c()
  colVec9 <- c()
  colVec10 <- c()
  colVec11 <- c()
  colVec12 <- c()
  colVec13 <- c()
  colVec14 <- c()
  colVec15 <- c()
  colVec16 <- c()
  colVec17 <- c()
  colVec18 <- c()

  print("check 3")
  print(K13)
  addCol <- NULL
  for(k in 1:length(rownames(K13))){
    if(K13[k,9] == "F"){
      addCol <- ""
    }else if(K13[k,9] == "x"){
      addCol <- "#A4D3EE"
    }else if(K13[k,9] == "N/A"){
      addCol <- ""
    }
    else{
      addCol <- "#FFA07A"
    }
    colVec6 <- c(colVec6, addCol)
  }

  addCol <- NULL
  for(k in 1:length(rownames(K13))){
    if(K13[k,10] == "N"){
      addCol <- ""
    }else if(K13[k,10] == "x"){
      addCol <- "#A4D3EE"
    } else if(K13[k,10] == "N/A"){
      addCol <- ""
    }
    else{
      addCol <- "#FFA07A"
    }
    colVec7 <- c(colVec7, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(K13))){
    if(K13[k,11] == "C"){
      addCol <- ""
    }else if(K13[k,11] == "x"){
      addCol <- "#A4D3EE"
    } else if(K13[k,11] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec8 <- c(colVec8, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(K13))){
    if(K13[k,12] == "M"){
      addCol <- ""
    }else if(K13[k,12] == "x"){
      addCol <- "#A4D3EE"
    } else if(K13[k,12] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec9 <- c(colVec9, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(K13))){
    if(K13[k,13] == "Y"){
      addCol <- ""
    }else if(K13[k,13] == "x"){
      addCol <- "#A4D3EE"
    } else if(K13[k,13] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec10 <- c(colVec10, addCol)
  }
  for(k in 1:length(rownames(K13))){
    if(K13[k,14] == "R"){
      addCol <- ""
    }else if(K13[k,14] == "x"){
      addCol <- "#A4D3EE"
    } else if(K13[k,14] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec11 <- c(colVec11, addCol)
  }
  for(k in 1:length(rownames(K13))){
    if(K13[k,15] == "I"){
      addCol <- ""
    }else if(K13[k,15] == "x"){
      addCol <- "#A4D3EE"
    } else if(K13[k,15] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec12 <- c(colVec12, addCol)
  }
  for(k in 1:length(rownames(K13))){
    if(K13[k,16] == "P"){
      addCol <- ""
    }else if(K13[k,16] == "x"){
      addCol <- "#A4D3EE"
    } else if(K13[k,16] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec13 <- c(colVec13, addCol)
  }
  for(k in 1:length(rownames(K13))){
    if(K13[k,17] == "R"){
      addCol <- ""
    }else if(K13[k,17] == "x"){
      addCol <- "#A4D3EE"
    } else if(K13[k,17] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec14 <- c(colVec14, addCol)
  }
  for(k in 1:length(rownames(K13))){
    if(K13[k,18] == "P"){
      addCol <- ""
    }else if(K13[k,18] == "x"){
      addCol <- "#A4D3EE"
    } else if(K13[k,18] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec15 <- c(colVec15, addCol)
  }
  for(k in 1:length(rownames(K13))){
    if(K13[k,19] == "C"){
      addCol <- ""
    }else if(K13[k,19] == "x"){
      addCol <- "#A4D3EE"
    } else if(K13[k,19] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec16 <- c(colVec16, addCol)
  }
  for(k in 1:length(rownames(K13))){
    if(K13[k,20] == "R"){
      addCol <- ""
    }else if(K13[k,20] == "x"){
      addCol <- "#A4D3EE"
    } else if(K13[k,20] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec17 <- c(colVec17, addCol)
  }
  for(k in 1:length(rownames(K13))){
    if(K13[k,21] == "A"){
      addCol <- ""
    }else if(K13[k,21] == "x"){
      addCol <- "#A4D3EE"
    } else if(K13[k,21] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec18 <- c(colVec18, addCol)
  }
  K13 <- K13[,c(1,3:21)]
  
  k13Table <- knitr::kable(K13, align = "c",caption = "Table 5: Genetic markers in <i>PfK13</i> associated with resistance to the antimalarial drug artemisinin ") %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12) %>%
    add_header_above(c("Sample Meatadata" = 7, "PfK13 SNPs " = 13),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", color = "black", background = "#CDCDC1", extra_css = "line-height: 20pt;") %>%
    column_spec(c(7), border_right = "2px solid gray") %>%
    column_spec(c(8), background = colVec6) %>%
    column_spec(c(9), background = colVec7) %>%
    column_spec(c(10), background = colVec8) %>%
    column_spec(c(11), background = colVec9) %>%
    column_spec(c(12), background = colVec10) %>%
    column_spec(c(13), background = colVec11) %>%
    column_spec(c(14), background = colVec12) %>%
    column_spec(c(15), background = colVec13) %>%
    column_spec(c(16), background = colVec14) %>%
    column_spec(c(17), background = colVec15) %>%
    column_spec(c(18), background = colVec16) %>%
    column_spec(c(19), background = colVec17) %>%
    column_spec(c(20), background = colVec18)

  #MDR1
  mdr1table_file <- paste(fullPath,"input_data/",stateID,"_mdr1Table.csv", sep ="")
  PFMDR1 = read.csv(mdr1table_file, header = T, stringsAsFactors=FALSE, colClasses = c("character"))
  PFMDR1 <- PFMDR1 %>% rename_with(~ str_replace(., 'Pfmdr1.', ''))
  print("checkmdr names")
  print(colnames(PFMDR1))
  PFMDR1<- PFMDR1[,c(3,2,6,5,7,9,11,14,17,18,19,20,21)]
  if(length(rownames(PFMDR1)) == 0){
    PFMDR1[1,c(1:14)] <- "N/A"
  }
  names(PFMDR1)[names(PFMDR1) == 'Patient_Travel_History_Interpretation'] <- 'Reported_Travel_History'
  names(PFMDR1)[names(PFMDR1) == 'GeoPrediction1'] <- 'GeoPrediction'
  names(PFMDR1)[names(PFMDR1) == 'hrpStatus'] <- 'HRP2/3 Deletions'
  names(PFMDR1)[names(PFMDR1) == 'Final_COI_predction..Mono.multi.strain.infection.'] <- 'Multiplicity of Infection'
  PFMDR1  <- PFMDR1[order(PFMDR1$LSDB_Specimen_ID),]
  rownames(PFMDR1) <- NULL
  # PFMDR1 <- PFMDR1[,c(2,1,4,5,3,6,7,9,10:14)]

  colVec6 <- c()
  colVec7 <- c()
  colVec8 <- c()
  colVec9 <- c()
  colVec10 <- c()  
  
  addCol <- NULL
  for(k in 1:length(rownames(PFMDR1))){
    if(PFMDR1[k,9] == "N"){
      addCol <- ""
    }else if(PFMDR1[k,9] == "x"){
      addCol <- "#A4D3EE"
    }else if(PFMDR1[k,9] == "N/A"){
      addCol <- ""
    }
    else{
      addCol <- "#FFA07A"
    }
    colVec6 <- c(colVec6, addCol)
  }

  addCol <- NULL
  for(k in 1:length(rownames(PFMDR1))){
    if(PFMDR1[k,10] == "Y"){
      addCol <- ""
    }else if(PFMDR1[k,10] == "x"){
      addCol <- "#A4D3EE"
    } else if(PFMDR1[k,10] == "N/A"){
      addCol <- ""
    }
    else{
      addCol <- "#FFA07A"
    }
    colVec7 <- c(colVec7, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(PFMDR1))){
    if(PFMDR1[k,11] == "S"){
      addCol <- ""
    }else if(PFMDR1[k,11] == "x"){
      addCol <- "#A4D3EE"
    } else if(PFMDR1[k,11] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec8 <- c(colVec8, addCol)
  }
   addCol <- NULL
  for(k in 1:length(rownames(PFMDR1))){
    if(PFMDR1[k,12] == "N"){
      addCol <- ""
    }else if(PFMDR1[k,12] == "x"){
      addCol <- "#A4D3EE"
    } else if(PFMDR1[k,12] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec9 <- c(colVec9, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(PFMDR1))){
    if(PFMDR1[k,13] == "D"){
      addCol <- ""
    }else if(PFMDR1[k,13] == "x"){
      addCol <- "#A4D3EE"
    } else if(PFMDR1[k,13] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec10 <- c(colVec10, addCol)
  }

  PFMDR1 <- PFMDR1[,c(1,3:13)]
  
  mdr1Table <- knitr::kable(PFMDR1, align = "c",caption = "Table 6: Genetic markers in <i>Pfmdr1</i> associated with resistance to the antimalarial drug lumefantrine") %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12) %>%
    add_header_above(c("Sample Metadata" = 7, "Pfmdr1 SNPs" = 5),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", color = "black", background = "#CDCDC1",extra_css = "line-height: 20pt;")  %>%
    column_spec(c(7), border_right = "2px solid gray") %>%
    column_spec(c(8), background = colVec6) %>%
    column_spec(c(9), background = colVec7) %>%
    column_spec(c(10), background = colVec8) %>%
    column_spec(c(11), background = colVec9) %>%
    column_spec(c(12), background = colVec10)



  #Cytb
  cytb_file <- paste(fullPath,"input_data/",stateID,"_cytbTable.csv", sep ="")
  pfcytb = read.csv(cytb_file, header = T, stringsAsFactors=FALSE, colClasses = c("character"))
  pfcytb <- pfcytb %>% rename_with(~ str_replace(., 'Pfcytb.', ''))
  pfcytb<- pfcytb[,c(3,2,6,5,7,9,11,14,17,18)]
  if(length(rownames(pfcytb)) == 0){
    pfcytb[1,c(1:12)] <- "N/A"
  }
  names(pfcytb)[names(pfcytb) == 'Patient_Travel_History_Interpretation'] <- 'Reported_Travel_History'
  names(pfcytb)[names(pfcytb) == 'GeoPrediction1'] <- 'GeoPrediction'
  names(pfcytb)[names(pfcytb) == 'hrpStatus'] <- 'HRP2/3 Deletions'
  names(pfcytb)[names(pfcytb) == 'Final_COI_predction..Mono.multi.strain.infection.'] <- 'Multiplicity of Infection'
  # pfcytb <- pfcytb[,c(2,1,4,5,3,6,7,9,10:12)]
  pfcytb  <- pfcytb[order(pfcytb$LSDB_Specimen_ID),]
  rownames(pfcytb) <- NULL

  colVec6 <- c()
  colVec7 <- c()
  colVec8 <- c()
  print(pfcytb)
  addCol <- NULL
  for(k in 1:length(rownames(pfcytb))){
    if(pfcytb[k,9] == "I"){
      addCol <- ""
    }else if(pfcytb[k,9] == "x"){
      addCol <- "#A4D3EE"
    }else if(pfcytb[k,9] == "N/A"){
      addCol <- ""
    }
    else{
      addCol <- "#FFA07A"
    }
    colVec6 <- c(colVec6, addCol)
  }

  addCol <- NULL
  for(k in 1:length(rownames(pfcytb))){
    if(pfcytb[k,10] == "Y"){
      addCol <- ""
    }else if(pfcytb[k,10] == "x"){
      addCol <- "#A4D3EE"
    } else if(pfcytb[k,10] == "N/A"){
      addCol <- ""
    }
    else{
      addCol <- "#FFA07A"
    }
    colVec7 <- c(colVec7, addCol)
  }
  # addCol <- NULL
  # for(k in 1:length(rownames(pfcytb))){
  #   if(pfcytb[k,11] == "Y"){
  #     addCol <- ""
  #   }else if(pfcytb[k,11] == "x"){
  #     addCol <- "#A4D3EE"
  #   } else if(pfcytb[k,11] == "N/A"){
  #     addCol <- ""
  #   }else{
  #     addCol <- "#FFA07A"
  #   }
  #   colVec8 <- c(colVec8, addCol)
  # }
  print(colnames(pfcytb))
  colnames(pfcytb)[10] <- "Y268S/C"
  pfcytb <- pfcytb[c(1,3,4,5,6,7,8,9,10)]
 
  cytbTable <- knitr::kable(pfcytb, align = "c",caption = "Table 7: Genetic markers in <i>Pfcytb</i> associated with resistance to the antimalarial drug atovaquone") %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12) %>%
    add_header_above(c("Sample Metadata" = 7, "Pfcytb SNPs" = 2),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", color = "black", background = "#CDCDC1", extra_css = "line-height: 20pt;") %>%
    column_spec(c(7), border_right = "2px solid gray") %>%
    column_spec(c(8), background = colVec6) %>%
    column_spec(c(9), background = colVec7)
    # column_spec(c(10), background = colVec8)

  #interpretation


  #dhps / dhfr
  #Table
  pfdhpsDhfrTable_file <- paste(fullPath,"input_data/",stateID,"_dhpsDhfrTable.csv", sep ="")
  pfdhps_dhfr_copy = read.csv(pfdhpsDhfrTable_file, header = T, stringsAsFactors=FALSE, colClasses = c("character"))
  print("after read dhps")

  print(pfdhps_dhfr_copy)
  pfdhps_dhfr_copy<- pfdhps_dhfr_copy[,c(3,2,6,5,7,9,11,14,17,18,19,20,21,22,23,24,25)]

  if(length(rownames(pfdhps_dhfr_copy)) == 0){
    pfdhps_dhfr_copy[1,c(1:18)] <- "N/A"
  }
  
  names(pfdhps_dhfr_copy)[names(pfdhps_dhfr_copy) == 'Patient_Travel_History_Interpretation'] <- 'Reported_Travel_History'
  names(pfdhps_dhfr_copy)[names(pfdhps_dhfr_copy) == 'GeoPrediction1'] <- 'GeoPrediction'
  names(pfdhps_dhfr_copy)[names(pfdhps_dhfr_copy) == 'hrpStatus'] <- 'HRP2/3 Deletions'
  names(pfdhps_dhfr_copy)[names(pfdhps_dhfr_copy) == 'Final_COI_predction..Mono.multi.strain.infection.'] <- 'Multiplicity of Infection'
  # pfdhps_dhfr_copy <- pfdhps_dhfr_copy[,c(2,1,4,5,3,6,7,9,10:18)]
  pfdhps_dhfr_copy  <- pfdhps_dhfr_copy[order(pfdhps_dhfr_copy$LSDB_Specimen_ID),]
  rownames(pfdhps_dhfr_copy) <- NULL
  # pfdhps_dhfr_copy <- pfdhps_dhfr_copy[,c(1,4,3,5,8:16)]

  colVec6 <- c()
  colVec7 <- c()
  colVec8 <- c()
  colVec9 <- c()
  colVec10 <- c()
  colVec11 <- c()
  colVec12 <- c()
  colVec13 <- c()
  colVec14 <- c()

  addCol <- NULL
  for(k in 1:length(rownames(pfdhps_dhfr_copy))){
    if(pfdhps_dhfr_copy[k,9] == "N"){
      addCol <- ""
    }else if(pfdhps_dhfr_copy[k,9] == "x"){
      addCol <- "#A4D3EE"
    }else if(pfdhps_dhfr_copy[k,9] == "N/A"){
      addCol <- ""
    }
    else{
      addCol <- "#FFA07A"
    }
    colVec6 <- c(colVec6, addCol)
  }

  addCol <- NULL
  for(k in 1:length(rownames(pfdhps_dhfr_copy))){
    if(pfdhps_dhfr_copy[k,10] == "C"){
      addCol <- ""
    }else if(pfdhps_dhfr_copy[k,10] == "x"){
      addCol <- "#A4D3EE"
    } else if(pfdhps_dhfr_copy[k,10] == "N/A"){
      addCol <- ""
    }
    else{
      addCol <- "#FFA07A"
    }
    colVec7 <- c(colVec7, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(pfdhps_dhfr_copy))){
    if(pfdhps_dhfr_copy[k,11] == "S"){
      addCol <- ""
    }else if(pfdhps_dhfr_copy[k,11] == "x"){
      addCol <- "#A4D3EE"
    } else if(pfdhps_dhfr_copy[k,11] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec8 <- c(colVec8, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(pfdhps_dhfr_copy))){
    if(pfdhps_dhfr_copy[k,12] == "I"){
      addCol <- ""
    }else if(pfdhps_dhfr_copy[k,12] == "x"){
      addCol <- "#A4D3EE"
    } else if(pfdhps_dhfr_copy[k,12] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec9 <- c(colVec9, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(pfdhps_dhfr_copy))){
    if(pfdhps_dhfr_copy[k,13] == "S"){
      addCol <- ""
    }else if(pfdhps_dhfr_copy[k,13] == "x"){
      addCol <- "#A4D3EE"
    } else if(pfdhps_dhfr_copy[k,13] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec10 <- c(colVec10, addCol)
  }
  for(k in 1:length(rownames(pfdhps_dhfr_copy))){
    if(pfdhps_dhfr_copy[k,14] == "A"){
      addCol <- ""
    }else if(pfdhps_dhfr_copy[k,14] == "x"){
      addCol <- "#A4D3EE"
    } else if(pfdhps_dhfr_copy[k,14] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec11 <- c(colVec11, addCol)
  }
  for(k in 1:length(rownames(pfdhps_dhfr_copy))){
    if(pfdhps_dhfr_copy[k,15] == "K"){
      addCol <- ""
    }else if(pfdhps_dhfr_copy[k,15] == "x"){
      addCol <- "#A4D3EE"
    } else if(pfdhps_dhfr_copy[k,15] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec12 <- c(colVec12, addCol)
  }
  for(k in 1:length(rownames(pfdhps_dhfr_copy))){
    if(pfdhps_dhfr_copy[k,16] == "A"){
      addCol <- ""
    }else if(pfdhps_dhfr_copy[k,16] == "x"){
      addCol <- "#A4D3EE"
    } else if(pfdhps_dhfr_copy[k,16] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec13 <- c(colVec13, addCol)
  }
  for(k in 1:length(rownames(pfdhps_dhfr_copy))){
    if(pfdhps_dhfr_copy[k,17] == "A"){
      addCol <- ""
    }else if(pfdhps_dhfr_copy[k,17] == "x"){
      addCol <- "#A4D3EE"
    } else if(pfdhps_dhfr_copy[k,17] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec14 <- c(colVec14, addCol)
  }
  # print("checck dhps colnames")
  # print(colnames(pfdhps_dhfr_copy))
  print("check dhfr colnames")
  print(colnames(pfdhps_dhfr_copy))
  # stop("dfhr colnames")
  # pfdhps_dhfr_copy <- pfdhps_dhfr_copy[,c(1,4,3,5,6:15)]
  # pfdhps_dhfr_copy <- pfdhps_dhfr_copy[,c(2,3,4,5,6,7,8,9:17)]
  # pfdhps_dhfr_copy <- pfdhps_dhfr_copy[,c(1,5,3,4,7,6,8,9:17)]
  pfdhps_dhfr_copy <- pfdhps_dhfr_copy[,c(1,3:17)]
  colnames(pfdhps_dhfr_copy)[12] <- "Pfdhps.S436A/H"
  colnames(pfdhps_dhfr_copy)[16] <- "Pfdhps.A613S/T"

  pfdhps_dhfr_copy  <- pfdhps_dhfr_copy  %>% rename_with(~ str_replace(., 'Pfdhps.', ''))
  pfdhps_dhfr_copy  <- pfdhps_dhfr_copy  %>% rename_with(~ str_replace(., 'Pfdhfr.', ''))
  
  pfdhps_dhfr_copyKbl <- knitr::kable(pfdhps_dhfr_copy, align = "c",caption ="Table 8: Genetic markers in <i>Pfdhfr</i> are associated with resistance to pyrimethamine and cycloguanil/proguanil and genetic markers in <i>Pfdhps</i> are associated with resistance to the antimalarial drug sulfadoxine.") %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12) %>%
    add_header_above(c("Sample Metadata" = 7, "Pfdhfr SNPs" = 3, "Pfdhps SNPs"= 6 ),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray",  background = "#CDCDC1", color = "black", extra_css = "line-height: 20pt;") %>%
    column_spec(c(7,10), border_right = "2px solid gray") %>%
    column_spec(c(8), background = colVec6) %>%
    column_spec(c(9), background = colVec7) %>%
    column_spec(c(10), background = colVec8) %>%
    column_spec(c(11), background = colVec9) %>%
    column_spec(c(12), background = colVec10) %>%
    column_spec(c(13), background = colVec11) %>%
    column_spec(c(14), background = colVec12) %>%
    column_spec(c(15), background = colVec13) %>%
    column_spec(c(16), background = colVec14) 

  
  ##plasmepsin placeholder
  plasmepsinTable_in <- PFCRT_copy
  plasmepsinTable_in[,c(8:12)] <- "N/A"
  colnames(plasmepsinTable_in)[8:12] <- "TBD"
  print("before plasmepsin")
  plasmepsinTable_in  <- plasmepsinTable_in[order(plasmepsinTable_in$LSDB_Specimen_ID),]
  rownames(plasmepsinTable_in) <- NULL
  plasmepsinTable <- knitr::kable(plasmepsinTable_in, align = "c",caption = "Table 9: Genetic markers in <i>plasmpesin 2</i> associated with resistance to the antimalarial drug piperaquine.", escape= FALSE) %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12, protect_latex = TRUE) %>%
      add_header_above(c("Sample Metadata" = 7, "Pfplasmepsin SNPs" = 5),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", color = "black", background = "#CDCDC1", extra_css = "line-height: 20pt;")  %>%
      column_spec(c(7), border_right = "2px solid gray")
      # column_spec(c(6), background = colVec6) %>%
      # column_spec(c(7), background = colVec7) %>%
      # column_spec(c(8), background = colVec8) %>%
      # column_spec(c(9), background = colVec9) %>%
      # column_spec(c(10), background = colVec10)
  print("before coronin")

  # ##coronin placeholder
  # coroninTable_in <- PFCRT_copy
  # coroninTable_in[,c(7:12)] <- "N/A"
  # colnames(coroninTable_in)[7:12] <- "TBD"
  # coroninTable_in  <- coroninTable_in[order(coroninTable_in$LSDB_Specimen_ID),]
  # rownames(coroninTable_in) <- NULL
  # coroninTable <- knitr::kable(coroninTable_in, align = "c",caption = "Table 10: Genetic markers in <i>Pfcoronin</i> associated with resistance to artemisinin based therapies.", escape= FALSE) %>%
  #     kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12, protect_latex = TRUE) %>%
  #     add_header_above(c("Sample Metadata" = 7, "Pfcoronin SNPs" = 5),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", color = "black", background = "#CDCDC1", extra_css = "line-height: 20pt;")  %>%
  #     column_spec(c(7), border_right = "2px solid gray")
  #     # column_spec(c(6), background = colVec6) %>%
  #     # column_spec(c(7), background = colVec7) %>%
  #     # column_spec(c(8), background = colVec8) %>%
  #     # column_spec(c(9), background = colVec9) %>%
  #     # column_spec(c(10), background = colVec10)
  # print("after coronin")
  # #MDR1

  coronintable_file <- paste(fullPath,"input_data/",stateID,"_coroninTable.csv", sep ="")
  PFCORONIN = read.csv(coronintable_file, header = T, stringsAsFactors=FALSE, colClasses = c("character"))
  PFCORONIN <- PFCORONIN %>% rename_with(~ str_replace(., 'Pfcoronin.', ''))

  if(length(rownames(PFCORONIN)) == 0){
    PFCORONIN[1,c(1:14)] <- "N/A"
  }
  names(PFCORONIN)[names(PFCORONIN) == 'Patient_Travel_History_Interpretation'] <- 'Reported_Travel_History'
  names(PFCORONIN)[names(PFCORONIN) == 'GeoPrediction1'] <- 'GeoPrediction'
  names(PFCORONIN)[names(PFCORONIN) == 'hrpStatus'] <- 'HRP2/3 Deletions'
  names(PFCORONIN)[names(PFCORONIN) == 'Final_COI_predction..Mono.multi.strain.infection.'] <- 'Multiplicity of Infection'

  PFCORONIN  <- PFCORONIN[order(PFCORONIN$LSDB_Specimen_ID),]
  rownames(PFCORONIN) <- NULL
  PFCORONIN <- PFCORONIN[,c(2,1,4,5,3,6,7,9,10:14)]

  colVec6 <- c()
  colVec7 <- c()
  colVec8 <- c()
  colVec9 <- c()
  colVec10 <- c()  
  
  addCol <- NULL
  for(k in 1:length(rownames(PFCORONIN))){
    if(PFCORONIN[k,9] == "G"){
      addCol <- ""
    }else if(PFCORONIN[k,9] == "x"){
      addCol <- "#A4D3EE"
    }else if(PFCORONIN[k,9] == "N/A"){
      addCol <- ""
    }
    else{
      addCol <- "#FFA07A"
    }
    colVec6 <- c(colVec6, addCol)
  }

  addCol <- NULL
  for(k in 1:length(rownames(PFCORONIN))){
    if(PFCORONIN[k,10] == "V"){
      addCol <- ""
    }else if(PFCORONIN[k,10] == "x"){
      addCol <- "#A4D3EE"
    } else if(PFCORONIN[k,10] == "N/A"){
      addCol <- ""
    }
    else{
      addCol <- "#FFA07A"
    }
    colVec7 <- c(colVec7, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(PFCORONIN))){
    if(PFCORONIN[k,11] == "P"){
      addCol <- ""
    }else if(PFCORONIN[k,11] == "x"){
      addCol <- "#A4D3EE"
    } else if(PFCORONIN[k,11] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec8 <- c(colVec8, addCol)
  }
   addCol <- NULL
  for(k in 1:length(rownames(PFCORONIN))){
    if(PFCORONIN[k,12] == "R"){
      addCol <- ""
    }else if(PFCORONIN[k,12] == "x"){
      addCol <- "#A4D3EE"
    } else if(PFCORONIN[k,12] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec9 <- c(colVec9, addCol)
  }
  addCol <- NULL
  for(k in 1:length(rownames(PFCORONIN))){
    if(PFCORONIN[k,13] == "E"){
      addCol <- ""
    }else if(PFCORONIN[k,13] == "x"){
      addCol <- "#A4D3EE"
    } else if(PFCORONIN[k,13] == "N/A"){
      addCol <- ""
    }else{
      addCol <- "#FFA07A"
    }
    colVec10 <- c(colVec10, addCol)
  }

  PFCORONIN <- PFCORONIN[,c(1,3:13)]
  
  coroninTable <- knitr::kable(PFCORONIN, align = "c",caption = "Table 10: Genetic markers in <i>Pfcoronin</i> associated with resistance to artemisinin based therapies.", escape= FALSE) %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12, protect_latex = TRUE) %>%
      add_header_above(c("Sample Metadata" = 7, "Pfcoronin SNPs" = 5),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", color = "black", background = "#CDCDC1", extra_css = "line-height: 20pt;")  %>%
      column_spec(c(7), border_right = "2px solid gray") %>%
      column_spec(c(8), background = colVec6) %>%
      column_spec(c(9), background = colVec7) %>%
      column_spec(c(10), background = colVec8) %>%
      column_spec(c(11), background = colVec9) %>%
      column_spec(c(12), background = colVec10)

  
  ##Now finished with the general / mars tables and can make the Pvivax report (still in the state loop)
  # setwd('..')
  knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, collapse= TRUE)


  ##Now start on the pvivax reports. Most of the preprocessing is done prior to the loop
  # merge3 <- merge2[grepl(stateID, merge2$Seq_ID),]
  #now need to find where genetic cluster number is

  print("check making pvivax table")
  print((stateID))
  print(head(resultsIn))
  # resultsForPvivax <- resultsIn[,c(3,4,5,6,8,1,13,9,11)]
  merge3<- resultsIn[,c(3,4,5,6,7,9,1,15,10,12,13)]
  print("check making pvivax table after selection")
  print(head(merge3))
  onlyVivax <- merge3[grep("vivax", merge3$Malaria_Species),c (1,3,4,5,6,7,11,10)]
  #make a basic table for vivax samples
  if(length(rownames(onlyVivax)) > 0){

  rownames(onlyVivax) <- NULL

  vivaxOnly_kable  <- knitr::kable(onlyVivax, "html", padding = 40, line_sep = 2, align = "c", caption = "Table 11: <i>P. vivax</i> samples processed by CDC") %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12) %>%
          add_header_above(c("CDC_IDs" = 1,"Sample Metadata" = 4, "Genotyping Results"= 3 ),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray",  background = "#CDCDC1", color = "black", extra_css = "line-height: 20pt;") %>%
          add_header_above(c("Samples Positive for P. vivax by PET-PCR"=8), italic = FALSE, font_size = 16, align = "center", line = FALSE, color = "black", background = "#FFD39B", extra_css = "line-height: 20pt;") %>%
          column_spec(c(1), border_right = "2px solid gray") %>%
          column_spec(c(5), border_right = "2px solid gray") 
  } else {

    onlyVivax<- data.frame(matrix(ncol = 8, nrow = 1))
    onlyVivax[1,1] <- "No P. vivax Samples"
    onlyVivax[1,c(2:8)] <- "N/A"
    colnames(onlyVivax) <- c("LSDB_Specimen_ID", "State_Lab_ID", "Case_ID", "Collection_Date", "Reported_Travel_History", "Sample_Status", "Genetic_Cluster_Number", "GeoPrediction")
    vivaxOnly_kable  <- knitr::kable(onlyVivax, "html", padding = 40, line_sep = 2, align = "c", caption = "Table 11: <i>P. vivax</i> samples processed by CDC") %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12) %>%
        add_header_above(c("No Samples Positive for P. vivax by PET-PCR"=8),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray",  background = "#CDCDC1", color = "black", extra_css = "line-height: 20pt;")
    

  }
    vivaxOnly_kable <- column_spec(vivaxOnly_kable, 1, width = "7cm")
    vivaxOnly_kable <- column_spec(vivaxOnly_kable, 2, width = "7cm")
    vivaxOnly_kable <- column_spec(vivaxOnly_kable, 3, width = "7cm")
    vivaxOnly_kable <- column_spec(vivaxOnly_kable, 4, width = "7cm")
    vivaxOnly_kable <- column_spec(vivaxOnly_kable, 5, width = "7cm")
    vivaxOnly_kable <- column_spec(vivaxOnly_kable, 6, width = "7cm")
    vivaxOnly_kable <- column_spec(vivaxOnly_kable, 7, width = "7cm")
    vivaxOnly_kable <- column_spec(vivaxOnly_kable, 8, width = "7cm")
  
  # colnames
  # merge3 <- merge(resultsForPvivax, parnasClusters, by = "LSDB_Sequence_ID", all.x = T)
  print("before complete cases")
  print(merge3$GeneticClusterNumber)
  # print(merge3)
  merge3 <- merge3[complete.cases(merge3[,'GeneticClusterNumber']),]
  merge3 <- merge3[!grepl("Not_Tested", merge3$GeneticClusterNumber),]
  merge3 <- merge3[!grepl("Incomplete", merge3$GeneticClusterNumber),]
  
  print("test merg3")
  print(merge3)
  if(length(rownames(merge3)) > 0){
    print("has cluster")
    #This loop is for states that have a pvivax genotype

    #use the tree from each state (if samples have been genotyped from that state)'
    # setwd("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/Plasmodium_reporting_folder/tree_outputs")
    # treeList <- list.files(pattern = paste("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/Plasmodium_reporting_folder/tree_outputs/",stateID, "_Pvivax_ampliseq_tree.pdf", sep = ""))
    # treeFile <- list.files(pattern = paste(stateID, "_Pvivax_ampliseq_tree.pdf", sep = ""))
    #print(treeList)
    # treeFile <- tail(treeList, n = 1)
    # if(length(treeList) > 0){ 
    #   treeFile <- tail(treeList, n = 1)
    # }else{
    #   treeList <- list.files(pattern = paste("allStates_Pvivax_ampliseq_tree.pdf", sep = ""))
    #   treeFile <- tail(treeList, n = 1)
    # }
    # print(treeFile)
    # setwd("..")
    # treeFile <- paste("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/Plasmodium_reporting_folder/tree_outputs/",treeFile, sep ="")
    # print(treeFile)
    # treeFile = "/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/Plasmodium_reporting_folder/tree_outputs/2024-12-31_1038_allStatesBlank_Pvivax_ampliseq_tree.pdf"
    treeFile <- list.files(pattern = paste(stateID, "_Pvivax_ampliseq_tree.pdf", sep = ""))

    # treeList <- list.files(treeFolder,full.names = T)
    # matching_files <- grepl(stateID, treeList)
    # if(length(matching_files)>0){
    #   treeFile <- tail(matching_files, n = 1)
    # }else{
    #   treeList <- list.files(pattern = paste(stateID, "_Pvivax_ampliseq_tree.pdf", sep = ""))
    #   treeFile <- tail(treeList, n = 1)
    # }
    # treeFile <- "2024-02-20_1501_Pvivax_allStates_tree.pdf" 

    #Make the rownames a new column so we can use the pack_rows command
    merge3 <- merge3[order(merge3$GeneticClusterNumber), ]
    row.names(merge3) <- 1:nrow(merge3)
    merge3$name_of_rows <- row.names(merge3)

    # print(colnames(merge3))
    #Make a DF that doesn't have the name of row columns (looks cleaner in the output)
    testNoNames <- merge3[,-12]
  
    print(colnames(testNoNames))
    print("before select colnames for rename and table")
    testNoNames <- testNoNames[,c(1,3,4,5,6,11,10)]
    colnames(testNoNames)[1] <- "LSDB_Specimen_ID"
    colnames(testNoNames)[2] <- "State_Lab_ID"
    colnames(testNoNames)[3] <- "Case_ID"
    colnames(testNoNames)[4] <- "Collection_date"
    colnames(testNoNames)[5] <- "Reported_Travel_History"
    colnames(testNoNames)[6] <- "Genetic_Cluster_Number"
    colnames(testNoNames)[7] <- "GeoPrediction"
    
    #make the kable table that will be in the rmarkdown report
    kable_out <- knitr::kable(testNoNames, "html", padding = 40, line_sep = 2, align = "c", caption = "Table 12: Genetic Clustering and Geographic Prediction Results for Samples Analyzed in CDC's <i>P. vivax</i> AmpliSeq Assay") %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center", font_size = 12) %>%
          add_header_above(c("CDC_IDs" = 1,"Sample Metadata" = 4, "Genotyping Results"= 2 ),italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray",  background = "#CDCDC1", color = "black", extra_css = "line-height: 20pt;") %>%
          add_header_above(c("Genetic Clusters Detected"= 7), italic = FALSE, font_size = 16, align = "center", line = FALSE, color = "black", background = "#FFEC8B", extra_css = "line-height: 20pt;") %>%
          column_spec(c(1), border_right = "2px solid gray") %>%
          column_spec(c(5), border_right = "2px solid gray") 
    #loop throug the genetic clusters to find all specimens in each genetic cluster
    #Use the merge3 DF that still has the name of rows column
    list_of_genetic_clusters <- as.numeric(unique(merge3$GeneticClusterNumber))
    print(testNoNames)
    #print(l)
    for(l in list_of_genetic_clusters) {
      # print(l)

      filtered_to_current_cluster <- filter(merge3, GeneticClusterNumber == l)
      # print(filtered_to_current_cluster)
      specimens_to_look_at <- as.numeric(filtered_to_current_cluster$name_of_rows)

      #This is where the individual row numbers are extracted and then the range is used to pack rows of the same cluster into the same subtable
      left_value <- min(specimens_to_look_at)
      right_value <- max(specimens_to_look_at)
      
      kable_out <- kable_out %>%
        pack_rows("Genetic cluster detected:", left_value, right_value, bold = TRUE, indent = FALSE,  label_row_css = "border-top: 2px solid; border-bottom: 2px solid; color:#704EA5")	
      
    }
    #Make column widths uniform (or can change if think it looks better)
    kable_out <- column_spec(kable_out, 1, width = "7cm")
    kable_out <- column_spec(kable_out, 2, width = "7cm")
    kable_out <- column_spec(kable_out, 3, width = "7cm")
    kable_out <- column_spec(kable_out, 4, width = "7cm")
    kable_out <- column_spec(kable_out, 5, width = "7cm")
    kable_out <- column_spec(kable_out, 6, width = "7cm")
    kable_out <- column_spec(kable_out, 7, width = "7cm")
    # kable_out <- column_spec(kable_out, 8, width = "7cm")
    # kable_out <- column_spec(kable_out, 9, width = "7cm")
  
  } else{
    print("has no cluster")
    #This loop is for the states that do not have a Pvivax genotype

    #Since no samples have been genotyped, use a tree that has no specimens labeled
    #Will need to update this once more of the code is finished
    # treeFile <- "2024-02-20_1501_Pvivax_allStates_tree.pdf"
    # setwd("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/Plasmodium_reporting_folder/tree_outputs")
    treeList <- list.files(pattern = paste(stateID, "_allStatesBlank_Pvivax_ampliseq_tree.pdf", sep = ""))
    treeFile <- tail(treeList, n = 1)
    # setwd("..")
    # treeFile <- paste("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/Plasmodium_reporting_folder/tree_outputs/",treeFile, sep ="")
    #update this once I figure out how to do all states into one report, use blank tree for states with no specimens
    # treeFile = "/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/Plasmodium_reporting_folder/tree_outputs/2024-12-31_1038_allStatesBlank_Pvivax_ampliseq_tree.pdf"

    testNoNames <- data.frame(matrix(ncol = 7, nrow = 1))
    testNoNames[1,1] <- "No Samples Tested"
    testNoNames[1,c(2:7)] <- "N/A"
		# colnames(merge3) <- colnames(resultsForPvivax)
    # rownames(merge3) <- 1:nrow(merge3)
    # merge3$name_of_rows <- rownames(merge3)
    # print(colnames(merge3))
    # testNoNames <- merge3[,-11]
    # # testNoNames <- testNoNames[,c(1,6,7,2,3,8,9)]
    # testNoNames <- testNoNames[,c(1,7,2,3,9,8)]
    #print(colnames(testNoNames))
    #print("before select colnames for rename and table")
    colnames(testNoNames)[1] <- "LSDB_Specimen_ID"
    colnames(testNoNames)[2] <- "State_Lab_ID"
    colnames(testNoNames)[3] <- "Case_ID"
    colnames(testNoNames)[4] <- "Collection_date"
    colnames(testNoNames)[5] <- "Reported_Travel_History"
    colnames(testNoNames)[6] <- "Genetic_Cluster_Number"
    colnames(testNoNames)[7] <- "GeoPrediction"
    kable_out <- knitr::kable(testNoNames, "html", padding = 40, line_sep = 2, align = "c", caption = "Table 12: Genetic Clustering and Geographic Prediction Results for Samples Analyzed in CDC's <i>P. vivax</i> AmpliSeq Assay") %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center", font_size = 12)
    #Change the =9 if you update the number of columns in the final merge2 
    kable_out <- add_header_above(kable_out, c("No P. vivax Samples with Sufficient Data for Genetic Clustering"= 7), italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", background = "#FF7F00", color = "black", extra_css = "line-height: 20pt;")
           
    #Make column widths uniform (or can change if think it looks better)
    kable_out <- column_spec(kable_out, 1, width = "7cm")
    kable_out <- column_spec(kable_out, 2, width = "7cm")
    kable_out <- column_spec(kable_out, 3, width = "7cm")
    kable_out <- column_spec(kable_out, 4, width = "7cm")
    kable_out <- column_spec(kable_out, 5, width = "7cm")
    kable_out <- column_spec(kable_out, 6, width = "7cm")
    kable_out <- column_spec(kable_out, 7, width = "7cm")
    # kable_out <- column_spec(kable_out, 7, width = "7cm")
    # kable_out <- column_spec(kable_out, 8, width = "7cm")
    # kable_out <- column_spec(kable_out, 9, width = "7cm")
    }

  #Now we're going to write a couple of different tables to csv files, so that they can be provided as downloadable links in the html.
  
  #First is the pvivax genetic cluster erport
  #This csv file will then be embedded in the rmarkdown report
  # csvName <- paste("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/Plasmodium_reporting_folder/reports/csv_outputs/",date_now, "_", time_now,"_", stateID, "_pvivax_clustering_report.csv", sep ="")
  csvName <- paste(fullPath,date_now, "_", time_now,"_", stateID, "_pvivax_clustering_report.csv", sep ="")
  
  write.csv(onlyVivax, file = csvName, quote = F, row.names =F)
  # write.csv(testNoNames, file = csvName, quote = F, row.names =F)

  #Next write the summary table of all submissions and the pfal MaRS snp results
  print(resultsForTable1_noGeo)
  # write.csv("test.csv", resultsForTable1_noGeo, quote =F, row.names = F)
  write.csv(resultsForTable1_noGeo, file = specimenProcessingCSV, quote = F, row.names =F)
  write.csv(completeSNPTable, file = completeMaRSName, quote = F, row.names =F)
  
  #make the rmarkdown report
  # knitr::knit_meta(class=NULL,clean=T)
  # rmarkdown::render("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/Plasmodium_reporting_folder/code_folder/DMS_report_v2.Rmd", "html_document",  run_pandoc = TRUE)
  # rmarkdown::render("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/Plasmodium_reporting_folder/code_folder/DMS_report_v2.Rmd", "html_document",  run_pandoc = TRUE)
  # rmarkdown::render("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/pvivax_nf_core/pvivax-pvivaxclustering/bin/DMS_report_v2.Rmd", "html_document",  run_pandoc = TRUE)
  rmarkdown::render(paste(fullPath,"DMS_report_portfolio_v2.Rmd", sep =""), "html_document",  run_pandoc = TRUE)

  htmlName <- paste(fullPath,date_now, "_", time_now,"_", stateID, "_plasmodium_genotyping_report.html", sep = "")
  # file.copy(from = "/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/pvivax_nf_core/pvivax-pvivaxclustering/bin/DMS_report_v2.html", to = paste(htmlName, sep = ""))
  file.copy(from = paste(fullPath,"DMS_report_portfolio_v2.html", sep  =""), to = paste(htmlName, sep = ""))

  #remove the tmp report files after the rmarkdown has been made, will make folders cleaner
  # file.remove("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/pvivax_nf_core/pvivax-pvivaxclustering/bin/DMS_report_v2.html")
  # file.remove(csvName)
  # file.remove(specimenProcessingCSV)
  # file.remove(completeMaRSName)

}