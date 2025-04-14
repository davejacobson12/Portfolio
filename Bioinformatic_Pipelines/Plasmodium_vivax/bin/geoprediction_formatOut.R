args <- commandArgs(trailingOnly = T)

#Create the merge metadata prediction file (Even though no meta is included, just a naming hold over)
gatk_country_predictionDF <- read.table(args[2], sep = "\t", header =T)
gatk_region_predictionDF <- read.table(args[3], sep = "\t", header =T)
gatk_barcodeCSV <- read.table(args[4], sep =",", header = T)

timestampIn <- strsplit(args[4],"_pvivax")[[1]][1]
print(timestampIn)

#reformat and rename columns for country and region prediction
gatk_country_predictionDF$Country_Prediction_1<- paste(gatk_country_predictionDF$Prediction_1, " (", round(gatk_country_predictionDF$Prob_1, digits = 5),")", sep = "")
gatk_country_predictionDF$Country_Prediction_2 <- paste(gatk_country_predictionDF$Prediction_2, " (", round(gatk_country_predictionDF$Prob_2, digits = 5),")", sep = "")
gatk_country_predictionDF$Country_Prediction_3 <- paste(gatk_country_predictionDF$Prediction_3, " (", round(gatk_country_predictionDF$Prob_3, digits = 5),")", sep = "")

gatk_region_predictionDF$Region_Prediction_1<- paste(gatk_region_predictionDF$Prediction_1, " (", round(gatk_region_predictionDF$Prob_1, digits = 5),")", sep = "")
gatk_region_predictionDF$Region_Prediction_2 <- paste(gatk_region_predictionDF$Prediction_2, " (", round(gatk_region_predictionDF$Prob_2, digits = 5),")", sep = "")
gatk_region_predictionDF$Region_Prediction_3 <- paste(gatk_region_predictionDF$Prediction_3, " (", round(gatk_region_predictionDF$Prob_3, digits = 5),")", sep = "")

#get the sample order (remove extra characters from name)
sampleOrder <- gatk_barcodeCSV$Sample
newVec <- c()
for(i in sampleOrder){
    newVar <- strsplit(i, split = "_")[[1]][1]
    tester <- substr(i, 1,15)
    checkChar <- substr(tester,15,15)
    if(checkChar == "_"){
        tester <- substr(tester,1,14)
    }else{
        tester <- tester
    }
    newVec <- c(newVec, tester)

}

#add sample names to file and write
# newVec2 <- gsub("_","-", newVec)
newVec2 <- gsub("_T1","", newVec)
gatk_country_predictionDF$Sample <- newVec2
gatk_region_predictionDF$Sample <- newVec2

gatk_country_predictionDF2 <- gatk_country_predictionDF[,c(15,1,2,3,4,5,6,7,8,9,10,11,12,13,14)]
gatk_region_predictionDF2 <- gatk_region_predictionDF[,c(15,1,2,3,4,5,6,7,8,9,10,11,12,13,14)]

write.table(gatk_country_predictionDF2, paste(timestampIn,"Country_prediction_withMetadata.txt", sep =""), quote = F, sep = "\t", row.names = F)
write.table(gatk_region_predictionDF2, paste(timestampIn,"Region_prediction_withMetadata.txt", sep =""), quote = F, sep = "\t", row.names = F)

#Read in DB of existing geopredictions for Pvivax samples already run
existingPredictions <- read.table(args[1], header = T, sep = "\t")

#format region predictions to get seq ID and top region prediction
newPredictions_clean <- gatk_region_predictionDF[,c(15,12)]
newPredictions_clean$Region_Prediction_1 <- gsub(" ", "_", newPredictions_clean$Region_Prediction_1)
newPredictions_clean$Sample <- gsub("-", "_", newPredictions_clean$Sample)
newPredictions_clean$Sample <- gsub("_T1", "", newPredictions_clean$Sample)
colnames(newPredictions_clean)[1] <- "Seq_ID"

#combine with existing predictions and write to table
combinedFile <- rbind(existingPredictions, newPredictions_clean)
combinedFile_final <- combinedFile[!duplicated(combinedFile$Seq_ID), ]

name_output <- paste(timestampIn,"_pvivax_geoPrediction_simple.txt", sep ="")
write.table(combinedFile_final, file = name_output, quote = F, sep  ="\t", row.names =F)