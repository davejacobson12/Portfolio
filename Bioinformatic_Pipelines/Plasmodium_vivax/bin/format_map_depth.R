

time_now <- format(Sys.time(), "%H%M")
date_now <- Sys.Date()
date_now <- gsub("-","", date_now)
timestamp <- paste0(date_now,"_",time_now)

myFiles <- list.files(pattern = "./*read_summary.txt")

#read in files for each sample and append to single df
fullDF <- data.frame(matrix(ncol=6,nrow = 0))
for(i in myFiles){
    inFile <- read.table(i, sep = "\t", header = F)
    inFile <- as.data.frame(inFile[,2])
    myOut <- t(inFile)
    fullDF <- rbind(fullDF, myOut)
}

#Format outputs and basic calculations
colnames(fullDF) <- c("LSDB_Sequence_ID", "Input_Reads", "AfterQCReads", "HumanMapped_Reads", "HumanUnmapped_Reads", "PercentRemainingReads_PvP01Mapped")

fullDF$AfterQCReads <- as.numeric(as.character(fullDF$AfterQCReads))
fullDF$Input_Reads <- as.numeric(as.character(fullDF$Input_Reads))
fullDF$HumanMapped_Reads<- as.numeric(as.character(fullDF$HumanMapped_Reads))
fullDF$HumanUnmapped_Reads <- as.numeric(as.character(fullDF$HumanUnmapped_Reads))
fullDF$PercentRemainingReads_PvP01Mapped <- gsub("%", "", fullDF$PercentRemainingReads_PvP01Mapped)
fullDF$ProportionRemainingReads_PvP01Mapped <- as.numeric(fullDF$PercentRemainingReads_PvP01Mapped) / 100

fullDF$Proportion_ReadsPassQC <- fullDF$AfterQCReads / fullDF$Input_Reads
fullDF$Proportion_CleanReads_Human <- fullDF$HumanMapped_Reads / fullDF$AfterQCReads
fullDF$Approx_PvP01_Mapped_Reads <- fullDF$ProportionRemainingReads_PvP01Mapped * fullDF$HumanUnmapped_Reads
fullDF$Proportion_CleanReads_PvP01 <- fullDF$Approx_PvP01_Mapped_Reads / fullDF$AfterQCReads

fullDF2 <- subset(fullDF, select = -c(PercentRemainingReads_PvP01Mapped))

write.csv(fullDF2, file = paste(timestamp, "_pvp01_swgs_mapStats.csv", sep = ""), row.names= F, quote =F)