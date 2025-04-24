
args <- commandArgs(TRUE)
library(ggplot2)

inFile <- read.table(args[1], header = T, sep = "\t")
pathway <- args[2]
name <- args[3]
inFile[,-1:-2] -> sampMatrix
length(inFile$Species) -> unclassIndex
unclassIndex - 1 -> speciesIndex
length(colnames(sampMatrix)) -> totalSamp
hadLen <- length(grep("TRUE",grepl("Had", colnames(inFile))))
bftmLen <- length(grep("TRUE",grepl("bftm", colnames(inFile))))
itaLen <- length(grep("TRUE",grepl("Ita", colnames(inFile))))
monLen <- length(grep("TRUE",grepl("SRR399", colnames(inFile))))
euroLen <- length(grep("TRUE",grepl("Ind", colnames(inFile))))
normanLen <- length(grep("TRUE",grepl("NO", colnames(inFile))))
matsesLen <- length(grep("TRUE",grepl("SM", colnames(inFile))))
fillMat <- matrix(ncol = 5, nrow = totalSamp )
for (i in 1:totalSamp){
    if(grepl("Had", colnames(sampMatrix)[i]) == "TRUE"){
        fillMat[i,1] <- colnames(sampMatrix)[i]
        fillMat[i,2] <- "hadza"
        fillMat[i,3] <- sum(sampMatrix[,i])
        fillMat[i,4] <- sampMatrix[unclassIndex,i]
        fillMat[i,5] <- (as.numeric(fillMat[i,4]) / as.numeric(fillMat[i,3]))
    } else if (grepl("bftm", colnames(sampMatrix)[i]) == "TRUE"){
        fillMat[i,1] <- colnames(sampMatrix)[i]
        fillMat[i,2] <- "bftm"
        fillMat[i,3] <- sum(sampMatrix[,i])
        fillMat[i,4] <- sampMatrix[unclassIndex,i]
        fillMat[i,5] <- (as.numeric(fillMat[i,4]) / as.numeric(fillMat[i,3]))
    } else if (grepl("Ind", colnames(sampMatrix)[i]) == "TRUE"){
        fillMat[i,1] <- colnames(sampMatrix)[i]
        fillMat[i,2] <- "euro"
        fillMat[i,3] <- sum(sampMatrix[,i])
        fillMat[i,4] <- sampMatrix[unclassIndex,i]
        fillMat[i,5] <- (as.numeric(fillMat[i,4]) / as.numeric(fillMat[i,3]))
    } else if (grepl("NO", colnames(sampMatrix)[i]) == "TRUE"){
        fillMat[i,1] <- colnames(sampMatrix)[i]
        fillMat[i,2] <- "norman"
        fillMat[i,3] <- sum(sampMatrix[,i])
        fillMat[i,4] <- sampMatrix[unclassIndex,i]
        fillMat[i,5] <- (as.numeric(fillMat[i,4]) / as.numeric(fillMat[i,3]))
    } else if (grepl("SM", colnames(sampMatrix)[i]) == "TRUE"){
        fillMat[i,1] <- colnames(sampMatrix)[i]
        fillMat[i,2] <- "matses"
        fillMat[i,3] <- sum(sampMatrix[,i])
        fillMat[i,4] <- sampMatrix[unclassIndex,i]
        fillMat[i,5] <- (as.numeric(fillMat[i,4]) / as.numeric(fillMat[i,3]))
    } else if (grepl("SRR399", colnames(sampMatrix)[i]) == "TRUE"){
        fillMat[i,1] <- colnames(sampMatrix)[i]
        fillMat[i,2] <- "mongolia"
        fillMat[i,3] <- sum(sampMatrix[,i])
        fillMat[i,4] <- sampMatrix[unclassIndex,i]
        fillMat[i,5] <- (as.numeric(fillMat[i,4]) / as.numeric(fillMat[i,3]))
    } else if (grepl("Ita", colnames(sampMatrix)[i]) == "TRUE"){
        fillMat[i,1] <- colnames(sampMatrix)[i]
        fillMat[i,2] <- "italy"
        fillMat[i,3] <- sum(sampMatrix[,i])
        fillMat[i,4] <- sampMatrix[unclassIndex,i]
        fillMat[i,5] <- (as.numeric(fillMat[i,4]) / as.numeric(fillMat[i,3]))
    }
}

colnames(fillMat) <- c("sample", "pop", "mapped_taxa", "mapped_unclass", "percent_unclass")
as.data.frame(fillMat) -> fillDF
as.numeric(as.character(fillDF$percent_unclass)) -> fillDF$percent_unclass

t(tapply(fillDF$percent_unclass, fillDF$pop, mean)) -> meanUnclass

row.names(meanUnclass) <- pathway

as.data.frame(meanUnclass) -> unclassDF

write.table(unclassDF, file = "percentUnclass.txt",  append =T, row.names =T, sep = "\t", quote = F, col.names =F)

ggUnclass <- ggplot(data = fillDF, mapping = aes(x = fillDF$pop, y = fillDF$percent_unclass)) + ggtitle(pathway, name)

ggUnclass + geom_boxplot()

ggsave(paste0(pathway, "percentUnclassified.png"), width = 10, height = 6)

