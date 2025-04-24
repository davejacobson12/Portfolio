#!/bin/sh
args <- commandArgs(TRUE)


library(vegan)
library(ggplot2)
library(reshape2)

#inputMat must be taxa summarize file, with # removed and taxa as rows and samples as columns. Values must be proportions (value 0 to 1).

inputMat <- read.table(args[1], sep = "\t", header =T, row.names = 1)
inputMeta <- read.table(args[2], sep = "\t", header =T)

sampleCat <- args[3]

#calculate renyi values and Hill numbers from they renyi

exp(renyi(t(inputMat))) -> inputHill
inputHill$`0` -> hillZero
inputHill$`0.25` -> hillQuarter
inputHill$`0.5` -> hillHalf
inputHill$`1` -> hillOne
inputHill$`2` -> hillTwo

cbind(hillZero, hillQuarter, hillHalf, hillOne, hillTwo) -> hillBind
row.names(hillBind) <- colnames(inputMat)
melt(hillBind) -> hillMelt
length(row.names(hillMelt)) -> hillLength

for(i in 1:hillLength){
    if(grepl("Zero", hillMelt$Var2[i]) == TRUE){
        hillMelt$order_q[i] <- 0
    } else if(grepl("Quarter", hillMelt$Var2[i]) == TRUE){
        hillMelt$order_q[i] <- 0.25
    } else if(grepl("Half", hillMelt$Var2[i]) == TRUE){
        hillMelt$order_q[i] <- 0.5
    } else if(grepl("One", hillMelt$Var2[i]) == TRUE){
        hillMelt$order_q[i] <- 1
    } else if(grepl("Two", hillMelt$Var2[i]) == TRUE){
        hillMelt$order_q[i] <- 2
    }
}

as.data.frame(hillMelt) -> hillDF
colnames(hillDF) <- c("SampleID", "hillVar", "value", "order_q")

merge(x = hillDF, y = inputMeta) -> metaMerge

#print(metaMerge)

sampleCatIndex <- grep(args[3], colnames(metaMerge))

hillPlot <- ggplot(metaMerge, mapping = aes(x = metaMerge$order_q, y = metaMerge$value, colour = metaMerge[,sampleCatIndex])) + geom_smooth(method = "auto")

hillPlot

ggsave(paste0(sampleCat, "hillNum.pdf"))
