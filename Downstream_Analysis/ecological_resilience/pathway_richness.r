args <- commandArgs(TRUE)

inFile <- read.table(args[1], sep = "\t", header = T)

suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gridExtra))


name <- args[4]

##Pathway Richness - calculate the number of species contributing to each function. This doesn't include the unclassified taxa, which may ultimately have an impact
melt(inFile, id.vars =c("Pathway", "Species")) -> melted
for (i in 1:nrow(melted)){
    if(grepl("bftm", melted$variable[i]) == TRUE){
        melted$pop[i] <- "bf"
        } else if (grepl("Ita", melted$variable[i]) == TRUE){
            melted$pop[i] <- "italy"
            } else if (grepl("Had", melted$variable[i]) == TRUE){
                melted$pop[i] <- "hadza"
                } else if (grepl("Ind", melted$variable[i]) == TRUE){
                    melted$pop[i] <- "euro"
                } else if (grepl("SRR399", melted$variable[i]) == TRUE){
                    melted$pop[i] <- "mongolia"
                } else if (grepl("NO", melted$variable[i]) == TRUE){
                    melted$pop[i] <- "norman"
                } else if (grepl("SM", melted$variable[i]) == TRUE){
                    melted$pop[i] <- "matses"
                }
}
melted[grep("TRUE", melted$value > 0),] -> pathwayPresent
pathwayPresent %>% group_by(variable, pop) %>% summarize(count=n()) -> pathwayRich
tapply(pathwayRich$count, pathwayRich$pop, mean) -> pathMean
tapply(pathwayRich$count, pathwayRich$pop, median) -> pathMedian
t(as.data.frame(pathMean)) -> meanDF
as.data.frame(meanDF) -> meanDF
pathwayRich$factor <- as.factor(pathwayRich$pop)
kruskal.test(pathwayRich$count, pathwayRich$factor) -> sig
meanDF$pval <- sig$p.value
row.names(meanDF) <- args[3]
#File needs to be written to temporary file, to append each value from loop. Final matrix is produced and written to file in ecol_plot.r
write.table(meanDF, file = "/Users/dave/Desktop/reference_based_mapping/tempRich.txt", sep = "\t", quote = F, row.names =T, col.names = F, append = T)
print(pathwayRich$count)
#boxplot of pathway richness between populations
pathway <- paste(args[3])
richGG <- ggplot(data = pathwayRich, mapping = aes(x = pathwayRich$pop, y = pathwayRich$count))
richGG + geom_boxplot( lier.colour = NA, aes(colour = pathwayRich$pop)) + xlab("Population") + ylab("Number of Species Contributing") + ggtitle(paste0(pathway, " Pathway Richness")) + theme(plot.title = element_text(hjust = 0.5)) 
ggsave(paste0(pathway, "_richness.png"))



## Calculate Gini-Simpson and Bray-Curtis seperated out by population, additionally calculate hill numbers and Pileu's evenness
suppressMessages(library(vegan))

#Ita and Had samples come from Italy and Hadza, Schnorr 2015
grep("Ita", names(inFile)) -> itaIndex
inFile[,c(2,itaIndex)] -> itaFile
itaMatrix <- itaFile[,-1]
rownames(itaMatrix) <- itaFile[,1]
diversity(t(itaMatrix), index = "simpson") -> itaSimpson
1 - itaSimpson -> itaGiniSimpson
vegdist(t(itaMatrix), method = "bray") -> itaBray
renyi(t(itaMatrix)) -> itaRyn
exp(itaRyn$`0`) -> itaHill0
exp(itaRyn$`1`) -> itaHill1
exp(itaRyn$`2`) -> itaHill2
itaH <- diversity(t(itaMatrix))
itaJ <- itaH/log(specnumber(t(itaMatrix)))
as.matrix(itaJ) -> itaJMat
colnames(itaJMat) <- "evenness"


grep("Had", names(inFile)) -> hadIndex
inFile[,c(2,hadIndex)] -> hadFile
hadMatrix <- hadFile[,-1]
rownames(hadMatrix) <- hadFile[,1]
diversity(t(hadMatrix), index = "simpson") -> hadSimpson
1 - hadSimpson -> hadGiniSimpson
vegdist(t(hadMatrix), method = "bray") -> hadBray
renyi(t(hadMatrix)) -> hadRyn
exp(hadRyn$`0`) -> hadHill0
exp(hadRyn$`1`) -> hadHill1
exp(hadRyn$`2`) -> hadHill2
hadH <- diversity(t(hadMatrix))
hadJ <- hadH/log(specnumber(t(hadMatrix)))
as.matrix(hadJ) -> hadJMat
colnames(hadJMat) <- "evenness"

#bftm are from Burkina Faso, samples collected by Centre Muraz, sequenced LMAMR
grep("bftm", names(inFile)) -> bftmIndex
inFile[,c(2,bftmIndex)] -> bftmFile
bftmMatrix <- bftmFile[,-1]
rownames(bftmMatrix) <- bftmFile[,1]
diversity(t(bftmMatrix), index = "simpson") -> bftmSimpson
1 - bftmSimpson -> bftmGiniSimpson
vegdist(t(bftmMatrix), method = "bray") -> bftmBray
renyi(t(bftmMatrix)) -> bftmRyn
exp(bftmRyn$`0`) -> bftmHill0
exp(bftmRyn$`1`) -> bftmHill1
exp(bftmRyn$`2`) -> bftmHill2
bftmH <- diversity(t(bftmMatrix))
bftmJ <- bftmH/log(specnumber(t(bftmMatrix)))
as.matrix(bftmJ) -> bftmJMat
colnames(bftmJMat) <- "evenness"


#Ind samples are healthy European women, from Karlsson 2013
grep("Ind", names(inFile)) -> euroIndex
inFile[,c(2,euroIndex)] -> euroFile
euroMatrix <- euroFile[,-1]
rownames(euroMatrix) <- euroFile[,1]
diversity(t(euroMatrix), index = "simpson") -> euroSimpson
1 - euroSimpson -> euroGiniSimpson
vegdist(t(euroMatrix), method = "bray") -> euroBray
renyi(t(euroMatrix)) -> euroRyn
exp(euroRyn$`0`) -> euroHill0
exp(euroRyn$`1`) -> euroHill1
exp(euroRyn$`2`) -> euroHill2
euroH <- diversity(t(euroMatrix))
euroJ <- euroH/log(specnumber(t(euroMatrix)))
as.matrix(euroJ) -> euroJMat
colnames(euroJMat) <- "evenness"

#SRR399 samples are from Mongolians
grep("SRR399", names(inFile)) -> monIndex
inFile[,c(2,monIndex)] -> monFile
monMatrix <- monFile[,-1]
rownames(monMatrix) <- monFile[,1]
diversity(t(monMatrix), index = "simpson") -> monSimpson
1 - monSimpson -> monGiniSimpson
vegdist(t(monMatrix), method = "bray") -> monBray
renyi(t(monMatrix)) -> monRyn
exp(monRyn$`0`) -> monHill0
exp(monRyn$`1`) -> monHill1
exp(monRyn$`2`) -> monHill2
monH <- diversity(t(monMatrix))
monJ <- monH/log(specnumber(t(monMatrix)))
as.matrix(monJ) -> monJMat
colnames(monJMat) <- "evenness"

##NO and SM samples are from Norman, OK and Matses HGs from Peru, respectively, sequenced LMAMR
grep("NO", names(inFile)) -> normanIndex
inFile[,c(2,normanIndex)] -> normanFile
normanMatrix <- normanFile[,-1]
rownames(normanMatrix) <- normanFile[,1]
diversity(t(normanMatrix), index = "simpson") -> normanSimpson
1 - normanSimpson -> normanGiniSimpson
vegdist(t(normanMatrix), method = "bray") -> normanBray
renyi(t(normanMatrix)) -> normanRyn
exp(normanRyn$`0`) -> normanHill0
exp(normanRyn$`1`) -> normanHill1
exp(normanRyn$`2`) -> normanHill2
normanH <- diversity(t(normanMatrix))
normanJ <- normanH/log(specnumber(t(normanMatrix)))
as.matrix(normanJ) -> normanJMat
colnames(normanJMat) <- "evenness"

grep("SM", names(inFile)) -> matsesIndex
inFile[,c(2,matsesIndex)] -> matsesFile
matsesMatrix <- matsesFile[,-1]
rownames(matsesMatrix) <- matsesFile[,1]
diversity(t(matsesMatrix), index = "simpson") -> matsesSimpson
1 - matsesSimpson -> matsesGiniSimpson
vegdist(t(matsesMatrix), method = "bray") -> matsesBray
renyi(t(matsesMatrix)) -> matsesRyn
exp(matsesRyn$`0`) -> matsesHill0
exp(matsesRyn$`1`) -> matsesHill1
exp(matsesRyn$`2`) -> matsesHill2
matsesH <- diversity(t(matsesMatrix))
matsesJ <- matsesH/log(specnumber(t(matsesMatrix)))
as.matrix(matsesJ) -> matsesJMat
colnames(matsesJMat) <- "evenness"

rbind(itaJMat, matsesJMat, hadJMat, normanJMat, bftmJMat, monJMat, euroJMat) -> fullJMat

as.data.frame(fullJMat) -> fullJDF
fullJDF$sample <- row.names(fullJDF)

for (i in 1:nrow(fullJDF)){
    if(grepl("bftm", fullJDF$sample[i]) == TRUE){
        fullJDF$pop[i] <- "bf"
    } else if (grepl("Ita", fullJDF$sample[i]) == TRUE){
        fullJDF$pop[i] <- "italy"
    } else if (grepl("Had", fullJDF$sample[i]) == TRUE){
        fullJDF$pop[i] <- "hadza"
    } else if (grepl("Ind", fullJDF$sample[i]) == TRUE){
       fullJDF$pop[i] <- "euro"
    } else if (grepl("SRR399", fullJDF$sample[i]) == TRUE){
        fullJDF$pop[i] <- "mongolia"
    } else if (grepl("NO", fullJDF$sample[i]) == TRUE){
        fullJDF$pop[i] <- "norman"
    } else if (grepl("SM", fullJDF$sample[i]) == TRUE){
        fullJDF$pop[i] <- "matses"
    }
}


evennessGG <- ggplot(fullJDF, mapping = aes(x= fullJDF$pop, y = fullJDF$evenness))

evennessGG + geom_boxplot()

ggsave("temp.pdf")

matrix(ncol = 3, nrow = 7) -> ecolMatrix
ecolMatrix[1,] <- c("Italy", round(mean(itaSimpson), 2), round(mean(itaBray), 2))
ecolMatrix[2,] <- c("Burkina_Faso", round(mean(bftmSimpson),2), round(mean(bftmBray),2))
ecolMatrix[3,] <- c("Hadza", round(mean(hadSimpson),2), round(mean(hadBray),2))
ecolMatrix[4,] <- c("Euro", round(mean(euroSimpson),2), round(mean(euroBray),2))
ecolMatrix[5,] <- c("Mongolian", round(mean(monSimpson),2), round(mean(monBray),2))
ecolMatrix[6,] <- c("Norman", round(mean(normanSimpson),2), round(mean(normanBray),2))
ecolMatrix[7,] <- c("Matses", round(mean(matsesSimpson),2), round(mean(matsesBray),2))
row.names(ecolMatrix) <- c(args[3], args[3], args[3], args[3], args[3], args[3], args[3])
write.table(ecolMatrix, "/Users/dave/Desktop/reference_based_mapping/delete.txt", append = T, col.names = F, row.names = T, sep  = "\t")

as.data.frame(ecolMatrix) -> ecolDF
colnames(ecolDF) <- c("Pop", "Mean_Gini_Simpson", "Mean_Bray_Curtis")
ecolDF$pathway <- c(args[3], args[3], args[3], args[3], args[3], args[3], args[3])
ecolGG <- ggplot(data = ecolDF, mapping = aes(x = as.numeric(as.character(ecolDF$Mean_Gini_Simpson)), y = as.numeric(as.character(ecolDF$Mean_Bray_Curtis)), shape = ecolDF$pathway))
ecolGG + geom_point(aes(colour = ecolDF$Pop), size = 4) + ylim(0,1) + xlim(0,1) + xlab("Gini-Simpson") + scale_shape_manual(values = c(19))+ ylab("Bray-Curtis") + scale_colour_manual(values = c(Burkina_Faso ="darkblue", Euro = "#fdbf6f", Hadza = "darkviolet", Italy = "#ff7f00", Mongolian = "gray54", Norman = "#e31a1c", Matses = "cornflowerblue")) + theme(panel.background = element_rect(fill = "whitesmoke"), panel.grid.minor = element_blank()) + ggtitle(pathway, name)

ggsave(paste0(pathway, "conDiv.pdf"), width = 10, height = 6)

merge(itaMatrix, hadMatrix, by = 0, all = T) -> temp
temp[,-1] -> firstMerge
temp[,1] -> rownames(firstMerge)
merge(firstMerge, bftmMatrix, by =0, all = T) -> temp2
temp2[,-1] -> secondMerge
temp2[,1] -> rownames(secondMerge)
merge(secondMerge, monMatrix, by =0, all = T) -> temp3
temp3[,-1] -> thirdMerge
temp3[,1] -> rownames(thirdMerge)
merge(thirdMerge, euroMatrix, by = 0, all =T) -> temp4
temp4[,-1] -> fourthMerge
temp4[,1] -> rownames(fourthMerge)
merge(fourthMerge, normanMatrix, by = 0, all =T) -> temp5
temp5[,-1] -> fifthMerge
temp5[,1] -> rownames(fifthMerge)
merge(fifthMerge, matsesMatrix, by = 0, all =T) -> temp6
temp6[,-1] -> sixthMerge
temp6[,1] -> rownames(sixthMerge)

write.table(sixthMerge, file = paste0(pathway, "com_matrix.txt"),  append =F, row.names =T, sep = "\t", quote = F)

max(pathwayRich$count) -> speciesMax
cbind(euroHill0, euroHill1, euroHill2) -> euroBind
row.names(euroBind) <- colnames(euroMatrix)
melt(euroBind) -> euroHillMelt
length(row.names(euroHillMelt)) -> euroNum
for( i in 1:euroNum){
    if(grepl("0", euroHillMelt$Var2[i]) == TRUE){
        euroHillMelt$hillNum[i] <- 0
    } else if(grepl("1", euroHillMelt$Var2[i]) == TRUE){
        euroHillMelt$hillNum[i] <- 1
    } else if(grepl("2", euroHillMelt$Var2[i]) == TRUE){
        euroHillMelt$hillNum[i] <- 2
    }
}
as.data.frame(euroHillMelt) -> euroHillDF
colnames(euroHillDF) <- c("SampleID", "hillVar", "value", "order_q")

p1 <- ggplot(euroHillMelt, aes(x = hillNum, y = value)) + geom_smooth(method = "auto", colour = "#fdbf6f") + ggtitle("European",paste0(pathway,"_",name)) + ylim(c(0, speciesMax)) + theme(panel.background = element_rect(fill = "whitesmoke"), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(0,speciesMax,10), limits = c(0, speciesMax)) + scale_x_continuous(breaks = seq(0,2,1))
#ggsave(paste0(pathway, "_euroHill.pdf"), width = 10, height = 6)

cbind(hadHill0, hadHill1, hadHill2) -> hadBind
row.names(hadBind) <- colnames(hadMatrix)
melt(hadBind) -> hadHillMelt
length(row.names(hadHillMelt)) -> hadNum
for( i in 1:hadNum){
    if(grepl("0", hadHillMelt$Var2[i]) == TRUE){
        hadHillMelt$hillNum[i] <- 0
    } else if(grepl("1", hadHillMelt$Var2[i]) == TRUE){
        hadHillMelt$hillNum[i] <- 1
    } else if(grepl("2", hadHillMelt$Var2[i]) == TRUE){
        hadHillMelt$hillNum[i] <- 2
    }
}
as.data.frame(hadHillMelt) -> hadHillDF
colnames(hadHillDF) <- c("SampleID", "hillVar", "value", "order_q")

p2 <- ggplot(hadHillMelt, aes(x = hillNum, y = value)) + geom_smooth(method = "auto", colour = "darkviolet") + ggtitle("Hadza",paste0(pathway,"_",name)) +  ylim(c(0, speciesMax))  + theme(panel.background = element_rect(fill = "whitesmoke"), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(0,speciesMax,10), limits = c(0, speciesMax)) + scale_x_continuous(breaks = seq(0,2,1))
#ggsave(paste0(pathway, "_hadzaHill.pdf"), width = 10, height = 6)

cbind(itaHill0, itaHill1, itaHill2) -> itaBind
row.names(itaBind) <- colnames(itaMatrix)
melt(itaBind) -> itaHillMelt
length(row.names(itaHillMelt)) -> itaNum
for( i in 1:itaNum){
    if(grepl("0", itaHillMelt$Var2[i]) == TRUE){
        itaHillMelt$hillNum[i] <- 0
    } else if(grepl("1", itaHillMelt$Var2[i]) == TRUE){
        itaHillMelt$hillNum[i] <- 1
    } else if(grepl("2", itaHillMelt$Var2[i]) == TRUE){
        itaHillMelt$hillNum[i] <- 2
    }
}
as.data.frame(itaHillMelt) -> itaHillDF
colnames(itaHillDF) <- c("SampleID", "hillVar", "value", "order_q")
p3 <- ggplot(itaHillMelt, aes(x = hillNum, y = value)) + geom_smooth(method = "auto", colour = "#ff7f00")  + ggtitle("Italy",paste0(pathway,"_",name))  + ylim(c(0, speciesMax))  + theme(panel.background = element_rect(fill = "whitesmoke"), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(0,speciesMax,10), limits = c(0, speciesMax)) + scale_x_continuous(breaks = seq(0,2,1))
#ggsave(paste0(pathway, "_italyHill.pdf"), width = 10, height = 6)

cbind(matsesHill0, matsesHill1, matsesHill2) -> matsesBind
row.names(matsesBind) <- colnames(matsesMatrix)
melt(matsesBind) -> matsesHillMelt
length(row.names(matsesHillMelt)) -> matsesNum
for( i in 1:matsesNum){
    if(grepl("0", matsesHillMelt$Var2[i]) == TRUE){
        matsesHillMelt$hillNum[i] <- 0
    } else if(grepl("1", matsesHillMelt$Var2[i]) == TRUE){
        matsesHillMelt$hillNum[i] <- 1
    } else if(grepl("2", matsesHillMelt$Var2[i]) == TRUE){
        matsesHillMelt$hillNum[i] <- 2
    }
}
as.data.frame(matsesHillMelt) -> matsesHillDF
colnames(matsesHillDF) <- c("SampleID", "hillVar", "value", "order_q")
p4 <- ggplot(matsesHillMelt, aes(x = hillNum, y = value)) + geom_smooth(method = "auto", colour = "cornflowerblue")  + ggtitle("Matses",paste0(pathway,"_",name)) + ylim(c(0, speciesMax))  + theme(panel.background = element_rect(fill = "whitesmoke"), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(0,speciesMax,10), limits = c(0, speciesMax)) + scale_x_continuous(breaks = seq(0,2,1))
#ggsave(paste0(pathway, "_matsesHill.pdf"), width = 10, height = 6)

cbind(normanHill0, normanHill1, normanHill2) -> normanBind
row.names(normanBind) <- colnames(normanMatrix)
melt(normanBind) -> normanHillMelt
length(row.names(normanHillMelt)) -> normanNum
for( i in 1:normanNum){
    if(grepl("0", normanHillMelt$Var2[i]) == TRUE){
        normanHillMelt$hillNum[i] <- 0
    } else if(grepl("1", normanHillMelt$Var2[i]) == TRUE){
        normanHillMelt$hillNum[i] <- 1
    } else if(grepl("2", normanHillMelt$Var2[i]) == TRUE){
        normanHillMelt$hillNum[i] <- 2
    }
}
as.data.frame(normanHillMelt) -> normanHillDF
colnames(normanHillDF) <- c("SampleID", "hillVar", "value", "order_q")
p5 <- ggplot(normanHillMelt, aes(x = hillNum, y = value)) + geom_smooth(method = "auto", colour = "#e31a1c")  + ggtitle("Norman",paste0(pathway,"_",name)) + ylim(c(0, speciesMax))  + theme(panel.background = element_rect(fill = "whitesmoke"), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(0,speciesMax,10), limits = c(0, speciesMax)) + scale_x_continuous(breaks = seq(0,2,1))
#ggsave(paste0(pathway, "_normanHill.pdf"), width = 10, height = 6)

cbind(bftmHill0, bftmHill1, bftmHill2) -> bftmBind
row.names(bftmBind) <- colnames(bftmMatrix)
melt(bftmBind) -> bftmHillMelt
length(row.names(bftmHillMelt)) -> bftmNum
for( i in 1:bftmNum){
    if(grepl("0", bftmHillMelt$Var2[i]) == TRUE){
        bftmHillMelt$hillNum[i] <- 0
    } else if(grepl("1", bftmHillMelt$Var2[i]) == TRUE){
        bftmHillMelt$hillNum[i] <- 1
    } else if(grepl("2", bftmHillMelt$Var2[i]) == TRUE){
        bftmHillMelt$hillNum[i] <- 2
    }
}
as.data.frame(bftmHillMelt) -> bftmHillDF
colnames(bftmHillDF) <- c("SampleID", "hillVar", "value", "order_q")
p6 <- ggplot(bftmHillMelt, aes(x = hillNum, y = value)) + geom_smooth(method = "auto", colour = "darkblue")   + ggtitle("Burkina Faso",paste0(pathway,"_",name)) + ylim(c(0, speciesMax))  + theme(panel.background = element_rect(fill = "whitesmoke"), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(0,speciesMax,10), limits = c(0, speciesMax)) + scale_x_continuous(breaks = seq(0,2,1))
#ggsave(paste0(pathway, "_bftmHill.pdf"), width = 10, height = 6)

cbind(monHill0, monHill1, monHill2) -> monBind
row.names(monBind) <- colnames(monMatrix)
melt(monBind) -> monHillMelt
length(row.names(monHillMelt)) -> monNum
for( i in 1:monNum){
    if(grepl("0", monHillMelt$Var2[i]) == TRUE){
        monHillMelt$hillNum[i] <- 0
    } else if(grepl("1", monHillMelt$Var2[i]) == TRUE){
        monHillMelt$hillNum[i] <- 1
    } else if(grepl("2", monHillMelt$Var2[i]) == TRUE){
        monHillMelt$hillNum[i] <- 2
    }
}
as.data.frame(monHillMelt) -> monHillDF
colnames(monHillDF) <- c("SampleID", "hillVar", "value", "order_q")
p7 <- ggplot(monHillMelt, aes(x = hillNum, y = value)) + geom_smooth(method = "auto", colour = "gray54")  + ggtitle("Mongolia",paste0(pathway,"_",name)) + ylim(c(0, speciesMax))  + theme(panel.background = element_rect(fill = "whitesmoke"), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(0,speciesMax,10), limits = c(0, speciesMax)) + scale_x_continuous(breaks = seq(0,2,1))
#ggsave(paste0(pathway, "_monHill.pdf"), width = 10, height = 6)

grid.arrange(p1, p3, p5, p7, p2, p4, p6, nrow =2)

g <- arrangeGrob(p1, p3, p5, p7, p2, p4, p6, nrow =2)

ggsave(file = paste0(pathway, "_hillNumbers.pdf"), g, width = 10, height =6)

rbind(monHillDF, bftmHillDF, normanHillDF, matsesHillDF, itaHillDF, hadHillDF, euroHillDF) -> popsHillBind

for (i in 1:nrow(popsHillBind)){
    if(grepl("bftm", popsHillBind$SampleID[i]) == TRUE){
        popsHillBind$pop[i] <- "Burkina_Faso"
    } else if (grepl("Ita", popsHillBind$SampleID[i]) == TRUE){
        popsHillBind$pop[i] <- "Italy"
    } else if (grepl("Had", popsHillBind$SampleID[i]) == TRUE){
        popsHillBind$pop[i] <- "Hadza"
    } else if (grepl("Ind", popsHillBind$SampleID[i]) == TRUE){
        popsHillBind$pop[i] <- "European"
    } else if (grepl("SRR399", popsHillBind$SampleID[i]) == TRUE){
        popsHillBind$pop[i] <- "Mongolia"
    } else if (grepl("NO", popsHillBind$SampleID[i]) == TRUE){
        popsHillBind$pop[i] <- "Norman"
    } else if (grepl("SM", popsHillBind$SampleID[i]) == TRUE){
        popsHillBind$pop[i] <- "Matses"
    }
}

popHillGG <- ggplot(popsHillBind, mapping = aes(x = popsHillBind$order_q, y = popsHillBind$value, colour = popsHillBind$pop))

popHillGG + geom_smooth(method = "auto", se = FALSE) + scale_colour_manual(values = c(Burkina_Faso ="darkblue", European = "#fdbf6f", Hadza = "darkviolet", Italy = "#ff7f00", Mongolia = "gray54", Norman = "#e31a1c", Matses = "cornflowerblue")) + theme(panel.background = element_rect(fill = "whitesmoke"), panel.grid.minor = element_blank()) + ggtitle(pathway, name)

ggsave(file = paste0(pathway, "_hillNumbers_1plot.pdf"), width = 10, height =6)
