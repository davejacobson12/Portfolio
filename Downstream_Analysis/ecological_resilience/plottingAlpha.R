args <- commandArgs(TRUE)
taxaLevel <- args[4]


library(ggplot2)
library(gridExtra)

inputAce <- read.table(args[1], header = T, sep = "\t", quote = "", fill = FALSE)
inputBut <- read.table(args[2], header = T, sep = "\t", quote = "", fill = FALSE)
inputProp <- read.table(args[3], header = T, sep = "\t", quote = "", fill = FALSE)

inputAce$lifestyleSpecific <- factor(inputAce$lifestyleSpecific, levels = c("westernIndustrial", "nonWesternIndustrial", "pastoral", "ruralAgriculture", "hunterGather"))
inputBut$lifestyleSpecific <- factor(inputBut$lifestyleSpecific, levels = c("westernIndustrial", "nonWesternIndustrial", "pastoral", "ruralAgriculture", "hunterGather"))
inputProp$lifestyleSpecific <- factor(inputProp$lifestyleSpecific, levels = c("westernIndustrial", "nonWesternIndustrial", "pastoral", "ruralAgriculture", "hunterGather"))

pdf(paste(taxaLevel, "3_SCFA_richness.pdf", sep = "_"))
par(mfrow = c(2,2))
boxplot(inputAce$richness ~ inputAce$lifestyleSpecific, main = paste("Acetate", taxaLevel, "Richness", sep = " "), cex.axis = 0.65) -> aceRich_plot
boxplot(inputBut$richness ~ inputBut$lifestyleSpecific, main = paste("Butyrate", taxaLevel, "Richness", sep = " "), cex.axis = 0.65) -> butRich_plot
boxplot(inputProp$richness ~ inputProp$lifestyleSpecific, main = paste("Propionate", taxaLevel, "Richness", sep = " "), cex.axis = 0.65) -> propRich_plot
dev.off()

pdf(paste(taxaLevel, "3_SCFA_PD.pdf", sep = "_"))
par(mfrow = c(2,2))
boxplot(inputAce$PD ~ inputAce$lifestyleSpecific, main = paste("Acetate", taxaLevel, "PD", sep = " "), cex.axis = 0.65) -> acePD_plot
boxplot(inputBut$PD ~ inputBut$lifestyleSpecific, main = paste("Butyrate", taxaLevel, "PD", sep = " "), cex.axis = 0.65) -> butPD_plot
boxplot(inputProp$PD ~ inputProp$lifestyleSpecific, main = paste("Propionate", taxaLevel, "PD", sep = " "), cex.axis = 0.65) -> propPD_plot
dev.off()

pdf(paste(taxaLevel, "3_SCFA_GS.pdf", sep = "_"))
par(mfrow = c(2,2))
boxplot(inputAce$GS ~ inputAce$lifestyleSpecific, main = paste("Acetate", taxaLevel, "GS", sep = " "), cex.axis = 0.65) -> aceGS_plot
boxplot(inputBut$GS ~ inputBut$lifestyleSpecific, main = paste("Butyrate", taxaLevel, "GS", sep = " "), cex.axis = 0.65) -> butGS_plot
boxplot(inputProp$GS ~ inputProp$lifestyleSpecific, main = paste("Propionate", taxaLevel, "GS", sep = " "), cex.axis = 0.65) -> propGS_plot
dev.off()
