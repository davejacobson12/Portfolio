args <- commandArgs(TRUE)
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))



read.table("/Users/dave/Desktop/reference_based_mapping/delete.txt", sep = "\t", header =F) -> concatMatrix
#plot Gini-Simpson by Bray-Curtis for all pathways (using mean values for each population)
#Write that information to a file
as.data.frame(concatMatrix) -> concatDF
colnames(concatDF) <- c("Pathway","Pop", "Mean_Gini_Simpson", "Mean_Bray_Curtis")
ecolGG <- ggplot(data = concatDF, mapping = aes(x = concatDF$Mean_Gini_Simpson, y = concatDF$Mean_Bray_Curtis, shape = concatDF$Pathway))

ecolGG + geom_point(aes(colour = concatDF$Pop), size = 2) + ylim(0,1) + xlim(0,1) + xlab("Gini-Simpson") + ylab("Bray-Curtis") + scale_shape_manual(values = c(3,4,6,15:18)) + scale_colour_manual(values = c(Burkina_Faso ="darkblue", Euro = "#fdbf6f", Hadza = "darkviolet", Italy = "#ff7f00", Mongolian = "gray54", Norman = "#e31a1c", Matses = "cornflowerblue")) + theme(panel.background = element_rect(fill = "whitesmoke"), panel.grid.minor = element_blank())

ggsave("/Users/dave/Desktop/reference_based_mapping/contributional_div.png", width = 10, height = 6)
write.table(concatDF, "/Users/dave/Desktop/reference_based_mapping/ecolMetrics.txt", col.names = T, row.names =F, sep  = "\t", quote = F)

#write the pathway richness matrix to a final file
richMat <- read.table("/Users/dave/Desktop/reference_based_mapping/tempRich.txt", header = F , sep  = "\t")
as.data.frame(richMat) -> richDF
colnames(richDF) <- c("Pathway", "MeanRichness_BF", "MeanRichness_Italy", "MeanRichness_Hadza", "pval")
richDF$fdf <- p.adjust(richDF$pval, method = "fdr")
write.table(richDF, "/Users/dave/Desktop/reference_based_mapping/meanRichness.txt", col.names =T, sep  = "\t", row.names = F, quote = F)
