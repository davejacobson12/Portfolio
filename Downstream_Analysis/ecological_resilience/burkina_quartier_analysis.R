#16s analysis for presentation in R

#Load into R L6 table with proportion values for abundance of taxa and metadata table

library(vegan)
library(reshape)
library(ggplot2)

#Calculate hill numbers from this data
exp(renyi(t(murazHill))) -> murazBind 
#remove columns above a q of 4
murazBind[,1:5] -> murazBind
#Rename columns
colnames(murazBind) <- c("murazHill0", "murazHillquart", "murazHillhalf", "murazHill1", "murazHill2")
#melt
melt(murazBind) -> murazHillMelt
as.data.frame(murazHillMelt) -> murazMeltdf
#merge with metadata
merge(x=murazMeltdf, y = bftm_adults_alpha) -> murazHillMerge_adults
merge(x= murazMeltdf, y = bftm_topQuartiers) -> murazHillMerge_quartiers
#ggplot for top quartiers by top quartier 
p9 <- ggplot(murazHillMerge_quartiers, mapping = aes(x = murazHillMerge_quartiers$order_q, y =murazHillMerge_quartiers$value, colour = murazHillMerge_quartiers$quartier_name))
ggSmoothOut <-  p9 + geom_smooth(se = FALSE) + scale_color_manual(values = c("#b2182b","#ef8a62","purple","grey19","#67a9cf","#2166ac")) + scale_y_continuous(breaks = seq(0,150,25)) + theme(legend.position = c(0.8, 0.8), panel.grid.minor = element_blank(), axis.title = element_text(size=18), axis.text =element_text(size  = 14)) + xlab("Order q") + ylab("Hill Number")
p10 <- ggplot(murazHillquart_integers, mapping = aes(x = murazHillquart_integers$order_factor, y =murazHillquart_integers$value, colour = murazHillquart_integers$quartier_name))
# p10+geom_boxplot()
ggBoxOut <- p10 + geom_boxplot() + scale_color_manual(values = c("#b2182b","#ef8a62","purple","grey19","#67a9cf","#2166ac")) + scale_y_continuous(breaks = seq(0,225,25)) + theme(legend.position = c(0.8, 0.8), panel.grid.minor = element_blank(), axis.title = element_text(size=18), axis.text =element_text(size  = 14)) + xlab("Order q") + ylab("Hill Number") + ylim(c(0,175))

#Kruskal and Dunn for q =0, 1, 2

################# q = 0 ########################################

# Kruskal-Wallis chi-squared = 11.735, df = 5, p-value = 0.03861


# Comparison           Z     P.unadj      P.adj
# 1      badama - koalayiri  0.17687821 0.859604058 0.92100435
# 2     badama - mougounsin  0.79895103 0.424318806 0.53039851
# 3  koalayiri - mougounsin  0.67723638 0.498255992 0.57491076
# 4         badama - sambin  1.73531862 0.082684386 0.20671096
# 5      koalayiri - sambin  1.68330803 0.092315476 0.19781888
# 6     mougounsin - sambin  0.88120404 0.378207395 0.51573736
# 7       badama - tagsyiri  2.68950261 0.007155859 0.05366894
# 8    koalayiri - tagsyiri  2.70910818 0.006746434 0.10119651
# 9   mougounsin - tagsyiri  1.85155570 0.064089650 0.32044825
# 10      sambin - tagsyiri  1.06280045 0.287872457 0.47978743
# 11      badama - yissouka  1.83583780 0.066381650 0.24893119
# 12   koalayiri - yissouka  1.79722235 0.072300328 0.21690098
# 13  mougounsin - yissouka  0.94258339 0.345894041 0.51884106
# 14      sambin - yissouka  0.02285988 0.981762040 0.98176204
# 15    tagsyiri - yissouka -1.09529019 0.273389525 0.51260536


################# q = 1 ########################################

# Kruskal-Wallis chi-squared = 18.671, df = 5, p-value = 0.002213

# Comparison          Z      P.unadj       P.adj
# 1      badama - koalayiri  0.7773224 0.4369685885 0.595866257
# 2     badama - mougounsin  0.4648791 0.6420180734 0.802522592
# 3  koalayiri - mougounsin -0.2803458 0.7792122310 0.834870247
# 4         badama - sambin  2.3924964 0.0167341935 0.050202581
# 5      koalayiri - sambin  1.7445873 0.0810567342 0.173693002
# 6     mougounsin - sambin  1.8955198 0.0580235761 0.145058940
# 7       badama - tagsyiri  3.4936781 0.0004764151 0.007146226
# 8    koalayiri - tagsyiri  2.9407171 0.0032745346 0.016372673
# 9   mougounsin - tagsyiri  3.0061089 0.0026461425 0.019846069
# 10      sambin - tagsyiri  1.2344771 0.2170251917 0.361708653
# 11      badama - yissouka  0.9808626 0.3266604870 0.489990731
# 12   koalayiri - yissouka  0.1828075 0.8549490577 0.854949058
# 13  mougounsin - yissouka  0.4611121 0.6447182150 0.743905633
# 14      sambin - yissouka -1.6561490 0.0976916549 0.183171853
# 15    tagsyiri - yissouka -2.9139686 0.0035686577 0.013382466


################# q = 2 ########################################

# Kruskal-Wallis chi-squared = 21.849, df = 5, p-value = 0.0005593

# Comparison          Z      P.unadj       P.adj
# 1      badama - koalayiri  0.5861775 0.5577562133 0.760576655
# 2     badama - mougounsin  0.3218393 0.7475744074 0.862585855
# 3  koalayiri - mougounsin -0.2421168 0.8086896532 0.808689653
# 4         badama - sambin  2.6027557 0.0092477823 0.027743347
# 5      koalayiri - sambin  2.1781534 0.0293946206 0.062988473
# 6     mougounsin - sambin  2.2586950 0.0239023651 0.059755913
# 7       badama - tagsyiri  3.5482313 0.0003878275 0.005817412
# 8    koalayiri - tagsyiri  3.2013799 0.0013677108 0.006838554
# 9   mougounsin - tagsyiri  3.2106834 0.0013241975 0.009931482
# 10      sambin - tagsyiri  1.0711039 0.2841227225 0.473537871
# 11      badama - yissouka  0.8369317 0.4026309592 0.603946439
# 12   koalayiri - yissouka  0.2437433 0.8074296023 0.865103145
# 13  mougounsin - yissouka  0.4771044 0.6332877972 0.791609747
# 14      sambin - yissouka -2.0522319 0.0401471306 0.075275870
# 15    tagsyiri - yissouka -3.1286171 0.0017563102 0.006586163



################# Specific Taxa by Quartier ########################################

##I'll just show 3 taxa

################# Oscillibacter ########################################

# Kruskal-Wallis chi-squared = 24.186, df = 5, p-value = 2e-04

# Comparison          Z      P.unadj        P.adj
# 1      badama - koalayiri  0.3870794 6.986974e-01 0.6986974253
# 2     badama - mougounsin -0.7211792 4.707993e-01 1.0000000000
# 3  koalayiri - mougounsin -1.1580524 2.468427e-01 1.0000000000
# 4         badama - sambin  2.3846003 1.709769e-02 0.1709768774
# 5      koalayiri - sambin  2.1575691 3.096135e-02 0.2786521201
# 6     mougounsin - sambin  3.1555732 1.601830e-03 0.0208237957
# 7       badama - tagsyiri  3.2579424 1.122232e-03 0.0157112437
# 8    koalayiri - tagsyiri  3.0991309 1.940892e-03 0.0232907066
# 9   mougounsin - tagsyiri  4.0143215 5.961703e-05 0.0008942555
# 10      sambin - tagsyiri  0.9889867 3.226696e-01 1.0000000000
# 11      badama - yissouka  0.7836464 4.332476e-01 1.0000000000
# 12   koalayiri - yissouka  0.4124168 6.800339e-01 1.0000000000
# 13  mougounsin - yissouka  1.5899492 1.118463e-01 0.7829237691
# 14      sambin - yissouka -1.8618607 6.262273e-02 0.5009818025
# 15    tagsyiri - yissouka -2.8565277 4.283026e-03 0.0471132884

################# Prevotella ########################################

# Kruskal-Wallis chi-squared = 26.686, df = 5, p-value = 6.567e-05

# Comparison           Z      P.unadj       P.adj
# 1      badama - koalayiri -0.45874764 6.464154e-01 1.000000000
# 2     badama - mougounsin  0.63175871 5.275446e-01 1.000000000
# 3  koalayiri - mougounsin  1.13412611 2.567417e-01 1.000000000
# 4         badama - sambin -2.51355478 1.195212e-02 0.107569094
# 5      koalayiri - sambin -2.21944537 2.645644e-02 0.211651517
# 6     mougounsin - sambin -3.18893325 1.427988e-03 0.017135862
# 7       badama - tagsyiri -3.23227734 1.228078e-03 0.015965014
# 8    koalayiri - tagsyiri -2.99578670 2.737378e-03 0.030111163
# 9   mougounsin - tagsyiri -3.89487147 9.825081e-05 0.001473762
# 10      sambin - tagsyiri -0.82512630 4.092999e-01 1.000000000
# 11      badama - yissouka -0.03198465 9.744843e-01 0.974484291
# 12   koalayiri - yissouka  0.48748669 6.259135e-01 1.000000000
# 13  mougounsin - yissouka -0.73831236 4.603247e-01 1.000000000
# 14      sambin - yissouka  2.82698753 4.698816e-03 0.046988158
# 15    tagsyiri - yissouka  3.62586027 2.880011e-04 0.004032015


Prevotella_gg <- ggplot(data = bftm_topQuartiers, aes(y = bftm_topQuartiers$Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotella, x = bftm_topQuartiers$quartier_name))

Prevotella_gg + geom_boxplot(mapping = aes(color = bftm_topQuartiers$quartier_name ))  + scale_color_manual(values = c("#b2182b","#ef8a62","purple","grey19","#67a9cf","#2166ac")) + theme(legend.position = c(0.125, 0.8), panel.grid.minor = element_blank(), axis.title = element_text(size=18), axis.text =element_text(size  = 14), title = element_text(hjust = 0.5)) + xlab("Compound") + ylab("Prevotella Abundance") + ggtitle("Prevotella Abundance")

################# Sporobacter ########################################

# Kruskal-Wallis chi-squared = 23.247, df = 5, p-value = 0.0003028

# Comparison          Z      P.unadj        P.adj
# 1      badama - koalayiri -0.1051367 9.162673e-01 0.9162673454
# 2     badama - mougounsin -0.7629294 4.455055e-01 1.0000000000
# 3  koalayiri - mougounsin -0.7104692 4.774133e-01 0.9548265047
# 4         badama - sambin  1.5627135 1.181200e-01 0.8268399114
# 5      koalayiri - sambin  1.8014841 7.162660e-02 0.6446394309
# 6     mougounsin - sambin  2.3783194 1.739175e-02 0.1913093000
# 7       badama - tagsyiri  3.3279678 8.748199e-04 0.0113726582
# 8    koalayiri - tagsyiri  3.6944997 2.203202e-04 0.0030844821
# 9   mougounsin - tagsyiri  4.1281349 3.657176e-05 0.0005485764
# 10      sambin - tagsyiri  1.9326134 5.328383e-02 0.5328383470
# 11      badama - yissouka  0.8596449 3.899848e-01 1.0000000000
# 12   koalayiri - yissouka  1.0555662 2.911664e-01 1.0000000000
# 13  mougounsin - yissouka  1.7126259 8.678138e-02 0.6942510310
# 14      sambin - yissouka -0.8433648 3.990245e-01 1.0000000000
# 15    tagsyiri - yissouka -2.8550924 4.302431e-03 0.0516291679

kruskal.test(bftm_topQuartiers$Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Oscillibacter ~ bftm_topQuartiers$quartier_name)

Sporobacter_gg <- ggplot(data = bftm_topQuartiers, aes(y = bftm_topQuartiers$Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Sporobacter, x = bftm_topQuartiers$quartier_name))

Sporobacter_gg + geom_boxplot(mapping = aes(color = bftm_topQuartiers$quartier_name ))  + scale_color_manual(values = c("#b2182b","#ef8a62","purple","grey19","#67a9cf","#2166ac")) + theme(legend.position = c(0.8, 0.8), panel.grid.minor = element_blank(), axis.title = element_text(size=18), axis.text =element_text(size  = 14), title = element_text(hjust = 0.5)) + xlab("Compound") + ylab("Sporobacter Abundance") + ggtitle("Sporobacter Abundance")
                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                        
########################### Make tree from strain tracking #######################                                                                                                                                                                                                                                                                                                                                                                             
args <- commandArgs(TRUE)
suppressMessages(library(ape))
suppressMessages(library(picante))
suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))
suppressMessages(library(reshape2))
suppressMessages(library(ggstance))
suppressMessages(library(RColorBrewer))

##make bootstrap tree in mega
#read in that newick tree
xdd <- read.tree("/Users/dave/Desktop/aapa2019_presentation/strain_analysis/E_rectale_MSA/Newick Export.nwk")
ggtree(xdd) -> da
###read in metadata with first column that has sample names as bftm####
###add row name for reference genome
as.character(df$treeName) -> df$treeName
df[120,1] <- "GCA_000020605.1_ASM2060v1_genomic"
daMerge <- da %<+% df
attach(df)
### start plotting based on metadata cats
daMerge + geom_tiplab(aes(fill = primary_water), geom = "label") + theme(legend.position = "right") 



####################### Heatmap of all pathways ###################################
#quartier is the only metadata category with significant differences by pathway abundance (90 pathways, vs 0 in other metadata)
#filtered the dataset to have pathways that have top 1/3 abundance

samples <- test1$treeName
correctOrder <- vector(mode ="numeric", length = 61)
for (i in 1:61){
grep(samples[i], colnames(quartierSig)) -> correctOrder[i]
}
quartierSig[,correctOrder] -> nsn
pathDepth <- vector(mode = "logical", length = 88)
for (i in 1:88){
sum(quartierSig[i,]) > 0.005 -> pathDepth[i]
}

nsn[grep("TRUE", pathDepth),] -> quartierSig_topPaths

mapDrugToColor <- function(annotations){
  colorsVector <- ifelse(annotations["quartier_name"]=="badama", "#b2182b", ifelse(annotations["quartier_name"]=="koalayiri", "#ef8a62", ifelse(annotations["quartier_name"]=="yissouka","#2166ac", ifelse(annotations["quartier_name"] == "mougounsin", "purple", ifelse(annotations["quartier_name"] == "tagsyiri","#67a9cf" , ifelse(annotations["quartier_name"] == "sambin", "grey19", "black"))))))
  return(colorsVector)
}
                                                                                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                              
mapDrugToColor(test1)
color <- colorRampPalette(c("red", "yellow", "white"))(1000)

testHeatmap<-function(paths, annotations) {    
sampleColors = mapDrugToColor(annotations)
heatmap.2(paths, margins=c(5,8), ColSideColors=sampleColors, col = color, trace ="none", Colv = F)
}

testHeatmap(test3, test2)

##test3 and test2 are the ordered versions that I want (I forgot to copy the exact code to get those orders, but it shouldn't be difficult to replicate if I need to)


### Taxa contributing to different pathways###
#two pathways that are high abundance in tagsyiri and sambin (2942 and 6700)
#two pathways that are low abundance in tagsyiri and sambin (5177, P124)
#showing top 9 taxa plus condensed down "other"
#undecided about whether I should show unclassified abundances (the "other" is abundance where taxa was identified with that pathway)


####################### Random Forest for 16S quartier membership ###################################
# only going to use the top 6 quartiers 
# at the L6 level, remove taxa that are 0 abundance in 10% of samples, which was 232 taxa
#collapase sambin and tagsyiri into non-diverse and other quartiers into diverse ecology 
                                                                                                                                                                                                                                                                                                                                                                              
library(randomForest)
library(rfUtilities)
library(caret)
library(plyr)
library(e1071)


#Normalize abundances and merge with metadata category 
quartierRF_norm <- sweep(quartierRFotu, 2, colSums(quartierRFotu) , '/') * 100
quartierRF_scale <- scale(quartierRF_norm, center = TRUE, scale = TRUE)  
quartierRF_scaled_state <- data.frame(t(quartierRF_scale))  
quartierRF_scaled_state$quartier <- samles[rownames(quartierRF_scaled_state), "ecol"]


# run random forest
set.seed(992)
quartierRF_state_classify <- randomForest( x=quartierRF_scaled_state[,1:(ncol(quartierRF_scaled_state)-1)] , y=quartierRF_scaled_state[ , ncol(quartierRF_scaled_state)] , ntree=501, importance=TRUE, proximities=TRUE )

#compare to random 
quartier_state_classify_sig <- rf.significance( x=quartierRF_state_classify ,  xdata= quartierRF_scaled_state[,1:(ncol(quartierRF_scaled_state)-1)] , nperm=1000 , ntree=501 )  
#results
# Number of permutations:  1000 
# p-value:  0 
# Model signifiant at p = 0 
# Model OOB error:  0.1617647 
# Random OOB error:  0.3970588 
# min random global error: 0.2205882 
# max random global error:  0.5147059 
# min random within class error: 0.6111111 
# max random within class error:  0.6111111

#when using top 6 quarter names independently, confusion matrix shows misclassification, mostly within the two groups
# Type of random forest: classification
# Number of trees: 501
# No. of variables tried at each split: 15

# OOB estimate of  error rate: 67.65%
# Confusion matrix:
# badama koalayiri mougounsin sambin tagsyiri yissouka class.error
# badama          0         3          0      0        0        6   1.0000000
# koalayiri       3         1          0      1        1        6   0.9166667
# mougounsin      1         4          0      0        1        3   1.0000000
# sambin          1         3          0      2        5        1   0.8333333
# tagsyiri        0         0          0      2        9        0   0.1818182
# yissouka        1         0          1      2        1       10   0.3333333

## to see what taxa are important features
quartier_state_classify_imp <- as.data.frame(quartierRF_state_classify$importance )
quartier_state_classify_imp$features <- rownames( quartier_state_classify_imp )
quartier_state_classify_imp_sorted <- arrange( quartier_state_classify_imp , desc(MeanDecreaseAccuracy)  )
barplot(quartier_state_classify_imp_sorted$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")
#top ten most important taxa
barplot(quartier_state_classify_imp_sorted[1:10,"MeanDecreaseAccuracy"], names.arg=quartier_state_classify_imp_sorted[1:10,"features"] , ylab="Mean Decrease in Accuracy (Variable Importance)", las=2, ylim=c(0,0.02), main="Classification RF")  
##can look at this table to see which taxa have contributions
quartier_state_classify_imp_sorted
                                                                                                                                                                                                                                                                                                                                                                              
##leave one out cross val
fit_control <- trainControl( method = "LOOCV" )    

Qu_state_classify_loocv <- train(quartierRF_scaled_state[,1:(ncol(quartierRF_scaled_state)-1)] , y=quartierRF_scaled_state[ , ncol(quartierRF_scaled_state)]  , method="rf", ntree=501 , tuneGrid=data.frame( mtry=15 ) , trControl=fit_control )

# Resampling: Leave-One-Out Cross-Validation 
# Summary of sample sizes: 67, 67, 67, 67, 67, 67, ... 
# Resampling results:

# Accuracy   Kappa    
# 0.8235294  0.5972359

# Tuning parameter 'mtry' was held constant at a value of 15



