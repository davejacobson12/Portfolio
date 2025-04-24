#Alpha Diversity plots of fecal samples. Comparing study groups as well as changes with sampling time period for New ovarian

#Alpha diversity metrics were generated in qiime and imported into R as ocm1.88_alpha

#Used grep within R to pull out all of the fecal samples

# Used ggplot to plot PD and OTUs vs Study Group - As separate plots
	#They need to be edited in illustrator to remove the NA column and Alter category names

gg_ocm1.88_alpha_fecal <- ggplot(data = ocm1.88_alpha_fecal, aes(x = ocm1.88_alpha_fecal$study_group, y = ocm1.88_alpha_fecal$PD, fill = ocm1.88_alpha_fecal$study_group)

gg_ocm1.88_alpha_fecal + geom_boxplot() + theme_classic() + scale_fill_manual(values = brewer.pal(n= 6, "Set2"), guide = F) + ylab("Faith's PD") + xlab("Study Group") + ggtitle("Study Group Impacts Fecal Alpha Diversity") + theme(plot.title = element_text(hjust = 0.5))

gg_ocm1.88_alpha_fecal <- ggplot(data = ocm1.88_alpha_fecal, aes(x = ocm1.88_alpha_fecal$study_group, y = ocm1.88_alpha_fecal$OTUs, fill = ocm1.88_alpha_fecal$study_group)

gg_ocm1.88_alpha_fecal + geom_boxplot() + theme_classic() + scale_fill_manual(values = brewer.pal(n= 6, "Set2"), guide = F) + ylab("Observed OTUs") + xlab("Study Group") + ggtitle("Study Group Impacts Fecal Alpha Diversity") + theme(plot.title = element_text(hjust = 0.5))

#Test for significant differences for OTU and PD between groups

> OTU.aov <- aov(ocm1.88_alpha_fecal$OTUs ~ ocm1.88_alpha_fecal$study_group)
> summary(OTU.aov)
                                Df Sum Sq Mean Sq F value Pr(>F)
ocm1.88_alpha_fecal$study_group  4  15183    3796   1.231  0.304
Residuals                       80 246730    3084               
7 observations deleted due to missingness

No significant results

> PD.aov <- aov(ocm1.88_alpha_fecal$PD ~ ocm1.88_alpha_fecal$study_group)
> summary(PD.aov)
                                Df Sum Sq Mean Sq F value Pr(>F)  
ocm1.88_alpha_fecal$study_group  4  308.7   77.16   3.248  0.016 *
Residuals                       80 1900.3   23.75                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
7 observations deleted due to missingness

There is a significant	result, use TukeyHSD post-hoc to investigate further. 

TukeyHSD(PD.aov)

PFI6 vs new_ovarian p value 0.0099 * the only significant value
PFI6 vs PFI24 p value 0.0879 * only other comparison below 0.1

# I repeated the same analysis for only new ovarian patients, using the three sample collection periods as factors.

gg_ocm1.88_alpha_fecal_new_ovarian <- ggplot(ocm1.88_alpha_fecal_new_ovarian, aes(x=ocm1.88_alpha_fecal_new_ovarian$Sampling_period, y = ocm1.88_alpha_fecal_new_ovarian$PD, fill = ocm1.88_alpha_fecal_new_ovarian$Sampling_period))

gg_ocm1.88_alpha_fecal_new_ovarian + geom_boxplot() + scale_fill_manual(values = brewer.pal(n = 3, name = "Set2"), guide = F) + theme_classic() + xlab("Sampling Period") + ylab("Faith's PD") + ggtitle("Platinum Treatment Cycle and Fecal Alpha Diversity") + theme(plot.title = element_text(hjust = 0.5))

ANOVA for both OTUs and PD had no significant p value. 

Made a jitter plot to show differential abundance of E. coli and Lactobacillus between PFI groups

> gg_ocm1.88_L6_vaginal_escherchia + geom_jitter(aes(colour = ocm1.88_vaginal_L6_topTaxa$study_group))  + theme(legend.position = "none", plot.title=element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_text(aes(label = ocm1.88_vaginal_L6_topTaxa$patient_num)) + ylab("Taxa Abundance")  + ggtitle("Escherichia coli Abundance Differs by PFI Status") + xlab("PFI Status")

Edited the graph in Illustrator to make legend for sample site (already denoted by shapes)and color by patient number (only for the ones above 10% E. coli

###To show plots where each sample is a point, colored by sample site and then grouped in to different panels depending	on the Lactobacillus abundance and study group.
###They are grouped in high or low lactobacillus based on if ONE of the vaginal samples is above 50% abundance in Lactobacillus. 
ovarian_cq_gg <- ggplot(data = ovarian_cq, aes(x = patient_num, y = cq, color =BodySite))
ovarian_cq_gg + geom_jitter(position = position_jitter(0.01)) + facet_wrap(High_lacto~study_group)

### doing the same with E. coli. 
 ovarian_cq_gg + geom_jitter(position = position_jitter(0.01)) + facet_wrap(High_ecoli~study_group)

### Boxplot looking at each individual sample to see if low/high cq can be correlated with high/low lactobacillus and e. coli

cq_taxa_abund <- ggplot(taxa_abund, aes(x=High_ecoli, y = cq))

### With different community types. I grouped samples in diverse community as those where the following taxa all summed up to be greater than 50% of the reads. Samples that didn't fit diverse/lactobacillus/ecoli were placed in unresolved
#Gardnerella, Atopobium, Corynebacterium, Bifidobacterium, dialister, prevotella, streptococcus, bacteroides, and megasphera.-Taken from albert 2015. 
cq_taxa_abund_gg <- ggplot(taxa_abund, aes(x=community, y = cq))
cq_taxa_abund_gg + geom_boxplot() + geom_jitter(position=position_jitter(0.2), aes(color=study_group))


