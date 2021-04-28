library(data.table)
library(DESeq2)
library(tidyr)
library(ggplot2)
library(viridis)
library(stringr)
library(dplyr)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")

#################
## add factors ##
#################
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$treatment <- factor(paste(asw_dds$Treatment))
asw_dds$pc1 <- factor(paste(asw_dds$PC1_sign))
##set reference levels
asw_dds$location <- relevel(asw_dds$location, ref="Invermay")
asw_dds$treatment <- relevel(asw_dds$treatment, ref="NC")
asw_dds_int <- copy(asw_dds)
##remove genes with less than 30 counts
keep <- rowSums(counts(asw_dds_int)) >= 30
asw_dds_int <- asw_dds_int[keep,]

################
## run DESeq2 ##
################
##control for all factors and test for location-specific responses to exposure
design(asw_dds_int) <- ~pc1+location+treatment+location:treatment
##want interaction between location:treatment
##WT as all factors only have 2 levels - LRT is for 3+
asw_dds_int <- DESeq(asw_dds_int, test="Wald")
saveRDS(asw_dds_int, "output/deseq2/asw/interaction/dds_interaction.rds")
