library(data.table)
library(DESeq2)

sample_data <- fread("data/sample_table.csv", header=TRUE)
asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")

asw_dds$group <- factor(paste(asw_dds$Treatment, asw_dds$Weevil_Location, sep="_"))
asw_dds$pc1 <- factor(paste(asw_dds$PC1_sign))


#PCA
asw_vst <- varianceStabilizingTransformation(asw_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(asw_vst, intgroup=c("group"))