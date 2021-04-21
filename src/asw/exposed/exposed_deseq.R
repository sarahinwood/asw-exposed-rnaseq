library(data.table)
library(DESeq2)
library(EnhancedVolcano)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = "")

asw_dds$treatment <- factor(paste(asw_dds$Treatment))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$pc1 <- factor(paste(asw_dds$PC1_sign))
asw_dds_exposed <- copy(asw_dds)
keep <- rowSums(counts(asw_dds_exposed)) >= 30
asw_dds_exposed <- asw_dds_exposed[keep,]

##control for pc1, compare treatment levels
design(asw_dds_exposed) <- ~pc1+location+treatment
asw_dds_exposed <- DESeq(asw_dds_exposed)

res_group <- results(asw_dds_exposed, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/deseq2/asw/exposed/res_group.csv")

##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw/exposed/sig_annots.csv")
saveRDS(asw_dds_exposed, "output/deseq2/asw/exposed/asw_dds_exposed.rds")

plotCounts(asw_dds_exposed, "TRINITY_DN16993_c0_g1", intgroup="treatment")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="", subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2, col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

##blast results
unann_blast <- fread("output/deseq2/asw/unann_degs/blastx_unann_res.csv")
unann_blast <- unann_blast[,c(1,3,11,13)]
sig_blast_annots <- merge(sig_annots, unann_blast, by="transcript_id", all.x=TRUE)
fwrite(sig_blast_annots, "output/deseq2/asw/exposed/sig_blast_annots.csv")
