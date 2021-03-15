library(tximport)
library(data.table)
library(DESeq2)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv")

asw_dds$pc1 <- factor(paste(asw_dds$PC1_sign))
asw_dds_pc1 <- copy(asw_dds)
design(asw_dds_pc1) <- ~pc1
asw_dds_pc1 <- DESeq(asw_dds_pc1)

res_group <- results(asw_dds_pc1, contrast = c("pc1", "positive", "negative"), lfcThreshold = 1, alpha = 0.1)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/deseq2/asw/pc1/res_group.csv")
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw/pc1/sig_annots.csv")
saveRDS(asw_dds_pc1, "output/deseq2/asw/pc1/asw_dds_pc1.rds")

plotCounts(asw_dds_pc1, "TRINITY_DN14750_c0_g1", intgroup = c("pc1"))

