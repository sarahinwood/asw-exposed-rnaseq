library(data.table)
library(DESeq2)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/best_annot_per_gene.csv")

asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds_location <- copy(asw_dds)
design(asw_dds_location) <- ~location
asw_dds_location <- DESeq(asw_dds_location)

res_group <- results(asw_dds_location, contrast = c("location", "Ruakura", "Invermay"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/deseq2/asw/location/res_group.csv")

##Make sig DEG data table and write to output
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw/location/sig_w_annots.csv")
saveRDS(asw_dds_location, "output/deseq2/asw/location/asw_dds_location.rds")

plotCounts(asw_dds_location, "", intgroup = c("location"), main="")

                    