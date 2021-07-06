library(data.table)
library(tidyr)
library(DESeq2)
library(EnhancedVolcano)

############
## set up ##
############
asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = "")
blast <- fread("output/deseq2/asw/unann_degs/interaction_blastx_unann_res.csv")

asw_dds$treatment <- factor(paste(asw_dds$Treatment))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$pc1 <- factor(paste(asw_dds$PC1_sign))
##set reference levels
asw_dds$location <- relevel(asw_dds$location, ref="Invermay")
asw_dds$treatment <- relevel(asw_dds$treatment, ref="NC")
asw_dds_exposed <- copy(asw_dds)

#keep <- rowSums(counts(asw_dds_exposed)) >= 30
#asw_dds_exposed <- asw_dds_exposed[keep,]

##control for pc1, compare treatment levels
design(asw_dds_exposed) <- ~pc1+location+treatment
asw_dds_exposed <- DESeq(asw_dds_exposed)

#############
## results ##
#############
res_group <- results(asw_dds_exposed, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
##Make sig DEG data table
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
##write output
fwrite(ordered_res_group_table, "output/deseq2/asw/all_exposed/res_group.csv")
fwrite(sig_annots, "output/deseq2/asw/all_exposed/sig_w_annots.csv")
saveRDS(asw_dds_exposed, "output/deseq2/asw/all_exposed/asw_dds_exposed.rds")
##nr blast annots
all_annots <- merge(sig_annots, blast, by="rn", all.x=TRUE)
fwrite(all_annots, "output/deseq2/asw/all_exposed/all_sig_annots.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 3, pCutoff = 0.05, colAlpha=0.5,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

##plot counts for genes of interest, sub in name
plotCounts(asw_dds_exposed, "TRINITY_DN60874_c0_g1", intgroup = c("treatment"))

##########
## plot ##
##########
asw_dds_location <- readRDS("output/deseq2/asw/all_exposed/asw_dds_exposed.rds")

##genes with some annotation
blastx_trino_annots <- subset(all_annots, !is.na(all_annots$sprot_Top_BLASTX_hit))
blast_nr_annots <- subset(all_annots, !is.na(all_annots$annotation))
degs_annots <- full_join(blastx_trino_annots, blast_nr_annots)                        
id_annot <- degs_annots[,c(1,9)]

##get gene counts
counts_table <- data.table(counts(asw_dds_location, normalized=TRUE), keep.rownames = TRUE)
annot_counts <- filter(counts_table, rn %in% degs_annots$rn)
##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:25], key="sample_name", value="normalized_counts")
plot_annots_counts <- full_join(plot_annots_counts, id_annot)
##sample group information
sample_table <- fread("data/sample_table.csv")
name_vs_group <- sample_table[,c(1,5,6)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)

##plot all annot DEGs using ggplot2
ggplot(plotting_counts) +
  geom_point(aes(x = Treatment, y = normalized_counts, colour=Treatment), alpha=0.5) +
  labs(colour="Weevil Location", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~plot_label, scales="free")
