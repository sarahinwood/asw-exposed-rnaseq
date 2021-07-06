library(data.table)
library(DESeq2)
library(tidyr)
library(ggplot2)
library(viridis)
library(stringr)
library(dplyr)
library(EnhancedVolcano)

##interaction term

asw_dds_int <- readRDS("output/deseq2/asw/interaction/dds_interaction.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = "")
blast <- fread("output/deseq2/asw/unann_degs/interaction_blastx_unann_res.csv", na.strings = "")

#############
## results ##
#############
resultsNames(asw_dds_int)
##effect of location in Invermay ASW
res_group_interaction <- results(asw_dds_int, alpha=0.05, lfcThreshold = 1)
mcols <- data.frame(mcols(res_group_interaction))
summary(res_group_interaction)
ordered_res_group_interaction <- res_group_interaction[order(res_group_interaction$padj),]
ordered_res_group_table_interaction <- data.table(data.frame(ordered_res_group_interaction), keep.rownames = TRUE)
##sig degs
ordered_sig_res_group_table <- subset(ordered_res_group_table_interaction, padj < 0.05)
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
all_annots <- merge(sig_annots, blast, by="rn", all.x=TRUE)
fwrite(ordered_res_group_table_interaction, "output/deseq2/asw/interaction/interaction/res_group.csv")
fwrite(sig_annots, "output/deseq2/asw/interaction/interaction/sig_w_annots.csv")
fwrite(all_annots, "output/deseq2/asw/interaction/interaction/all_sig_annots.csv")

EnhancedVolcano(ordered_res_group_table_interaction, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 3, pCutoff = 0.05, colAlpha=0.5,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

##plot counts for genes of interest, sub in name
plotCounts(asw_dds_int, "TRINITY_DN22325_c0_g2", intgroup = c("location", "treatment"))

###############################
## plot multiple gene counts ##
###############################
all_annots <- fread("output/deseq2/asw/interaction/interaction/all_sig_annots_plot.csv", na.strings="")
gene_annot <- all_annots[,c(1,9)]
##filter for annot genes
all_annotated <- subset(all_annots, !is.na(all_annots$plot_annot))
##get gene counts
counts_table <- data.table(counts(asw_dds_int, normalized=TRUE), keep.rownames = TRUE)
annot_counts <- filter(counts_table, rn %in% all_annotated$rn)
##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:25], key="sample_name", value="normalized_counts")
##sample group information
sample_table <- fread("data/sample_table.csv")
sample_table$group <- paste(sample_table$Weevil_Location, sample_table$Treatment, sep=" ")
name_vs_group <- sample_table[,c(1,17)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)
plotting_counts <- inner_join(plotting_counts, gene_annot)
plotting_counts$group <- str_replace_all(plotting_counts$group, "Exposed", "exposed")
plotting_counts$group <- str_replace_all(plotting_counts$group, "NC", "negative control")
group_order <- c("Invermay negative control", "Invermay exposed", "Ruakura negative control", "Ruakura exposed")
plotting_counts$group <- factor(plotting_counts$group, levels=group_order)
##add alphabetical label to each plot
plotting_counts$gene_label <- paste(plotting_counts$rn)
##plot all annot DEGs using ggplot2
ggplot(plotting_counts) +
  geom_point(aes(x = group, y = normalized_counts, colour=group)) +
  labs(colour="Sample group", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~plot_annot, scales="free", ncol=3)


