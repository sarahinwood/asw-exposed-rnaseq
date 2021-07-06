library(data.table)
library(DESeq2)
library(tidyr)
library(ggplot2)
library(viridis)
library(stringr)
library(dplyr)
library(EnhancedVolcano)

##difference between exposed Ruakura and Invermay

asw_dds_int <- readRDS("output/deseq2/asw/interaction/dds_interaction.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = "")
blast <- fread("output/deseq2/asw/unann_degs/interaction_blastx_unann_res.csv", na.strings = "")

#############
## results ##
#############
resultsNames(asw_dds_int)
##effect of location in Invermay ASW
res_group_location_exposed <- results(asw_dds_int, alpha=0.05, lfcThreshold = 1, list(c("location_Ruakura_vs_Invermay", "locationRuakura.treatmentExposed")))
mcols <- data.frame(mcols(res_group_location_exposed))
summary(res_group_location_exposed)
ordered_res_group_location_exposed <- res_group_location_exposed[order(res_group_location_exposed$padj),]
ordered_res_group_table_location_exposed <- data.table(data.frame(ordered_res_group_location_exposed), keep.rownames = TRUE)
##sig degs
ordered_sig_res_group_table <- subset(ordered_res_group_table_location_exposed, padj < 0.05)
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
all_annots <- merge(sig_annots, blast, by="rn", all.x=TRUE)
fwrite(ordered_res_group_table_location_exposed, "output/deseq2/asw/interaction/location_exposure/res_group.csv")
fwrite(sig_annots, "output/deseq2/asw/interaction/location_exposure/sig_w_annots.csv")
fwrite(all_annots, "output/deseq2/asw/interaction/location_exposure/all_sig_annots.csv")

EnhancedVolcano(ordered_res_group_table_location_exposed, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 3, pCutoff = 0.05, colAlpha=0.5,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

##plot counts for genes of interest, sub in name
plotCounts(asw_dds_int, "TRINITY_DN13158_c0_g1", intgroup = c("location", "treatment"))

###############################
## plot multiple gene counts ##
###############################
##filter for annot genes
annotated_genes_trin_blastx <- subset(all_annots, !is.na(all_annots$sprot_Top_BLASTX_hit))
annotated_genes_trin_blastp <- subset(all_annots, !is.na(all_annots$sprot_Top_BLASTP_hit))
annotated_genes_nr_blast <- subset(all_annots, !is.na(all_annots$annotation))
all_annotated_trin <- full_join(annotated_genes_trin_blastx, annotated_genes_trin_blastp)
all_annotated <- full_join(all_annotated_trin, annotated_genes_nr_blast)
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
plotting_counts$group <- str_replace_all(plotting_counts$group, "Exposed", "exposed")
plotting_counts$group <- str_replace_all(plotting_counts$group, "NC", "negative control")
##order so that two ex groups side by side in plot as these are samples being compared
group_order <- c("Invermay negative control", "Ruakura negative control", "Invermay exposed", "Ruakura exposed")
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
  facet_wrap(~gene_label, scales="free")


