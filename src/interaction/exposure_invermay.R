library(data.table)
library(DESeq2)
library(tidyr)
library(ggplot2)
library(viridis)
library(stringr)
library(dplyr)

##effect of exposure on Invermay ASW

asw_dds_int <- readRDS("output/deseq2/asw/location_exposure_pc1_int/dds_interaction.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = "")

#############
## results ##
#############
resultsNames(asw_dds_int)
##effect of exposure in Invermay ASW
res_group_exposure <- results(asw_dds_int, alpha=0.05, lfcThreshold = 1, contrast=c("treatment", "Exposed", "NC"))
summary(res_group_exposure)
ordered_res_group_exposure <- res_group_exposure[order(res_group_exposure$padj),]
ordered_res_group_table_exposure <- data.table(data.frame(ordered_res_group_exposure), keep.rownames = TRUE)
##sig degs
ordered_sig_res_group_table <- subset(ordered_res_group_table_exposure, padj < 0.05)
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(ordered_res_group_table_exposure, "output/deseq2/asw/interaction/exposure_invermay/res_group.csv")
fwrite(sig_annots, "output/deseq2/asw/interaction/exposure_invermay/sig_w_annots.csv")

##plot counts for genes of interest, sub in name
plotCounts(asw_dds_int, "TRINITY_DN14886_c0_g1", intgroup = c("location", "treatment"))

###############################
## plot multiple gene counts ##
###############################

##get gene counts
counts_table <- data.table(counts(asw_dds_int, normalized=TRUE), keep.rownames = TRUE)
annot_counts <- filter(counts_table, rn %in% sig_annots$rn)
##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:25], key="sample_name", value="normalized_counts")
##sample group information
sample_table <- fread("data/sample_table.csv")
sample_table$group <- paste(sample_table$Weevil_Location, sample_table$Treatment, sep=" ")
name_vs_group <- sample_table[,c(1,17)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)
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
  facet_wrap(~gene_label, scales="free", ncol=3)


