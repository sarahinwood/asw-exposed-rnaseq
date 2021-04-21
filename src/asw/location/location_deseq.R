library(data.table)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = "")

asw_dds$treatment <- factor(paste(asw_dds$Treatment))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$pc1 <- factor(paste(asw_dds$PC1_sign))
asw_dds_location <- copy(asw_dds)
keep <- rowSums(counts(asw_dds_location)) >= 30
asw_dds_location <- asw_dds_location[keep,]

##control for pc1, compare treatment levels
design(asw_dds_location) <- ~pc1+treatment+location
asw_dds_location <- DESeq(asw_dds_location)

res_group <- results(asw_dds_location, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
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

plotCounts(asw_dds_location, "TRINITY_DN39424_c0_g1", intgroup = c("location"), main="")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

#########################
## plot cooks distance ##
#########################
cooks_distance <- data.frame(log10(assays(asw_dds_location)[["cooks"]]))
melted_cooks_dist <- cooks_distance %>% gather(colnames(cooks_distance)[1:24], key="sample_name", value="cooks_distance")
sample_table <- fread("data/sample_table.csv")
sample_table$group <- paste(sample_table$Weevil_Location, sample_table$Treatment, sep=" ")
name_to_group <- sample_table[,c(1, 17)]
cooks_plot <- merge(melted_cooks_dist, name_to_group, all.x=TRUE)
##plot
ggplot(cooks_plot) +
  geom_boxplot(aes(x = sample_name, y = cooks_distance, colour=group), outlier.alpha=0.08) +
  labs(y="Cook's distance", x="", colour="Sample group")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())

boxplot(log10(assays(asw_dds_location)[["cooks"]]), range=0, las=2)

##blast results
unann_blast <- fread("output/deseq2/asw/unann_degs/blastx_unann_res.csv")
unann_blast <- unann_blast[,c(1,3,11,13)]
sig_blast_annots <- merge(sig_annots, unann_blast, by="transcript_id", all.x=TRUE)
fwrite(sig_blast_annots, "output/deseq2/asw/location/sig_blast_annots.csv")

##########
## plot ##
##########

##genes with some annotation
blastx_trino_annots <- subset(sig_blast_annots, !is.na(sig_blast_annots$sprot_Top_BLASTX_hit))
blast_nr_annots <- subset(sig_blast_annots, !is.na(sig_blast_annots$annotation))
degs_annots <- full_join(blastx_trino_annots, blast_nr_annots)                        

##get gene counts
counts_table <- data.table(counts(asw_dds_location, normalized=TRUE), keep.rownames = TRUE)
annot_counts <- filter(counts_table, rn %in% degs_annots$rn)
##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:25], key="sample_name", value="normalized_counts")
##sample group information
sample_table <- fread("data/sample_table.csv")
name_vs_group <- sample_table[,c(1,5)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)

##plot all annot DEGs using ggplot2
ggplot(plotting_counts) +
  geom_point(aes(x = Weevil_Location, y = normalized_counts, colour=Weevil_Location)) +
  labs(colour="Weevil Location", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~rn, scales="free", ncol=2)



                    