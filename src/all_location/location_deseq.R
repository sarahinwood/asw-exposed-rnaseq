library(data.table)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(viridis)
library(tidyverse)
library(pheatmap)

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
asw_dds_location <- copy(asw_dds)

##control for pc1, compare treatment levels
design(asw_dds_location) <- ~pc1+treatment+location
asw_dds_location <- DESeq(asw_dds_location)

#############
## results ##
#############
res_group <- results(asw_dds_location, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
##Make sig DEG data table
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
##write output
fwrite(ordered_res_group_table, "output/deseq2/asw/all_location/res_group.csv")
fwrite(sig_annots, "output/deseq2/asw/all_location/sig_w_annots.csv")
saveRDS(asw_dds_location, "output/deseq2/asw/all_location/asw_dds_location.rds")
##nr blast annots
all_annots <- merge(sig_annots, blast, by="rn", all.x=TRUE)
fwrite(all_annots, "output/deseq2/asw/all_location/all_sig_annots.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 3, pCutoff = 0.05, colAlpha=0.5,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

##plot counts for genes of interest, sub in name
plotCounts(asw_dds_location, "TRINITY_DN17206_c1_g1", intgroup = c("location"))

##########
## plot ##
##########
asw_dds_location <- readRDS("output/deseq2/asw/all_location/asw_dds_location.rds")

##genes with some annotation
blastx_trino_annots <- subset(all_annots, !is.na(all_annots$sprot_Top_BLASTX_hit))
blast_nr_annots <- subset(all_annots, !is.na(all_annots$annotation))
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
  geom_point(aes(x = Weevil_Location, y = normalized_counts, colour=Weevil_Location), alpha=0.5) +
  labs(colour="Weevil Location", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~rn, scales="free")

##low counts transcript lengths
transcript_info <- fread("data/asw-transcriptome/output/trinity_abundance/RSEM.isoforms.results")
transcript_length <- transcript_info[,c(1,3,4)]
sig_length <- merge(sig_annots, transcript_length, by="transcript_id", all.x=TRUE)
lowcount <- subset(sig_length, baseMean<1)
mean(lowcount$length)

highcount<- subset(sig_length, baseMean>=5)
mean(highcount$length)
sum(!(highcount$prot_id==""))

#############
## heatmap ##
#############

asw_dds_location <- readRDS("output/deseq2/asw/all_location/asw_dds_location.rds")
degs <- fread("output/deseq2/asw/all_location/sig_w_annots.csv")

##vst transform
vst <- varianceStabilizingTransformation(asw_dds_location, blind=TRUE)
vst_assay_dt <- data.table(assay(vst), keep.rownames=TRUE)
##subset for DEGs
vst_degs <- subset(vst_assay_dt, rn %in% degs$rn)
##turn first row back to row name
vst_degs_plot <- vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")

##get location label info
sample_to_label <- data.table(data.frame(colData(asw_dds_location)[,c("Weevil_Location", "sample_name")]))
sample_to_label <- sample_to_label %>% remove_rownames %>% column_to_rownames(var="sample_name")

location_colours <- list(Weevil_Location = c(Invermay="#3B0F70FF", Ruakura="#B63679FF"))
##plot
##not clustered by sample
pheatmap(vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_label, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))