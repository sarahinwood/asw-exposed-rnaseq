library(data.table)
library(DESeq2)
library(tidyr)
library(ggplot2)
library(viridis)
library(stringr)
library(dplyr)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = "")

#################
## add factors ##
#################
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$treatment <- factor(paste(asw_dds$Treatment))
asw_dds$pc1 <- factor(paste(asw_dds$PC1_sign))
##set reference levels
asw_dds$location <- relevel(asw_dds$location, ref="Invermay")
asw_dds$treatment <- relevel(asw_dds$treatment, ref="NC")
asw_dds_int <- copy(asw_dds)
##remove genes with less than 30 counts
keep <- rowSums(counts(asw_dds_int)) >= 30
asw_dds_int <- asw_dds_int[keep,]

################
## run DESeq2 ##
################
##control for all factors and test for location-specific responses to exposure
design(asw_dds_int) <- ~pc1+location+treatment+location:treatment
##want interaction between location:treatment
##WT as all factors only have 2 levels - LRT is for 3+
asw_dds_int <- DESeq(asw_dds_int, test="Wald")
saveRDS(asw_dds_int, "output/deseq2/asw/location_exposure_pc1_int/dds_interaction.rds")

#############
## results ##
#############
resultsNames(asw_dds_int)
##first level of results should be interaction
res_group <- results(asw_dds_int, alpha=0.05)
summary(res_group)
ordered_res_group <- res_group[order(res_group$padj),]
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/deseq2/asw/location_exposure_pc1_int/res_group.csv")
##sig degs
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
sig_gene_names <- ordered_sig_res_group_table[,c(1)]
fwrite(sig_gene_names, "output/deseq2/asw/location_exposure_pc1_int/sig_deg_names.csv")
sig_annots <- merge(trinotate, ordered_sig_res_group_table, by.x="#gene_id", by.y="rn", all.y=TRUE)
fwrite(sig_annots, "output/deseq2/asw/location_exposure_pc1_int/sig_w_annots.csv")

##plot counts for genes of interest, sub in name
plotCounts(asw_dds_int, "TRINITY_DN14828_c0_g1", intgroup = c("location", "treatment"))

##blast results
unann_blast <- fread("output/deseq2/asw/unann_degs/blastx_unann_res.csv")
unann_blast <- unann_blast[,c(1,3,11,13)]
sig_blast_annots <- merge(sig_annots, unann_blast, by="transcript_id", all.x=TRUE)
fwrite(sig_blast_annots, "output/deseq2/asw/location_exposure_pc1_int/sig_blast_annots.csv")

##genes with some annotation
blastx_trino_annots <- subset(sig_blast_annots, !is.na(sig_blast_annots$sprot_Top_BLASTX_hit))
blast_nr_annots <- subset(sig_blast_annots, !is.na(sig_blast_annots$annotation))
degs_annots <- full_join(blastx_trino_annots, blast_nr_annots)                        

##get gene counts
counts_table <- data.table(counts(asw_dds_int, normalized=TRUE), keep.rownames = TRUE)
annot_counts <- filter(counts_table, rn %in% degs_annots$`#gene_id`)
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
plotting_counts$label <- c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t")
plotting_counts$gene_label <- paste(plotting_counts$label, plotting_counts$rn, sep=") ")
##plot all annot DEGs using ggplot2
ggplot(plotting_counts) +
  geom_point(aes(x = group, y = normalized_counts, colour=group)) +
  labs(colour="Sample group", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~gene_label, scales="free", ncol=3)


