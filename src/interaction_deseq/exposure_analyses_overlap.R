library(data.table)
library(VennDiagram)
library(ggplot2)
library(DESeq2)
library(dplyr)
library(tidyr)
library(stringr)
library(viridis)

###############
## DEG lists ##
###############

exposure_ruakura <-  fread("output/deseq2/asw/interaction/exposure_ruakura/sig_w_annots.csv")
exposure_invermay <- fread("output/deseq2/asw/interaction/exposure_invermay/sig_w_annots.csv")
location_exposure <- fread("output/deseq2/asw/interaction/location_exposure/sig_w_annots.csv")
interaction <- fread("output/deseq2/asw/interaction/interaction/sig_w_annots.csv")

trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv")
blast <- fread("output/deseq2/asw/unann_degs/interaction_blastx_unann_res.csv")
trinotate_blast <- merge(trinotate, blast, by.x="#gene_id", by.y="rn", all.x=TRUE)

#############
## overlap ##
#############

vd2 <- venn.diagram(x = list("effect of exposure on Ruakura ASW"=exposure_ruakura$rn,
                             "effect of exposure on Invermay ASW"=exposure_invermay$rn,
                             "location specific repsonse to exposure"=location_exposure$rn),
                    fill=c("#440154FF", "#FDE725FF", "#21908CFF"), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd2)
ggsave(file="output/deseq2/asw/interaction/venn_diagrams/locex_interaction.svg", plot=vd2, width=6, height=4)

######################
##overlap gene lists##
######################

##location_exposure specific
a <- setdiff(location_exposure$rn, exposure_ruakura$rn)
location_exposure_sp <- setdiff(a, exposure_invermay$rn)
location_exposure_sp_annots <- subset(trinotate_blast, `#gene_id` %in% location_exposure_sp)

#exposure_ruakura specific
b <- setdiff(exposure_ruakura$rn, location_exposure$rn)
ex_ru_sp <- setdiff(b, exposure_invermay$rn)
ex_ru_sp_annots <- subset(trinotate_blast, `#gene_id` %in% ex_ru_sp)

#exposure_invermay specific
c <- setdiff(exposure_invermay$rn, location_exposure$rn)
ex_inv_sp <- setdiff(c, exposure_ruakura$rn)
ex_inv_sp_annots <- subset(trinotate_blast, `#gene_id` %in% ex_inv_sp)


##location nc vs ruakura exposed
locnc_ruex <- intersect(location_nc$rn, exposure_ruakura$rn)
locnc_ruex_annots <- subset(trinotate_blast, `#gene_id` %in% locnc_ruex)
##location nc vs invermay exposed
locnc_invex <- intersect(location_nc$rn, exposure_invermay$rn)
locnc_invex_annots <- subset(trinotate_blast, `#gene_id` %in% locnc_invex)
##ruakura vs invermay exposed
ruex_inex <- intersect(exposure_invermay$rn, exposure_ruakura$rn)
ruex_inex_annots <- subset(trinotate_blast, `#gene_id` %in% ruex_inex)
##all
all <- intersect(locnc_ruex_annots$`#gene_id`, locnc_invex_annots$`#gene_id`)
all_annots <- subset(trinotate_blast, `#gene_id` %in% all)



###############################
## plot multiple gene counts ##
###############################

asw_dds_int <- readRDS("output/deseq2/asw/interaction/dds_interaction.rds")
##get gene counts
counts_table <- data.table(counts(asw_dds_int, normalized=TRUE), keep.rownames = TRUE)

##change depending on gene list
annot_counts <- filter(counts_table, rn %in% ex_ru_sp)

##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:25], key="sample_name", value="normalized_counts")
##sample group information
sample_table <- fread("data/sample_table.csv")
sample_table$group <- paste(sample_table$Weevil_Location, sample_table$Treatment, sep=" ")
name_vs_group <- sample_table[,c(1,17)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)
plotting_counts$group <- str_replace_all(plotting_counts$group, "Exposed", "exposed")
plotting_counts$group <- str_replace_all(plotting_counts$group, "NC", "negative control")
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

