library(data.table)
library(VennDiagram)

##DEG lists
location_int <- fread("output/deseq2/asw/interaction/location_nc/sig_w_annots.csv")
location_old <- fread("output/deseq2/asw/location/sig_w_annots.csv")

##venn diagram
vd1 <- venn.diagram(x = list("location_int"=location_int$rn, "location_old"=location_old$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)

##location_exposure specific
overlapped <- intersect(location_int$rn, location_old$rn)
interaction_sp <- setdiff(location_int$rn, location_old$rn)
old_sp <- setdiff(location_old$rn, location_int$rn)

###############################
## plot multiple gene counts ##
###############################
##plot counts for analysis specific DEGs
asw_dds_int <- readRDS("output/deseq2/asw/interaction/dds_interaction.rds")

##get gene counts
counts_table <- data.table(counts(asw_dds_int, normalized=TRUE), keep.rownames = TRUE)

##change depending on gene list
annot_counts <- filter(counts_table, rn %in% old_sp)

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
