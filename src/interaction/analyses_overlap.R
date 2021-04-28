library(data.table)
library(VennDiagram)

##DEG lists
exposure_ruakura <-  fread("output/deseq2/asw/interaction/exposure_ruakura/sig_w_annots.csv")
exposure_invermay <- fread("output/deseq2/asw/interaction/exposure_invermay/sig_w_annots.csv")
location_nc <- fread("output/deseq2/asw/interaction/location_nc/sig_w_annots.csv")
location_exposure <- fread("output/deseq2/asw/interaction/location_exposure/sig_w_annots.csv")
interaction <- fread("output/deseq2/asw/interaction/interaction/sig_w_annots.csv")

##venn diagrams
vd1 <- venn.diagram(x = list("exposure_ruakura"=exposure_ruakura$rn, "exposure_invermay"=exposure_invermay$rn, "location_exposure"=location_exposure$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)

vd2 <- venn.diagram(x = list("exposure_ruakura"=exposure_ruakura$rn, "exposure_invermay"=exposure_invermay$rn, 
                            "location_exposure"=location_exposure$rn,
                            "interaction"=interaction$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd2)

##plot counts for analysis specific DEGs
asw_dds_int <- readRDS("output/deseq2/asw/interaction/dds_interaction.rds")

##location_exposure specific
a <- setdiff(location_exposure$rn, exposure_ruakura$rn)
loc_ex_sp <- setdiff(a, exposure_invermay$rn)

#exposure_ruakura specific
b <- setdiff(exposure_ruakura$rn, location_exposure$rn)
ex_ru_sp <- setdiff(b, exposure_invermay$rn)

#exposure_invermay specific
c <- setdiff(exposure_invermay$rn, location_exposure$rn)
ex_inv_sp <- setdiff(c, exposure_ruakura$rn)

##overlaps
locex_ruex <- intersect(location_exposure$rn, exposure_ruakura$rn)
locex_invex <- intersect(location_exposure$rn, exposure_invermay$rn)
ruex_inex <- intersect(exposure_invermay$rn, exposure_ruakura$rn)
all <- intersect(ruex_inex, location_exposure$rn)

##interaction analysis
##specific
d <- setdiff(interaction$rn, location_exposure$rn)
e <- setdiff(d, exposure_ruakura$rn)
int_specific <- setdiff(e, exposure_invermay$rn)
##not interaction specific
shared_int <- setdiff(interaction$rn, int_specific)

###############################
## plot multiple gene counts ##
###############################

##get gene counts
counts_table <- data.table(counts(asw_dds_int, normalized=TRUE), keep.rownames = TRUE)

##change depending on gene list
annot_counts <- filter(counts_table, rn %in% shared_int)

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
  facet_wrap(~gene_label, scales="free", ncol=4)
