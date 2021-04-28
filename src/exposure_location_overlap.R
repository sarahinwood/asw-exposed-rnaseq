library(data.table)
library(VennDiagram)

location_degs <- fread("output/deseq2/asw/location/sig_w_annots.csv")
exposed_degs <- fread("output/deseq2/asw/exposed/sig_annots.csv")

vd <- venn.diagram(x = list("Location"=location_degs$rn, "Parasitoid\nexposure"=exposed_degs$rn), filename=NULL,
                   fill=c("#440154FF", "#FDE725FF"), alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd)

shared_degs <- data.table(intersect(location_degs$rn, exposed_degs$rn))
loc_shared <- merge(shared_degs, location_degs, by.x="V1", by.y="rn")
loc_ex_shared <- merge(loc_shared, exposed_degs, by.x="V1", by.y="rn")
shared_degs_table <- loc_ex_shared[,c(1, 2, 3, 7, 25, 29)]
fwrite(shared_degs_table, "output/deseq2/asw/location-exposed-overlap/shared_degs_table.csv")

##location specific DEGs
location_specific <- data.table(setdiff(location_degs$rn, exposed_degs$rn))
location_specific_full <- merge(location_degs, location_specific, by.x="rn", by.y="V1")
fwrite(location_specific_full, "output/deseq2/asw/location-exposed-overlap/location_specific_degs.csv")

##exposed specific DEGs
exposed_specific <- data.table(setdiff(exposed_degs$rn, location_degs$rn))
exposed_specific_full <- merge(exposed_degs, exposed_specific, by.x="rn", by.y="V1")
fwrite(exposed_specific_full, "output/deseq2/asw/location-exposed-overlap/exposed_specific_degs.csv")

asw_dds_location <- readRDS("output/deseq2/asw/location/asw_dds_location.rds")
asw_dds_exposed <- readRDS("output/deseq2/asw/exposed/asw_dds_exposed.rds")

plotCounts(asw_dds_location, "TRINITY_DN24317_c4_g1", intgroup = c("location"), main="")
plotCounts(asw_dds_exposed, "TRINITY_DN28910_c0_g1", intgroup = c("treatment"), main="")

##########
## plot ##
##########

asw_dds_location <- readRDS("output/deseq2/asw/location/asw_dds_location.rds")

##get gene counts
counts_table <- data.table(counts(asw_dds_location, normalized=TRUE), keep.rownames = TRUE)
annot_counts <- filter(counts_table, rn %in% shared_degs$V1)
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
  facet_wrap(~rn, scales="free", ncol=3)

##samples with counts while rest have zero
outlier_samples <- subset(plotting_counts, normalized_counts > 0)
fwrite(outlier_samples, "output/deseq2/asw/location-exposed-overlap/outlier_samples.csv")
length(unique(outlier_samples$sample_name))
