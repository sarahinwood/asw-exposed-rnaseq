library(data.table)
library(dplyr)
library(VennDiagram)

###############
## DEG lists ##
###############

location <- fread("output/deseq2/asw/all_location/sig_w_annots.csv")
exposure <- fread("output/deseq2/asw/all_exposed/sig_w_annots.csv")

vd1 <- venn.diagram(x = list("Location DEGs"=location$rn,
                             "Exposure DEGs"=exposure$rn),
                    fill=c("#440154FF", "#FDE725FF"), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)
ggsave(file="output/deseq2/asw/venn_diagrams/location_exposure_venn.svg", plot=vd1, width=6, height=2)

ex_sp <- setdiff(exposure$rn, location$rn)
ex_sp_annots <- subset(exposure, rn %in% ex_sp)

loc_sp <- setdiff(location$rn, exposure$rn)
loc_sp_annots <- subset(location, rn %in% loc_sp)

intersect <- intersect(exposure$rn, location$rn)
intersect_annots <- subset(location, rn %in% intersect)

##table of counts for intersecting DEGs
asw_dds_location <- readRDS("output/deseq2/asw/all_exposed/asw_dds_exposed.rds")
counts_table <- data.table(counts(asw_dds_location, normalized=TRUE), keep.rownames = TRUE)
intersect_counts <- filter(counts_table, rn %in% intersect)
fwrite(intersect_counts, "output/deseq2/asw/venn_diagrams/loc_ex_intersect_counts.csv")
