library(data.table)
library(VennDiagram)

dual_location <- fread("output/deseq2/asw_dual/location_pairwise/sig_w_annots.csv")
dual_location$edited_rn <- tstrsplit(dual_location$rn, "ASW_", keep=c(2))
dual_location_5 <- subset(dual_location, padj<0.05)
location <- fread("output/deseq2/asw/location_pairwise/sig_w_annots.csv")
location_5 <- subset(location, padj<0.05)

##location overlap padj<0.1
vd <- venn.diagram(x = list("Location"=location$`#gene_id`, "Dual Location"=dual_location$edited_rn), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1,)
grid.newpage()
grid.draw(vd)
##location overlap padj<0.05
vd <- venn.diagram(x = list("Location padj<0.05"=location_5$`#gene_id`, "Dual Location padj<0.05"=dual_location_5$edited_rn), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1,)
grid.newpage()
grid.draw(vd)

dual_int <- fread("output/deseq2/asw_dual/location_exposure_int/sig_w_annots.csv")
dual_int$edited_rn <- tstrsplit(dual_int$rn, "ASW_", keep=c(2))
dual_int_5 <- subset(dual_int, padj<0.05)
int <- fread("output/deseq2/asw/location_exposure_int/sig_w_annots.csv")
int_5 <- subset(int, padj<0.05)

##interaction overlap padj<0.1
vd2 <- venn.diagram(x = list("int"=int$rn, "Dual int"=dual_int$edited_rn), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1,)
grid.newpage()
grid.draw(vd2)
##interaction overlap padj<0.05
vd2 <- venn.diagram(x = list("int"=int_5$rn, "Dual int"=dual_int_5$edited_rn), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1,)
grid.newpage()
grid.draw(vd2)

##dual location vs dual interaction:
vd2 <- venn.diagram(x = list("dual int"=dual_int$edited_rn, "dual location"=dual_location$edited_rn), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1,)
grid.newpage()
grid.draw(vd2)

##location vs interaction:
vd2 <- venn.diagram(x = list("int"=int$rn, "location"=location$`#gene_id`), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1,)
grid.newpage()
grid.draw(vd2)
