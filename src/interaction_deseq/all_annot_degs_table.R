library(data.table)
library(dplyr)
library(VennDiagram)
library(viridis)

###############
## DEG lists ##
###############

exposure_ruakura <-  fread("output/deseq2/asw/interaction/exposure_ruakura/all_sig_annots.csv", na.strings="")
exposure_ruakura$analysis <- paste("exposure_ruakura")

exposure_invermay <- fread("output/deseq2/asw/interaction/exposure_invermay/all_sig_annots.csv", na.strings="")
exposure_invermay$analysis <- paste("exposure_invermay")

location_exposure <- fread("output/deseq2/asw/interaction/location_exposure/all_sig_annots.csv", na.strings="")
location_exposure$analysis <- paste("location_exposure")

interaction <- fread("output/deseq2/asw/interaction/interaction/all_sig_annots.csv", na.strings="")
interaction$analysis <- paste("interaction")

location_nc <- fread("output/deseq2/asw/interaction/location_nc/all_sig_annots.csv", na.strings="")
location_nc$analysis <- paste("location_nc")

a <- full_join(exposure_ruakura, exposure_invermay)
b <- full_join(a, location_exposure)
c <- full_join(b, interaction)
full_degs_table <- full_join(c, location_nc)
length(unique(full_degs_table$rn))

##returns other DEGs with no annot but not set to NA
trin_blastx <- subset(full_degs_table, !is.na(full_degs_table$sprot_Top_BLASTX_hit))
##returns other DEGs with no annot but not set to NA
trin_blastp <- subset(full_degs_table, !is.na(full_degs_table$sprot_Top_BLASTP_hit))
##is working
nr_blastx <- subset(full_degs_table, !is.na(full_degs_table$annotation))

##table of all DEGs with annots
trin_annots <- full_join(trin_blastx, trin_blastp)
all_annots <- full_join(trin_annots, nr_blastx)
length(unique(all_annots$rn))
all_annots_table <- all_annots[,c(1,9,13,36,37)]
fwrite(all_annots_table, "output/deseq2/asw/interaction/all_degs_tables/all_annot_degs.csv")

##venn diagram of annotated DEGs
ex_ru <- subset(all_annots_table, analysis=="exposure_ruakura")
ex_inv <- subset(all_annots_table, analysis=="exposure_invermay")
loc_nc <- subset(all_annots_table, analysis=="location_nc")
loc_ex <- subset(all_annots_table, analysis=="location_exposure")
int <- subset(all_annots_table, analysis=="interaction")

vd <- venn.diagram(x = list("effect of exposure on Ruakura ASW"=ex_ru$rn,
                             "effect of exposure on Invermay ASW"=ex_inv$rn,
                             "effect of location on unexposed ASW"=loc_nc$rn,
                             "location specific repsonse to exposure"=loc_ex$rn,
                             "location-exposure interaction"=int$rn),
                    fill=c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd)
