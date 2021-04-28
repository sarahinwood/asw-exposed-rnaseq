library(data.table)
library(dplyr)

exposure_ruakura <-  fread("output/deseq2/asw/interaction/exposure_ruakura/sig_w_annots.csv", na.strings = "")
exposure_invermay <- fread("output/deseq2/asw/interaction/exposure_invermay/sig_w_annots.csv", na.strings = "")
location_nc <- fread("output/deseq2/asw/interaction/location_nc/sig_w_annots.csv", na.strings = "")
location_exposure <- fread("output/deseq2/asw/interaction/location_exposure/sig_w_annots.csv", na.strings = "")
interaction <- fread("output/deseq2/asw/interaction/interaction/sig_w_annots.csv", na.strings = "")

all_degs <- full_join(exposure_ruakura, exposure_invermay)
all_degs <- full_join(all_degs, location_nc)
all_degs <- full_join(all_degs, location_exposure)
all_degs <- full_join(all_degs, interaction)
all_degs <- all_degs[,c(1,8,9,13)]
##blank cells with no blastx are not all set to NA
##none of the genes without BlastX have BlastP either
unann_degs <- subset(all_degs, is.na(all_degs$sprot_Top_BLASTX_hit))
list_unann_degs <- list(unique(unann_degs$rn))
fwrite(list_unann_degs, "output/deseq2/asw/unann_degs/interaction_unann_deg_list.txt")

