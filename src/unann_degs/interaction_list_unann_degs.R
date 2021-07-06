library(data.table)
library(dplyr)

exposure_ruakura <-  fread("output/deseq2/asw/interaction/exposure_ruakura/sig_w_annots.csv", na.strings = "")
exposure_invermay <- fread("output/deseq2/asw/interaction/exposure_invermay/sig_w_annots.csv", na.strings = "")
location <- fread("output/deseq2/asw/location/sig_w_annots.csv", na.strings = "")
location_exposure <- fread("output/deseq2/asw/interaction/location_exposure/sig_w_annots.csv", na.strings = "")
interaction <- fread("output/deseq2/asw/interaction/interaction/sig_w_annots.csv", na.strings = "")
exposure <- fread("output/deseq2/asw/all_exposed/sig_w_annots.csv", na.strings = "")

all_degs <- full_join(exposure_ruakura, exposure_invermay)
all_degs <- full_join(all_degs, location)
all_degs <- full_join(all_degs, location_exposure)
all_degs <- full_join(all_degs, interaction)
all_degs <- full_join(all_degs, exposure)

all_degs <- all_degs[,c(1,8,9,13)]

##blank cells with no blastx are not all set to NA
##none of the genes without BlastX have BlastP either
unann_degs <- subset(all_degs, is.na(all_degs$sprot_Top_BLASTX_hit))
##add TRINITY_DN514_c0_g1 - gets plant cell wall annot
DN514 <- subset(all_degs, all_degs$rn=="TRINITY_DN514_c0_g1")
all_unann_degs <- full_join(unann_degs, DN514)
list_unann_degs <- list(unique(all_unann_degs$rn))
fwrite(list_unann_degs, "output/deseq2/asw/unann_degs/interaction_unann_deg_list.txt")

