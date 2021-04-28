library(data.table)
library(dplyr)

interaction_degs <- fread("output/deseq2/asw/location_exposure_pc1_int/sig_w_annots.csv", na.strings="")
setnames(interaction_degs, old=c("#gene_id"), new=c("rn"))
location_degs <- fread("output/deseq2/asw/location/sig_w_annots.csv", na.strings="")
exposed_degs <- fread("output/deseq2/asw/exposed/sig_annots.csv", na.strings="")

all_degs <- full_join(interaction_degs, location_degs)
all_degs <- full_join(all_degs, exposed_degs)
all_degs <- all_degs[,c(1,2,3,7)]
##blank cells with no blastx are not all set to NA
##none of the genes without BlastX have BlastP either
unann_degs <- subset(all_degs, is.na(all_degs$sprot_Top_BLASTX_hit))
list_unann_degs <- list(unique(unann_degs$rn))
fwrite(list_unann_degs, "output/deseq2/asw/unann_degs/unann_deg_list.txt")

        