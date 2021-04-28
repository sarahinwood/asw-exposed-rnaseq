library("data.table")
library("dplyr")

unann_degs_blast <- fread("output/deseq2/asw/unann_degs/nr_blastx.outfmt6")
setnames(unann_degs_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("transcript_id", "blast_hit", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##remove uncharacterized
fil_blast <- data.table(dplyr::filter(unann_degs_blast, !grepl('uncharacterized', annotation, ignore.case=TRUE)))
fil_blast <- data.table(dplyr::filter(fil_blast, !grepl('hypothetical', annotation, ignore.case=TRUE)))
setorder(fil_blast, transcript_id, evalue, -bit_score)
blast_res <- fil_blast[,.SD[which.min(evalue)], by=transcript_id]

fwrite(blast_res, "output/deseq2/asw/unann_degs/blastx_unann_res.csv")
