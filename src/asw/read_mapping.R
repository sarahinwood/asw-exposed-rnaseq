library(data.table)
library(DESeq2)

##asw reads mapped
dds_group_asw <- readRDS("output/deseq2/asw/asw_dds.rds")
counts_table_asw <- (data.table(counts(dds_group_asw)))
counts_colSums_asw <- setDT(data.frame(colSums(counts_table_asw, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_asw, old=c("rn", "colSums.counts_table_asw..na.rm...TRUE."), new=c("Sample_name", "readpairs_mapped_ASW"))

##both reads mapped
bbduk_reads_out <- fread("output/bbduk_trim/bbduk_reads_out.csv")
full_read_mapping <- merge(counts_colSums_asw, bbduk_reads_out, by="Sample_name")
full_read_mapping$bbduk_halved <- (full_read_mapping$bbduk_reads_out)/2
full_read_mapping$`%_ofall_ASW` <- (full_read_mapping$readpairs_mapped_ASW/full_read_mapping$bbduk_halved)*100
##slight differences from salmon - due to salmon filtering out some extra reads first
fwrite(full_read_mapping, "output/deseq2/read_mapping.csv")
