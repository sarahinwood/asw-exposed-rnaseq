library(DESeq2)
library(data.table)
library(ggplot2)
library(dplyr)
library(Mfuzz)
library(viridis)
library(stringr)

seed <- 6
set.seed(seed)

##function to calculate geometric mean
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

##dds object
asw_dds_int <- readRDS("output/deseq2/asw/location_exposure_pc1_int/dds_pc1sign.rds")
##list of sig DEG names
sig_gene_names <- fread("output/deseq2/asw/location_exposure_pc1_int/sig_deg_names.csv")[,unique(sig_gene_names)]
sig_gene_names_list <- list(sig_gene_names$rn)

##vst log transform data - absolute counts now changed so cannot compare between genes (only between samples for 1 gene)
vst <- varianceStabilizingTransformation(asw_dds_int, blind = FALSE)
##make matrix of transformed data
vst_matrix <- data.table(as.matrix(assay(vst)), keep.rownames = TRUE)
##melt to make long table rather than wide
long_vst_data <- melt(vst_matrix, id.vars = "rn", variable.name = "Sample_name", value.name = "vst")
##make table of colData from dds object
long_coldata <- data.table(as.data.frame(colData(asw_dds_int)))
long_coldata$group <- paste(long_coldata$Weevil_Location, long_coldata$Treatment, sep="_")

##merge sample data with vst data
merged_exp_coldata <- merge(long_vst_data, long_coldata, all.x = TRUE, all.y = FALSE, by.x="Sample_name", by.y="sample_name")
##generate table of mean vst values
mean_vst <- merged_exp_coldata[,.(vst_mean = gm_mean(vst)),by=.(rn, group)]
##make long table a wide table instead
mean_vst_wide <- dcast(mean_vst, rn~group)
##make matrix with gene names as row names
expression_matrix <- as.matrix(data.frame(mean_vst_wide, row.names = "rn"))
sig_expression_matrix <- expression_matrix[sig_gene_names_dt$rn,]

##linking treatment labels to treatments
pheno_data <- data.frame(row.names = colnames(expression_matrix), treatment = colnames(expression_matrix))
vg <- ExpressionSet(assayData = sig_expression_matrix, phenoData = new('AnnotatedDataFrame', data = pheno_data))

##standardise expression values to have mean of 0 and st dev of 1 - necessary for mfuzz
vg_s <- standardise(vg)
#optimise parameters - mestimate(vg_s)?
mestimate(vg_s)
##m determines influence of noise on cluster analysis - increasing m reduces the influence of genes with low membership values
#m prevents clustering of random data
##mestimate gave 3.29 - but don't want it this high???
m <- 2
#use to determine no. clusters - Dmin = mindist between clusters, should decline slower after reaching optimal no. of clusters
x <- Dmin(vg_s, m, crange = seq(4, 16, 1), repeats = 1)
diff(x)
##7 looks best - declines slowly afterwards

#c=no. clusters
c1 <- mfuzz(vg_s, c = 5, m=m)
clusters<- acore(vg_s, c1, min.acore = 0.7)
##centre=TRUE, centre.col="black", centre.lwd=4 - puts a thick black line through middle of each cluster showing average pattern
pdf("output/deseq2/asw/location_exposure_pc1_int/mfuzz/mfuzz_plot.pdf", height=7, width=8)
mfuzz.plot2(vg_s, c1, mfrow = c(3, 2), min.mem = 0.7, x11=FALSE, time.labels=c("Invermay Exposed", "Invermay NC", "Ruakura Exposed", "Ruakura NC"), xlab="Sample group")
dev.off()

cluster_membership <- rbindlist(clusters, idcol = "cluster")

cluster_expr_wide <- data.table(exprs(vg_s), keep.rownames = TRUE)
setnames(cluster_expr_wide, "rn", "NAME")
cluster_expr <- melt(cluster_expr_wide,
                     id.vars = "NAME",
                     variable.name = "time",
                     value.name = "scaled_vst")

##66 genes clustered
cluster_pd <- merge(cluster_membership,
                    cluster_expr,
                    by = "NAME",
                    all.x = TRUE,
                    all.y = FALSE)
fwrite(cluster_pd, "output/deseq2/asw/location_exposure_pc1_int/mfuzz/gene_clusters.csv")

##merge with annots
sig_blast_annots <- fread("output/deseq2/asw/location_exposure_pc1_int/sig_blast_annots.csv")
cluster_annotations <- merge(cluster_pd, sig_blast_annots, by.x="NAME", by.y="#gene_id", all.x=TRUE)
fwrite(cluster_annotations, "output/deseq2/asw/location_exposure_pc1_int/mfuzz/clusters_annots.csv")

##format for plot
cluster_pd$cluster_label <- paste("Cluster", cluster_pd$cluster, sep=" ")
cluster_pd$time <- str_replace_all(cluster_pd$time, "_NC", " negative control")
cluster_pd$time <- str_replace_all(cluster_pd$time, "_Exposed", " exposed")
##reorder groups
group_order <- c("Invermay negative control", "Invermay exposed", "Ruakura negative control", "Ruakura exposed")
cluster_pd$time <- factor(cluster_pd$time, levels=group_order)

##dot plot makes it more clear that not looking across time
##BUT also loses link between dots for each gene - add in geom_line
ggplot(cluster_pd, aes(x = time, y = scaled_vst, colour = `MEM.SHIP`, group = NAME)) +
  theme_bw(base_size = 8) +
  xlab("Sample group") + ylab("Scaled, mapped reads") +
  facet_wrap(~ cluster_label, scales="fixed") +
  geom_line(color="grey", alpha=0.3, size=0.5)+
  geom_point()+
  scale_color_viridis(name="Cluster\nmembership")
