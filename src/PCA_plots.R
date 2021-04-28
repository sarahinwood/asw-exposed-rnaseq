library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
asw_dds$RQN <- factor(paste(asw_dds$RQN))
asw_dds$concentration <- factor(paste(asw_dds$Concentration))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$treatment <- factor(paste(asw_dds$Treatment))

asw_vst <- varianceStabilizingTransformation(asw_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
pca_plot <- plotPCA(asw_vst, intgroup=c("location"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot, "percentVar")) 

##location
##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(pca_plot, aes(x=PC1, y=PC2, color=location))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="ASW source\nlocation")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()

pca_plot <- plotPCA(asw_vst, intgroup=c("treatment"), returnData=TRUE)
##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(pca_plot, aes(x=PC1, y=PC2, color=treatment))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Treatment")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()

##PCA with concentration and RQN
asw_pca_table <- data.table(plotPCA(asw_vst, intgroup=c("concentration"), returnData=TRUE))

asw_pca_table$concentration <- sapply(asw_pca_table$concentration, as.character)
asw_pca_table$concentration <- sapply(asw_pca_table$concentration, as.numeric)

##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(asw_pca_table, aes(x=PC1, y=PC2, color=concentration))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis()+
  labs(colour="RNA\nconcentration\n(ng/uL)")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()

##RQN
asw_pca_table <- data.table(plotPCA(asw_vst, intgroup=c("RQN"), returnData=TRUE))
asw_pca_table$RQN <- sapply(asw_pca_table$RQN, as.character)
asw_pca_table$RQN <- sapply(asw_pca_table$RQN, as.numeric)

##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(asw_pca_table, aes(x=PC1, y=PC2, color=RQN))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis()+
  labs(colour="RNA quality\nnumber")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()
