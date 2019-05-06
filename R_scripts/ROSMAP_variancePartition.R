library(WGCNA)
library(flashClust)
library(gplots)
library(cluster)
library(igraph);
library(RColorBrewer);
library(variancePartition)
library(doParallel)
library(dplyr)

options(stringsAsFactors=TRUE)
enableWGCNAThreads()

# load mayo data
load("/home/vivek/AMP_AD/Mayo/Analysis/Step04_WGCNA/Discovery_Set/AD/rWGCNA_Mayo_ForPreservation.rda")

# load ROSMAP metadata and DLPFC expression data:
rosmap_metadata <- read.csv("/home/smorabit/data/rosmap_RNAseq_Metadata_Complete_updated.csv")
load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/ROSMAP/ROSMAP_DLPFC_expression_metaData.rda')
datExpr.ROSMAP=as.data.frame(t(normExpr.ROSMAP))

# intersect genes from Mayo and ROSMAP datasets
gnS <- intersect(colnames(datExpr.Ref),colnames(datExpr.ROSMAP))
datExpr.Ref1 <- datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.ROSMAP <- datExpr.ROSMAP[,match(gnS,colnames(datExpr.ROSMAP))]
moduleColors.match <- moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]

# match ROSMAP samples with ROSMAP metadata object
ids <- intersect(rosmap_metadata$SampleID, targets.ROSMAP$SampleID)
datExpr.ROSMAP <- subset(datExpr.ROSMAP, rownames(datExpr.ROSMAP) %in% ids)
rosmap_metadata <- subset(rosmap_metadata, rosmap_metadata$SampleID %in% ids)

# compute MEs for ROSMAP:
MEs.ROSMAP <- moduleEigengenes(datExpr.ROSMAP, colors = moduleColors.match)$eigengenes

# add column for pathological AD in ROSMAP rosmap metadata
rosmap_metadata$pathologic_AD <- ifelse(
  ((rosmap_metadata$BRAAK > 2) &
  (rosmap_metadata$CERAD < 3)),
  1,0)

# formula for variance partition
form <- ~ age_death + RINcontinuous + (1|BRAAK) + (1|CERAD) + (1|pathologic_AD) + (1|Batch) + (1|msex)

# select metadata for variance partition:
meta <- as.data.frame(cbind(rosmap_metadata$age_death, rosmap_metadata$RINcontinuous,
                      rosmap_metadata$BRAAK, rosmap_metadata$CERAD, rosmap_metadata$pathologic_AD,
                      rosmap_metadata$Batch, rosmap_metadata$msex))
colnames(meta) <- c("age_death", "RINcontinuous", "BRAAK", "CERAD", "pathologic_AD", "Batch", "msex")
meta$age_death <- as.numeric(meta$age_death)
meta$RINcontinuous <- as.numeric(meta$RINcontinuous)
meta$msex <- as.factor(meta$msex)
meta$BRAAK <- as.factor(meta$BRAAK)
meta$CERAD <- as.factor(meta$CERAD)
meta$Batch <- as.factor(meta$Batch)
meta$pathologic_AD <- as.factor(meta$pathologic_AD)

# run variance partition
varPart <- fitExtractVarPartModel(t(MEs.ROSMAP), form, meta)
vp <- sortCols( varPart )

# Bar plot of variance fractions for the first 10 genes
bar_plot <- plotPercentBars( vp )
ggsave("rosmap_variance_fraction.pdf", bar_plot)

# violin plot of contribution of each variable to total variance
violin_plot <- plotVarPart( vp)
ggsave("rosmap_variance_partition.pdf", violin_plot)
