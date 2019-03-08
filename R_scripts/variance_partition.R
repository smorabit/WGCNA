library(WGCNA)
library(flashClust)
library(gplots)
library(cluster)
library(igraph);
library(RColorBrewer);
library(variancePartition)
library(doParallel)

options(stringsAsFactors=FALSE)
enableWGCNAThreads()


load("/home/vivek/AMP_AD/Mayo/Analysis/Step04_WGCNA/Discovery_Set/AD/rWGCNA_Mayo_ForPreservation.rda")

#=====================================================================================
#  Part 1: variance explained Mayo MEs
#=====================================================================================

# specify variables
# note the different format for continuous vs categorical variables
form <- ~ AgeAtDeath + RIN + Seq.PC1 + Seq.PC2 + (1|Diagnosis) + (1|Gender) + (1|Source)

# reformat metadata:
meta <- as.data.frame(cbind(targets.Ref$AgeAtDeath, targets.Ref$RIN, targets.Ref$Seq.PC1,
                      targets.Ref$Seq.PC2, targets.Ref$Diagnosis, targets.Ref$Gender,
                      targets.Ref$Source))
colnames(meta) <- c("AgeAtDeath", "RIN", "Seq.PC1", "Seq.PC2", "Diagnosis", "Gender", "Source")
meta$AgeAtDeath <- as.numeric(meta$AgeAtDeath)
meta$RIN <- as.numeric(meta$RIN)
meta$Seq.PC1 <- as.numeric(meta$Seq.PC1)
meta$Seq.PC2 <- as.numeric(meta$Seq.PC2)

varPart <- fitExtractVarPartModel(t(MEs.Mayo), form, meta)

# sort variables (i.e. columns) by median fraction of variance explained
vp <- sortCols( varPart )

# Bar plot of variance fractions for the first 10 genes
bar_plot <- plotPercentBars( vp )
ggsave("mayo_variance_fraction.pdf", bar_plot)

# violin plot of contribution of each variable to total variance
violin_plot <- plotVarPart( vp )
ggsave("mayo_variance_partition.pdf", violin_plot)

################################################################################
# Part 2: Project MAYO modules onto ROSMAP dataset:
################################################################################

# load ROSMAP metadata:
rosmap_metadata <- read.csv("/home/smorabit/data/rosmap_RNAseq_Metadata_Complete_updated.csv")


############ROSMAP Dataset
library(WGCNA)

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/ROSMAP/ROSMAP_DLPFC_expression_metaData.rda')

datExpr.ROSMAP=as.data.frame(t(normExpr.ROSMAP))

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')


setLabels=c("Mayo.AD","ROSMAP")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.ROSMAP))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.ROSMAP=datExpr.ROSMAP[,match(gnS,colnames(datExpr.ROSMAP))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), ROSMAP =list(data= datExpr.ROSMAP))
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_ROSMAP_Mayo_AD.rda"
save(list=ls(),file=filename)
