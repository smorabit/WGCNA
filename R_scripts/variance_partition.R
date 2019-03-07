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
#  Part 1: variance explained
#=====================================================================================

cl <- makeCluster(4)
registerDoParallel(cl)

# specify variables
# note the different format for continuous vs categorical variables
form <- ~ AgeAtDeath + RIN + (1|Diagnosis) + (1|Gender)
# data(varPartData)

# reformat metadata:
meta <- as.data.frame(cbind(targets.Ref$AgeAtDeath, targets.Ref$RIN, targets.Ref$Diagnosis,
                      targets.Ref$Gender))
colnames(meta) <- c("AgeAtDeath", "RIN", "Diagnosis", "Gender")
meta$AgeAtDeath <- as.numeric(meta$AgeAtDeath)
meta$RIN <- as.numeric(meta$RIN)

varPart <- fitExtractVarPartModel(t(datExpr.Ref), form, meta)

# testing example
data(varPartData)
form <- ~ Age + (1|Individual) + (1|Tissue) + (1|Batch)
varPart <- fitExtractVarPartModel( geneExpr, form, info )
