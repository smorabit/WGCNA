##################Consensus Network WGCNA####################
# Consensus performed between human brain RNA-seq expression data from three independent consortia,
# Mayo clinic, Mt. Sinai, and ROSMAP.
setwd("~/WGCNA")
getwd()
library(WGCNA)
library(flashClust)
library(gplots)
library(cluster)
library(igraph); #for part 8
library(RColorBrewer); #for part 8
options(stringsAsFactors=FALSE)
#Only when not on server, comment out:
enableWGCNAThreads()
#=====================================================================================
#
#  Part 1: Data Input, cleaning, and pre-processing
#  Note: The only thing that I am doing so far is loading the data and then
#        moving on to step 2.
#
#=====================================================================================

#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO
gset <- getGEO("GSE59630", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL5175", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("00000000000000000000000000000000000000000000000000",
               "00000000111111111111111111111111111111111111111111",
               "1111111111111111")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
#exprs(gset) is expression matrix
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_LIST","SPOT_ID","RANGE_GB","RANGE_STRAND","RANGE_START"))
write.table(tT, file="DE_GEO.txt", row.names=F, sep="\t")

##########################Clean up data######################################
targets.ALL = pData(gset);
datExpr = as.data.frame(exprs(gset));

#remove hippocampus regions - 6 total so 110 patients left
targets.ALL.select  = targets.ALL[!(targets.ALL$region == "HIP"),]
#remove extra columns
targets.ALL.select = as.data.frame(targets.ALL.select)
targets.ALL.select2 = targets.ALL.select[ , which(names(targets.ALL.select) %in% c("geo_accession","disease status:ch1",
                                                              "age:ch1", "region:ch1", "Sex:ch1"))]
#remove hippocampus regions from expression
targets.ALL.HIP  = targets.ALL[(targets.ALL$region == "HIP"),]
GSM = targets.ALL.HIP$geo_accession;
datExpr.select = datExpr[ , -which(names(datExpr) %in% c(GSM))]

#gene IDs are affymetrix, convert to ensembl:
annot=read.csv("/home/erebboah/WGCNA/ensembl_GEO.csv") #could not do this cause I don't have this file fml
expr = as.integer(rownames(datExpr.select));
ensembl = annot$From;
ind=intersect(ensembl, expr)
#16597 x 110
datExpr.select.subset=datExpr.select[na.omit(match(ind,rownames(datExpr.select))),]
annot1= annot[na.omit(match(ind,annot$From)),]
#16077 x 110
datExpr_ENSGID=collapseRows(datExpr.select.subset, annot1$To, as.character(annot1$From))$datETcollapsed  # changing to ENSG_IDs using collapseRows
save(datExpr_ENSGID, targets.ALL.select2, file="datExpr_GSE59630.rda")

#########################################################################
#Mayo Clinic data:
load(file="~/data/MAYO_TCX_Expression_metaData.rda")
datExpr.Ref.MAYO <- t(normExpr.TCX)
targets.Ref.MAYO <- targets.TCX

#Mt. Sinai data:
load(file = "~/data/MSSM_FP_STG_PHG_IFG_Expression_metaData.rda")
datExpr.Ref.MSSM <- t(normExpr.MSSM)
targets.Ref.MSSM <- targets.MSSM

#ROSMAP data:
load(file = "~/data/ROSMAP_DLPFC_expression_metaData.rda")
datExpr.Ref.ROSMAP <- t(normExpr.ROSMAP)
targets.Ref.ROSMAP <- targets.ROSMAP

#only interested in genes that intersect across our datasets (14425 genes):
gnS <- Reduce(intersect, list(colnames(datExpr.Ref.MAYO), colnames(datExpr.Ref.MSSM), colnames(datExpr.Ref.ROSMAP)))
datExpr.Ref.MAYO <- datExpr.Ref.MAYO[,match(gnS,colnames(datExpr.Ref.MAYO))]
datExpr.Ref.MSSM <- datExpr.Ref.MSSM[,match(gnS,colnames(datExpr.Ref.MSSM))]
datExpr.Ref.ROSMAP <- datExpr.Ref.ROSMAP[,match(gnS,colnames(datExpr.Ref.ROSMAP))]

#make a list of our 3 datasets:
nSets <- 3
setLabels <- c("Mayo", "Mt. Sinai", "ROSMAP")
multiExpr <- vector(mode="list", length=length(setLabels))

i <- 1
for (dat in list(datExpr.Ref.MAYO, datExpr.Ref.MSSM, datExpr.Ref.ROSMAP)){
  multiExpr[[i]] <- list(data=dat)
  i <- i+1
}

#check that data is in the correct format using WCGNA checkSets() function:
exprSize <- checkSets(multiExpr)

#$nGenes
#[1] 14426
#$nSamples
#[1] 262 753 632

#Make multiMeta of traits for all sets - list of lists
multiMeta <- list(MAYO = list(data=targets.Ref.MAYO), MSSM = list(data=targets.Ref.MSSM), ROSMAP = list(data=targets.Ref.ROSMAP))

#Check that all genes and samples have sufficiently low numbers of missing values
gsg <- goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK #this was FALSE

#If NOT all ok -
if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

#update data based on multiExpr:
datExpr.Ref.MAYO <- multiExpr[[1]]$data
datExpr.Ref.MSSM <- multiExpr[[2]]$data
datExpr.Ref.MSSM <- multiExpr[[3]]$data

#Save everything
nGenes <- exprSize$nGenes;
nSamples <- exprSize$nSamples;
save(multiExpr, multiMeta, nGenes, nSamples, setLabels,
     targets.Ref.MAYO, targets.Ref.MSSM, targets.Ref.ROSMAP,
     datExpr.Ref.MAYO, datExpr.Ref.MSSM, datExpr.Ref.ROSMAP, file = "wgcna_consensus_01.10.19.rda");

#===============================================================================
#
#  Part 2: Choose soft thresholding power
#
#===============================================================================
# Load the data saved in the first part
load(file = "wgcna_consensus_01.10.19.rda");

# Get the number of sets in the multiExpr structure.
nSets <- checkSets(multiExpr)$nSets

# Choose a set of soft-thresholding powers
powers <- c(seq(1,10,by=1), seq(12,30, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables <- vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] <- list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 5, networkType="signed", corFnc="bicor")[[2]]);

#save the powerTables object:
save(powerTables, file="powerTables.rda")

# Plot the results
pdf("1_Power__01.10.19.pdf", height=10, width=18)
colors = c("blue", "red", "black")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()

#=====================================================================================
#
#  Part 3: Network Construction
#
#=====================================================================================
#Set soft thresholding power to number based on plots
softpower=12;
# Auto network
load(file = "wgcna_consensus_01.10.19.rda");
net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                              maxBlockSize = 30000, ## This should be set to a smaller size if the user has limited RAM
                              randomSeed = 12345,
                              corType = "pearson", ## no use for bicor
                              power = softpower,
                              consensusQuantile = 0.2,
                              networkType = "signed",
                              TOMType = "unsigned",
                              TOMDenom = "min",
                              scaleTOMs = TRUE, scaleQuantile = 0,
                              sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                              useDiskCache = TRUE, chunkSize = NULL,
                              deepSplit = 2,
                              detectCutHeight = 0.99999999, minModuleSize = 20,
                              mergeCutHeight = 0.2,
                              saveConsensusTOMs = TRUE,
                              consensusTOMFilePattern = "ConsensusTOM-block.%b.rda")

save(list=ls(),file="wgcna_consensus_network_01.11.19.rda")

#load data
load("wgcna_consensus_network_01.11.19.rda")

consMEs = net$multiMEs;
moduleColors = net$colors;
table(moduleColors)

load("ConsensusTOM-block.1.rda") # consensus TOM
consTree= hclust(1-consTomDS,method="average");

# Relate modules with traits for each dataset
# (can't be done as a loop b/c field names are not the same)
################################## MAYO ################################################

Diagnosis <- as.numeric(relevel(factor(as.character(targets.Ref.MAYO$Diagnosis)), 'CONTROL'))
Age <- as.numeric(targets.Ref.MAYO$AgeAtDeath)
RIN <- as.numeric(targets.Ref.MAYO$RIN)
Gender <- as.numeric(factor(targets.Ref.MAYO$Gender))
PCT_PF_READS_ALIGNED <- as.numeric(targets.Ref.MAYO$PCT_PF_READS_ALIGNED)

#this for loop generates a bunch of warnings
geneSigsMAYO <- matrix(NA,nrow=5,ncol=ncol(datExpr.Ref.MAYO))
for(i in 1:ncol(geneSigsMAYO)) {
  exprvec=as.numeric(datExpr.Ref.MAYO[,i])
  ager=bicor(Age,exprvec,use="pairwise.complete.obs")
  sexr=bicor(exprvec, Gender,use="pairwise.complete.obs")
  conditionr=bicor(exprvec, Diagnosis,use="pairwise.complete.obs")
  rinr=bicor(RIN,exprvec,use="pairwise.complete.obs")
  pctr=bicor(PCT_PF_READS_ALIGNED, exprvec, use="pairwise.complete.obs")
  geneSigsMAYO[,i]=c(ager, sexr, conditionr, rinr, pctr)
  cat('Done for gene...',i,'\n')
}

#set colors
for (i in 1:nrow(geneSigsMAYO)){
  geneSigsMAYO[i,] <- numbers2colors(as.numeric(geneSigsMAYO[i,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
}

rownames(geneSigsMAYO) <- c("Mayo Age","Mayo Gender","Mayo Diagnosis","Mayo RIN","Mayo Reads Aligned")

################################## MSSM ################################################
Diagnosis <- as.numeric(relevel(factor(as.character(targets.Ref.MSSM$Diagnosis)), 'CONTROL'))
Age <- as.numeric(targets.Ref.MSSM$AOD)
RIN <- as.numeric(targets.Ref.MSSM$RIN)
Gender <- as.numeric(factor(targets.Ref.MSSM$SEX))
PCT_PF_READS_ALIGNED <- as.numeric(targets.Ref.MSSM$PCT_PF_READS_ALIGNED)
CDR <- as.numeric(targets.Ref.MSSM$CDR)

#this for loop also generates a bunch of warnings
geneSigsMSSM <- matrix(NA,nrow=6,ncol=ncol(datExpr.Ref.MSSM))
for(i in 1:ncol(geneSigsMSSM)) {
  exprvec=as.numeric(datExpr.Ref.MSSM[,i])
  ager=bicor(Age,exprvec,use="pairwise.complete.obs")
  sexr=bicor(exprvec, Gender,use="pairwise.complete.obs")
  conditionr=bicor(exprvec, Diagnosis,use="pairwise.complete.obs")
  rinr=bicor(exprvec, RIN, use="pairwise.complete.obs")
  pctr=bicor(exprvec, PCT_PF_READS_ALIGNED, use="pairwise.complete.obs")
  cdrr=bicor(exprvec, CDR, use="pairwise.complete.obs")
  geneSigsMSSM[,i]=c(ager, sexr,conditionr,rinr,pctr,cdrr)
  cat('Done for gene...',i,'\n')
}

for (i in 1:nrow(geneSigsMSSM)){
  geneSigsMSSM[i,] <- numbers2colors(as.numeric(geneSigsMSSM[i,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
}
rownames(geneSigsMSSM) <- c("MSSM Age","MSSM Gender","MSSM Diagnosis", "MSSM RIN", "MSSM Reads Aligned", "MSSM CDR")

############################### ROSMAP ########################################

Diagnosis <- as.numeric(relevel(factor(as.character(targets.Ref.ROSMAP$Diagnosis)), 'CONTROL'))
Age <- as.numeric(targets.Ref.ROSMAP$age_death)
Gender <- as.numeric(factor(targets.Ref.ROSMAP$msex))
RIN <- as.numeric(targets.Ref.ROSMAP$RINcontinuous)
cogdx <- as.numeric(targets.Ref.ROSMAP$cogdx)
PCT_PF_READS_ALIGNED <- as.numeric(targets.Ref.ROSMAP$PCT_PF_READS_ALIGNED)

geneSigsROSMAP <- matrix(NA,nrow=6,ncol=ncol(datExpr.Ref.ROSMAP))
for(i in 1:ncol(geneSigsROSMAP)) {
  exprvec=as.numeric(datExpr.Ref.ROSMAP[,i])
  ager=bicor(Age,exprvec,use="pairwise.complete.obs")
  sexr=bicor(exprvec, Gender,use="pairwise.complete.obs")
  conditionr=bicor(exprvec, Diagnosis,use="pairwise.complete.obs")
  rinr=bicor(exprvec, RIN, use="pairwise.complete.obs")
  cogdxr=bicor(exprvec, cogdx, use="pairwise.complete.obs")
  pctr=bicor(exprvec, PCT_PF_READS_ALIGNED, use="pairwise.complete.obs")
  geneSigsROSMAP[,i]=c(ager, sexr, conditionr, rinr, cogdxr, pctr)
  cat('Done for gene...',i,'\n')
}

for (i in 1:nrow(geneSigsROSMAP)){
  geneSigsROSMAP[i,] <- numbers2colors(as.numeric(geneSigsROSMAP[i,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
}
rownames(geneSigsROSMAP)=c("ROSMAP Age","ROSMAP Gender","ROSMAP Diagnosis", "ROSMAP RIN", "ROSMAP cogdx", "ROSMAP Reads Aligned")

save(list=ls(),file="wgcna_consensus_network_01.11.19.rda")

############# Calculate modules for each set of parameters #####################
mColorh <- mLabelh <- colorLabels <- NULL
for (minModSize in c(40,100,160)) {
  for (dthresh in c(0.1,0.2,0.25)) {
    for (ds in c(2,4)) {
      print("Trying parameters:")
      print(c(minModSize,dthresh,ds))
      tree = cutreeHybrid(dendro = consTree, pamStage=FALSE,
                          minClusterSize = minModSize, cutHeight = 0.99999999,
                          deepSplit = ds, distM = as.matrix(1-consTomDS))

      merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels,
                                  cutHeight = dthresh)
      mColorh <- cbind(mColorh,labels2colors(merged$colors))
      mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
    }
  }
}

# Plotting modules for each set of params and traits
mColorh1=cbind(mColorh,geneSigsMAYO[1,],geneSigsMAYO[2,],geneSigsMAYO[3,],geneSigsMAYO[4,],geneSigsMAYO[5,],
               geneSigsMSSM[1,], geneSigsMSSM[2,], geneSigsMSSM[3,], geneSigsMSSM[4,], geneSigsMSSM[5,], geneSigsMSSM[6,],
               geneSigsROSMAP[1,], geneSigsROSMAP[2,], geneSigsROSMAP[3,], geneSigsROSMAP[4,], geneSigsROSMAP[5,], geneSigsROSMAP[6,])
rownames_geneSigs = c(rownames(geneSigsMAYO), rownames(geneSigsMSSM), rownames(geneSigsROSMAP))
mLabelh1=c(mLabelh,rownames_geneSigs)

pdf("ConsensusTOM_MultiDendro_DSAD_01.14.19.pdf",height=25,width=20)
plotDendroAndColors(consTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed consensus network with power = 12"));
dev.off()

##############Choose parameters for final dendrogram#########################
load("wgcna_consensus_network_01.11.19.rda");

#select dendrogram parameters
mms <- 100
ds <- 4
dthresh <- 0.2

tree <- cutreeHybrid(dendro = consTree, pamStage=FALSE,
                    minClusterSize = mms, cutHeight = 0.99999999,
                    deepSplit = ds, distM = as.matrix(1-consTomDS))

merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels,
                            cutHeight = dthresh)

# Eigengenes of the new merged modules:
MEs_MAYO <- merged$newMEs[[1]]$data;
MEs_MSSM <- merged$newMEs[[2]]$data;
MEs_ROSMAP <- merged$newMEs[[3]]$data;

# Module color associated with each gene
moduleColor.cons <- labels2colors(merged$colors)

mColorh <- cbind(labels2colors(merged$colors))
mLabelh <- c("Merged Colors")

mColorh1 <- cbind(mColorh,geneSigsMAYO[1,],geneSigsMAYO[2,],geneSigsMAYO[3,],geneSigsMAYO[4,],geneSigsMAYO[5,],
               geneSigsMSSM[1,], geneSigsMSSM[2,], geneSigsMSSM[3,], geneSigsMSSM[4,], geneSigsMSSM[5,], geneSigsMSSM[6,],
               geneSigsROSMAP[1,], geneSigsROSMAP[2,], geneSigsROSMAP[3,], geneSigsROSMAP[4,], geneSigsROSMAP[5,], geneSigsROSMAP[6,])
rownames_geneSigs <- c(rownames(geneSigsMAYO), rownames(geneSigsMSSM), rownames(geneSigsROSMAP))
mLabelh1 <- c(mLabelh,rownames_geneSigs)

pdf("ConsensusTOM_FinalDendro_MAYO_MSSM_ROSMAP_01.15.19.pdf",height=10,width=16)
plotDendroAndColors(consTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power = 12, mms=",mms,"ds=",ds,"dthresh=",dthresh,"cquant=0.2"));
dev.off()

#Just look at diagnosis traits
mColorh1 <- cbind(mColorh,geneSigsMAYO[3,], geneSigsMSSM[3,], geneSigsROSMAP[3,])
rownames_geneSigs <- c("Mayo Diagnosis", "MSSM Diagnosis", "ROSMAP Diagnosis")
mLabelh1 <- c(mLabelh,rownames_geneSigs)

pdf("ConsensusTOM_FinalDendro_MAYO_MSSM_ROSMAP_diagnosis_01.15.19.pdf",height=10,width=16)
plotDendroAndColors(consTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power = 12, mms=",mms,"ds=",ds,"dthresh=",dthresh,"cquant=0.2"));
dev.off()

# Convert numerical lables to colors for labeling of modules
MEColors <- labels2colors(as.numeric(substring(names(MEs_MAYO), 3)));
MEColorNames <- paste("ME", MEColors, sep="");

colnames(MEs_MAYO)=c(MEColorNames)
colnames(MEs_MSSM)=c(MEColorNames)
colnames(MEs_ROSMAP)=c(MEColorNames)

#Calculate kMEs and p value, Z score of kMEs for each gene (connectivity of each gene to a module eigengene) across datasets
consKME1=consensusKME(multiExpr=multiExpr, moduleColor.cons,
                      multiEigengenes = NULL,
                      consensusQuantile = 0.2,
                      signed = TRUE)
#Consensus kMEs:
consensus.KMEs=consKME1[,regexpr('consensus.kME',names(consKME1))>0]

save(consensus.KMEs, consKME1, multiExpr, MEColorNames, moduleColor.cons, consTree,
     MEs_MAYO, MEs_MSSM, MEs_ROSMAP, file="Consensus_MEs_01.15.19.rda")

#=====================================================================================
#
#  Part 4: Get annotated gene lists and their associated kMEs
#
#=====================================================================================
load("Consensus_MEs_01.15.19.rda")
ensembl <- read.csv("/home/vivek/FTD_Seeley/Analysis_Nov2017/ENSG85_Human.csv.gz") # Convert Ensembl gene ID to gene names
consensus.KMEs$Ensembl.Gene.ID <- paste(rownames(consensus.KMEs))

merged <- merge(consensus.KMEs,ensembl,by.x="Ensembl.Gene.ID",by.y="Ensembl.Gene.ID",all.x=T)
ind <- match(consensus.KMEs$Ensembl.Gene.ID,merged$Ensembl.Gene.ID)
merged1 <- merged[ind,]
consensus.KMEs.annot <- merged1

geneInfo.cons <- as.data.frame(cbind(consensus.KMEs.annot$Ensembl.Gene.ID,consensus.KMEs.annot$Associated.Gene.Name,
                                  moduleColor.cons,consensus.KMEs))
geneInfo.cons <- geneInfo.cons[,-ncol(geneInfo.cons)] # check if last column is Ensembl gene id

colnames(geneInfo.cons)[1]= "Ensembl.Gene.ID"
colnames(geneInfo.cons)[2]= "GeneSymbol"
colnames(geneInfo.cons)[3]= "Initially.Assigned.Module.Color"

write.csv(geneInfo.cons,'geneInfo.cons.MAYO_MSSM_ROSMAP_01.17.19.csv') #Final annotated geneInfo file is input to the rest of the analysis
save(list=ls(),file='geneInfo.cons.01.17.19.rda')

#=====================================================================================
#
#  Part 5: Module-trait relationships
#
#=====================================================================================

load("geneInfo.cons.112418.rda")

################################## Mayo ################################################
nSamples = nrow(datExpr.Ref.MAYO);
nGenes = ncol(datExpr.Ref.MAYO);

Diagnosis <- as.numeric(relevel(factor(as.character(targets.Ref.MAYO$Diagnosis)), 'CONTROL'))
Age <- as.numeric(targets.Ref.MAYO$AgeAtDeath)
RIN <- as.numeric(targets.Ref.MAYO$RIN)
Gender <- as.numeric(factor(targets.Ref.MAYO$Gender))
PCT_PF_READS_ALIGNED <- as.numeric(targets.Ref.MAYO$PCT_PF_READS_ALIGNED)

factors1_MAYO <- cbind(Diagnosis, Age, Gender, RIN, PCT_PF_READS_ALIGNED)
PCvalues <- MEs_MAYO[,-ncol(MEs_MAYO)] #exclude grey

moduleTraitCor_MAYO = cor(PCvalues, factors1_MAYO, use = "p");
moduleTraitPvalue_MAYO = corPvalueStudent(moduleTraitCor_MAYO, nSamples);
colnames(moduleTraitPvalue_MAYO) = paste("p.value.", colnames(moduleTraitCor_MAYO), sep="");

## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g.
txtMat <-  signif(moduleTraitPvalue_MAYO, 2)
txtMat[txtMat>=0.05] <- ""
txtMat[txtMat <0.05&txtMat >0.01] <- "*"
txtMat[txtMat <0.01&txtMat >0.005] <- "**"
#txtMat[txtMat <0.005&txtMat >0] <- "***"

txtMat1 <- signif( moduleTraitCor_MAYO,2)
#we only want to look at pearson correlations in certain range
txtMat1[txtMat1 > -0.3&txtMat1<0.2] <- ""
textMatrix1 = paste( txtMat1, '\n', '(',txtMat ,')', sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( moduleTraitPvalue_MAYO),nrow=nrow(moduleTraitPvalue_MAYO))


#Plot heatmap
pdf(paste('NetworkPlot_MAYO_consensus_01.16.19.pdf'),width=16,height=30)
par( mar = c(8, 12, 3, 3) );
labeledHeatmap(Matrix = moduleTraitCor_MAYO,
                xLabels = colnames(factors1_MAYO),
                yLabels = rownames(moduleTraitPvalue_MAYO),
                ySymbols = rownames(moduleTraitPvalue_MAYO),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix1,
                setStdMargins = FALSE,
                cex.text = 1.5,
                zlim = c(-1, 1),
                cex.lab.x = 1.2,
                main = paste("Module-trait relationships")
)

#Plot eigengene heatmap
par(cex = 1.0)
plotEigengeneNetworks(MEs_MAYO, "Eigengene Network", marHeatmap = c(3,4,2,2),
                      marDendro = c(0,4,1,2),cex.adjacency = 0.3,plotDendrograms = TRUE,
                      xLabelsAngle = 90,heatmapColors=blueWhiteRed(100)[51:100])

#Plot boxplots, scatterplots
toplot=t(MEs_MAYO)
cols=substring(colnames(MEs_MAYO),3,20)
par(mfrow=c(4,4))
par(mar=c(5,6,4,2))
for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(as.vector(as.factor(targets.Ref.MAYO$Diagnosis)),c('CONTROL','AD')),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
  verboseScatterplot(x=as.numeric(targets.Ref.MAYO$AgeAtDeath),y=toplot[i,],xlab="Age",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
  boxplot(toplot[i,]~factor(targets.Ref.MAYO$Gender),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
}

dev.off()

################################## MSSM ################################################
nSamples = nrow(datExpr.Ref.MSSM);
nGenes = ncol(datExpr.Ref.MSSM);

Diagnosis <- as.numeric(relevel(factor(as.character(targets.Ref.MSSM$Diagnosis)), 'CONTROL'))
Age <- as.numeric(targets.Ref.MSSM$AOD)
RIN <- as.numeric(targets.Ref.MSSM$RIN)
Gender <- as.numeric(factor(targets.Ref.MSSM$SEX))
PCT_PF_READS_ALIGNED <- as.numeric(targets.Ref.MSSM$PCT_PF_READS_ALIGNED)
CDR <- as.numeric(targets.Ref.MSSM$CDR)

factors1_MSSM=cbind(Diagnosis, Age, Gender, RIN, PCT_PF_READS_ALIGNED, CDR)
PCvalues<-MEs_MSSM[,-ncol(MEs_MSSM)] #exclude grey


moduleTraitCor_MSSM= cor(PCvalues, factors1_MSSM, use = "p");
moduleTraitPvalue_MSSM = corPvalueStudent(moduleTraitCor_MSSM, nSamples);
colnames(moduleTraitPvalue_MSSM) = paste("p.value.", colnames(moduleTraitCor_MSSM), sep="");

## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g.
txtMat <-  signif(moduleTraitPvalue_MSSM,2)
txtMat[txtMat>=0.05] <- ""
txtMat[txtMat <0.05&txtMat >0.01] <- "*"
txtMat[txtMat <0.01&txtMat >0.005] <- "**"
#txtMat[txtMat <0.005&txtMat >0] <- "***"
txtMat1 <- signif( moduleTraitCor_MSSM,2)

#we only want to look at pearson correlations in certain range
txtMat1[txtMat1> -0.3&txtMat1<0.2] <- ""
textMatrix1 = paste( txtMat1, '\n', '(',txtMat ,')', sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( moduleTraitPvalue_MSSM),nrow=nrow(moduleTraitPvalue_MSSM))

#Plot heatmap
pdf(paste('NetworkPlot_MSSM_consensus_01.16.19.pdf'),width=16,height=30)
par( mar = c(8, 12, 3, 3) );
labeledHeatmap( Matrix = moduleTraitCor_MSSM,
                xLabels = colnames(factors1_MSSM),
                yLabels = rownames(moduleTraitPvalue_MSSM),
                ySymbols = rownames(moduleTraitPvalue_MSSM),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix1,
                setStdMargins = FALSE,
                cex.text = 1.5,
                zlim = c(-1, 1),
                cex.lab.x = 1.2,
                main = paste("Module-trait relationships")
);

#Plot eigengene heatmap
par(cex = 1.0)
plotEigengeneNetworks(MEs_MSSM, "Eigengene Network", marHeatmap = c(3,4,2,2),
                      marDendro = c(0,4,1,2),cex.adjacency = 0.3,plotDendrograms = TRUE,
                      xLabelsAngle = 90,heatmapColors=blueWhiteRed(100)[51:100])

#Plot boxplots, scatterplots
toplot=t(MEs_MSSM)
cols=substring(colnames(MEs_MSSM),3,20)
par(mfrow=c(4,4))
par(mar=c(5,6,4,2))
for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(as.vector(as.factor(targets.Ref.MSSM$Diagnosis)),c('CONTROL','AD')),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
  verboseScatterplot(x=as.numeric(targets.Ref.MSSM$AOD),y=toplot[i,],xlab="Age",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
  boxplot(toplot[i,]~factor(targets.Ref.MSSM$SEX),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
}

dev.off()

############################### ROSMAP ########################################
nSamples = nrow(datExpr.Ref.ROSMAP);
nGenes = ncol(datExpr.Ref.ROSMAP);

Diagnosis <- as.numeric(relevel(factor(as.character(targets.Ref.ROSMAP$Diagnosis)), 'CONTROL'))
Age <- as.numeric(targets.Ref.ROSMAP$age_death)
Gender <- as.numeric(factor(targets.Ref.ROSMAP$msex))
RIN <- as.numeric(targets.Ref.ROSMAP$RINcontinuous)
cogdx <- as.numeric(targets.Ref.ROSMAP$cogdx)
PCT_PF_READS_ALIGNED <- as.numeric(targets.Ref.ROSMAP$PCT_PF_READS_ALIGNED)

factors1_ROSMAP <- cbind(Diagnosis, Age, Gender, RIN, PCT_PF_READS_ALIGNED, cogdx)
PCvalues <- MEs_ROSMAP[,-ncol(MEs_ROSMAP)] #exclude grey

moduleTraitCor_ROSMAP <- cor(PCvalues, factors1_ROSMAP, use = "p");
moduleTraitPvalue_ROSMAP <- corPvalueStudent(moduleTraitCor_ROSMAP, nSamples);
colnames(moduleTraitPvalue_ROSMAP) = paste("p.value.", colnames(moduleTraitCor_ROSMAP), sep="");

## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g.
txtMat <-  signif(moduleTraitPvalue_ROSMAP,2)
txtMat[txtMat>=0.05] <- ""
txtMat[txtMat <0.05&txtMat >0.01] <- "*"
txtMat[txtMat <0.01&txtMat >0.005] <- "**"
#txtMat[txtMat <0.005&txtMat >0] <- "***"
txtMat1 <- signif( moduleTraitCor_ROSMAP,2)

#we only want to look at pearson correlations in certain range
txtMat1[txtMat1> -0.3&txtMat1<0.2] <- ""

textMatrix1 = paste( txtMat1, '\n', '(',txtMat ,')', sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( moduleTraitPvalue_ROSMAP),nrow=nrow(moduleTraitPvalue_ROSMAP))

#plotting
pdf(paste('NetworkPlot_ROSMAP_consensus_01.16.19.pdf'),width=16,height=30)
par( mar = c(8, 12, 3, 3) );
labeledHeatmap( Matrix = moduleTraitCor_ROSMAP,
                xLabels = colnames(factors1_ROSMAP),
                yLabels = rownames(moduleTraitPvalue_ROSMAP),
                ySymbols = rownames(moduleTraitPvalue_ROSMAP),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix1,
                setStdMargins = FALSE,
                cex.text = 1.5,
                zlim = c(-1, 1),
                cex.lab.x = 1.2,
                main = paste("Module-trait relationships")
);

#Plot eigengene heatmap
par(cex = 1.0)
plotEigengeneNetworks(MEs_ROSMAP, "Eigengene Network", marHeatmap = c(3,4,2,2), marDendro = c(0,4,1,2),
                      cex.adjacency = 0.3,plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(100)[51:100])

#Plot boxplots, scatterplots
toplot=t(MEs_ROSMAP)
cols=substring(colnames(MEs_ROSMAP),3,20)
par(mfrow=c(4,4))
par(mar=c(5,6,4,2))
for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(as.vector(as.factor(targets.Ref.ROSMAP$Diagnosis)),c('CONTROL','AD')),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
  verboseScatterplot(x=as.numeric(targets.Ref.ROSMAP$age_death),y=toplot[i,],xlab="Age",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
  boxplot(toplot[i,]~factor(targets.Ref.ROSMAP$msex),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
}

dev.off()

moduleTraitPval=cbind(moduleTraitPvalue_MAYO, moduleTraitPvalue_MSSM, moduleTraitPvalue_ROSMAP)
moduleTraitCor=cbind(moduleTraitCor_MAYO, moduleTraitCor_MSSM, moduleTraitCor_ROSMAP)
write.csv(moduleTraitPval,'moduleTraitPvalue_MAYO_MSSM_ROSMAP_01.16.19.csv')
write.csv(moduleTraitCor,'moduleTraitCor_MAYO_MSSM_ROSMAP_01.16.19.csv')

#=====================================================================================
#
#  Part 6: GO analysis
#
#=====================================================================================
# Makes bar plots of top enriched GO terms for each module

load(file='geneInfo.cons.01.17.19.rda')

dir.create("./geneInfo")
dir.create("./geneInfo/background/")
dir.create("./geneInfo/input/")
dir.create("./geneInfo/output/")

geneInfo.cons$SystemCode =rep("En",length=nrow(geneInfo.cons))
background=geneInfo.cons[,"Ensembl.Gene.ID"]
background=as.data.frame(background)

## Output files for GO elite
background <- cbind(background,rep("En",length=length(background)))
colnames(background) <- c("Source Identifier","SystemCode")
write.table(background,"./geneInfo/background/denominator.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

uniquemodcolors=unique(moduleColor.cons)
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']

# i = Number of modules
for(i in 1:length(uniquemodcolors)){
  thismod= uniquemodcolors[i]
  ind=which(colnames(geneInfo.cons)==paste("consensus.kME",thismod,sep=""))
  thisInfo=geneInfo.cons[geneInfo.cons$Initially.Assigned.Module.Color==thismod, c(1, 26, ind)] ##18=Ensembl.ID, 21="SystemCode",ind=kME value what tf is going on on this line. the c(1,19, ind) is confusing.
  colnames(thisInfo) <- c("Source Identifier","SystemCode","kME")
  write.table(thisInfo,file=paste("./geneInfo/input/",thismod,"_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}

# Run GO elite as nohupped shell script:
codedir <- "/home/vivek/bin/GO-Elite_v.1.2.5-Py"
pathname <- "~/WGCNA/geneInfo"
nperm=10000
system(paste("nohup python ",codedir,"/GO_Elite.py --species Hs --mod Ensembl --permutations ",
             nperm,"  --method \"z-score\" --zscore 1.96 --pval 0.01 --num 5 --input ",pathname,
             "/input --denom ",pathname,"/background --output ",pathname,"/output &",sep=""))

# Plotting the GO Output
pathname <- "~/WGCNA/geneInfo/output/GO-Elite_results/CompleteResults"

uniquemodcolors=uniquemodcolors[-c(2, 14)] # For some reason sometimes modules are not run correctly, therefore won't be able to be plotted so they are excluded

#manually set uniquemodcolors:
uniquemodcolors = c("black", "brown", "darkgreen", "darkred", "darkturquoise", "greenyellow",
                    "lightcyan", "lightgreen", "lightyellow", "magenta", "midnightblue", "orange",
                    "pink", "purple", "royalblue", "tan", "grey60")

#for some reason the royalblue module gives errors so I removed that for now
uniquemodcolors = c("black", "brown", "darkgreen", "darkred", "darkturquoise", "greenyellow",
                    "lightcyan", "lightgreen", "lightyellow", "magenta", "midnightblue", "orange",
                    "pink", "purple", "tan", "grey60")


pdf("GOElite_plot_Modules_01.23.19.pdf",height=8,width=12)
for(i in 1:length(uniquemodcolors)){
  thismod = uniquemodcolors[i]
  tmp=read.csv(file=paste(pathname,"/ORA_pruned/",thismod,"_Module-GO_z-score_elite.txt",sep=""),sep="\t")
  tmp=subset(tmp,Ontology.Type!='cellular_component')
  tmp=tmp[,c(2,9)] ## Select GO-terms and Z-score
  tmp=tmp[order(tmp$Z.Score,decreasing=T),] #
  if (nrow(tmp)<10){
    tmp1=tmp ## Take top 10 Z-score
    tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
    par(mar=c(5,40,5,2))
    barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
    abline(v=2,col="red")
  } else {
    tmp1=tmp[c(1:10),] ## Take top 10 Z-score
    tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
    par(mar=c(5,40,5,2))
    barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
    abline(v=2,col="red")
    }

  cat('Done ...',thismod,'\n')
}

dev.off()

#=====================================================================================
#
#  Part 7: Cell type enrichment
#
#=====================================================================================
geneInfo=read.csv('geneInfo.cons.MAYO_MSSM_ROSMAP_01.17.19.csv')

datKME <- geneInfo[,c("Ensembl.Gene.ID","Initially.Assigned.Module.Color")] ## Get a list of genes to test for enrichment, e.g. genes with modules defined
testbackground <- as.character(geneInfo$Ensembl.Gene.ID) # background list
datKME=subset(datKME,Initially.Assigned.Module.Color!="grey")
namestestlist <- names(table(datKME[,2])) ## module
multiTest <- vector(mode = "list", length = length(namestestlist))
names(multiTest) <- namestestlist

for (i in 1:length(multiTest))
{
  multiTest[[i]] <- datKME[datKME[,2]==namestestlist[i],1]
}

## From Zhang et al., 2014 - ## CSV file with reference genes in 1st column, annotated list with Category in 2nd
datCells <- read.csv("/home/vivek//bin/ZhangEtAlCellTypeList_humanENSG.csv")

## Set up reference lists
namesreflist <- names(table(datCells[,2])) ## category or module color
multiRef <- vector(mode = "list", length = length(namesreflist))
names(multiRef) <- namesreflist
for (i in 1:length(multiRef))
{
  multiRef[[i]] <- datCells[datCells[,2]==namesreflist[i],3]
}

refbackground<- testbackground

source('/home/vivek//AD/Zhang/ORA.R')

ORA.OR = matrix(NA,nrow=length(multiTest),ncol=length(multiRef));
colnames(ORA.OR) = names(multiRef);
rownames(ORA.OR) = names(multiTest);
ORA.P = matrix(NA,nrow=length(multiTest),ncol=length(multiRef));
colnames(ORA.P) = names(multiRef);
rownames(ORA.P) = names(multiTest);

for (i in 1:length(multiRef)) {
  for (j in 1:length(multiTest)) {
    result = ORA(multiTest[[j]],multiRef[[i]],testbackground,refbackground);
    ORA.OR[j,i] = result[1];
    ORA.P[j,i] = result[2];
  }
}

ORA.OR<-apply(ORA.OR,2,as.numeric)
dim(ORA.OR)<-dim(ORA.P)

##FDR correct

FDRmat.Array <- matrix(p.adjust( ORA.P,method="fdr"),nrow=nrow( ORA.P),ncol=ncol( ORA.P))
rownames(  FDRmat.Array)=rownames(ORA.P)
colnames(  FDRmat.Array)=colnames(ORA.P)

ORA.P=matrix(as.numeric(ORA.P),nrow=nrow( ORA.P),ncol=ncol( ORA.P))
ORA.OR=matrix(as.numeric(ORA.OR),nrow=nrow( ORA.OR),ncol=ncol( ORA.OR))
rownames(ORA.P) <- rownames(ORA.OR) <- rownames(  FDRmat.Array)
colnames(ORA.P) <- colnames(ORA.OR) <- colnames(  FDRmat.Array)

dispMat <- ORA.OR ## You can change this to be just log2(Bmat) if you want the color to reflect the odds ratios
#Use the text function with the FDR filter in labeledHeatmap to add asterisks
txtMat <-  ORA.OR
txtMat[FDRmat.Array >0.05] <- ""
txtMat[FDRmat.Array <0.05&FDRmat.Array >0.01] <- "*"
txtMat[FDRmat.Array <0.01&FDRmat.Array >0.005] <- "**"
txtMat[FDRmat.Array <0.005] <- "***"

txtMat1 <- signif( ORA.OR,2)
txtMat1[txtMat1<2] <- ""

textMatrix1 = paste( txtMat1, '\n', txtMat , sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( ORA.P),nrow=nrow( ORA.P))

# Make heatmap of modules and cell types: neurons, microglia, myelinating oligodendrocytes, astrocytes, and endothelial cells
#Got a warning message upon running this below block of code
pdf("CellTypeEnrich_WGCNAMods_01.23.19.pdf", width=6,height=10)
labeledHeatmap(Matrix=dispMat,
               yLabels=rownames(dispMat),
               yColorLabels=TRUE,
               xLabels= colnames(dispMat),
               colors=blueWhiteRed(40),
               textMatrix = textMatrix1,
               cex.lab.x=1.0,
               zlim=c(-0.1,3),
               main="Cell-type enrichment Heatmap")
dev.off()


#=====================================================================================
#
#  Part 8: TOM network plot
#
#=====================================================================================
# Will make node and edge plots with hub genes in the center surrounded by all other genes in each module

load("geneInfo.cons.01.17.19.rda")
load("wgcna_consensus_01.10.19.rda")

#Get the top connected genes in the module
uniquemodcolors = unique(moduleColor.cons);
uniquemodcolors <- uniquemodcolors[!uniquemodcolors %in% "grey"]
TOM.matrix = as.matrix(consTomDS);

pdf("ModuleNetworks.pdf",height=9,width=10);
for (mod in uniquemodcolors)  {
  numgenesingraph = 100;
  numconnections2keep = 1500;
  cat('module:',mod,'\n');
  geneInfo.cons=geneInfo.cons[geneInfo.cons$GeneSymbol!="NA",]
  colind = which(colnames(geneInfo.cons)==paste("consensus.kME",mod, sep=""));
  rowind = which(geneInfo.cons[,3]==mod);
  cat(' ',length(rowind),'probes in module\n');
  submatrix = geneInfo.cons[rowind,];
  orderind = order(submatrix[,colind],decreasing=TRUE);
  if (length(rowind) < numgenesingraph) {
    numgenesingraph = length(rowind);
    numconnections2keep = numgenesingraph * (numgenesingraph - 1);
  }
  cat('Making network graphs, using top',numgenesingraph,'probes and',numconnections2keep,'connections of TOM\n');
  submatrix = submatrix[orderind[1:numgenesingraph],];
  #Identify the columns in the TOM that correspond to these hub probes
  matchind = match(submatrix$Ensembl.Gene.ID,colnames(multiExpr));
  reducedTOM = TOM.matrix[matchind,matchind];

  orderind = order(reducedTOM,decreasing=TRUE);
  connections2keep = orderind[1:numconnections2keep];
  reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
  reducedTOM[connections2keep] = 1;

  g0 <- graph.adjacency(as.matrix(reducedTOM[1:10,1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMata <- layout.circle(g0)    

  g0 <- graph.adjacency(as.matrix(reducedTOM[11:50,11:50]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMatb <- layout.circle(g0)

  g0 <- graph.adjacency(as.matrix(reducedTOM[51:ncol(reducedTOM),51:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMatc <- layout.circle(g0)
  g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.8, layoutMatc)

  plot(g1,edge.color="grey",vertex.color=mod,vertex.label=as.character(submatrix$GeneSymbol),vertex.label.cex=0.7,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat,vertex.size=submatrix[,colind]^2*8,main=paste(mod,"module"))

}
dev.off();

################################################################################
#                                 bookmark                                     #
################################################################################
