library(WGCNA)
library(flashClust)
library(gplots)
library(cluster)
library(igraph);
library(RColorBrewer);
library(variancePartition)
library(doParallel)
library(dplyr)
library(readxl)

options(stringsAsFactors=TRUE)
enableWGCNAThreads()
#########################################################################
# load / configure datasets
#########################################################################

# load ROSMAP metadata and DLPFC expression data:
rosmap_metadata <- read.csv("/home/smorabit/data/rosmap_RNAseq_Metadata_Complete_updated.csv")
load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/ROSMAP/ROSMAP_DLPFC_expression_metaData.rda')
actical_meta <- as.data.frame(read_excel("collab.xlsx"))


datExpr.Ref.ROSMAP <- t(normExpr.ROSMAP)
targets.Ref.ROSMAP <- targets.ROSMAP
ids <- intersect(targets.Ref.ROSMAP$SampleID, actical_meta$mrna_id)
datExpr.actical <- subset(datExpr.ROSMAP, rownames(datExpr.ROSMAP) %in% ids)
ROSMAP_meta <- subset(targets.Ref.ROSMAP, targets.Ref.ROSMAP$SampleID %in% ids)
actical_meta <- subset(actical_meta, actical_meta$mrna_id %in% ids )

nSets <- 1
setLabels <- c("actical")
multiExpr <- vector(mode="list", length=length(setLabels))
multiExpr[[1]] <- list(data=datExpr.actical)
exprSize <- checkSets(multiExpr)

#########################################################################
# Select soft power threshold
#########################################################################

# Choose a set of soft-thresholding powers
powers <- c(seq(1,10,by=1), seq(12,30, by=2));
powerTables <- vector(mode = "list", length = 1);
# Call the network topology analysis function for each set in turn
powerTables[[1]] <- list(data = pickSoftThreshold(datExpr.actical, powerVector=powers,verbose = 5, networkType="signed", corFnc="bicor")[[2]]);

# Plot the results
pdf("power_actical.pdf", height=10, width=18)
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

#########################################################################
# Network construction
#########################################################################

#Set soft thresholding power to number based on plots
softpower=20;

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

#save progress up to this point:
save(list=ls(),file="wgcna_actical.rda")

consMEs = net$multiMEs;
moduleColors = net$colors;
table(moduleColors)

load("ConsensusTOM-block.1.rda") # consensus TOM
consTree= hclust(1-consTomDS,method="average");

# correlate modules to traits

# cogdx, msex, age_death, pmi, last_tactivity_d_avg_new, tactivity_acth_avgnew, cogn_global_lv
cogdx <- as.numeric(actical_meta$cogdx)
age <- as.numeric(actical_meta$age_death)
gender <- as.numeric(factor(actical_meta$msex))
pmi <- as.numeric(actical_meta$pmi)
last_tactivity_d_avg_new <- as.numeric(actical_meta$last_tactivity_d_avgnew)
tactivity_acth_avgnew <- as.numeric(actical_meta$tactivity_acth_avgnew)
cogn_global_lv <- as.numeric(actical_meta$cogn_global_lv)

geneSigsActical <- matrix(NA,nrow=7,ncol=ncol(datExpr.actical))
for(i in 1:ncol(geneSigsActical)) {
  exprvec=as.numeric(datExpr.actical[,i])
  ager=bicor(age,exprvec,use="pairwise.complete.obs")
  sexr=bicor(exprvec, gender,use="pairwise.complete.obs")
  cogdxr=bicor(exprvec, cogdx ,use="pairwise.complete.obs")
  cognr=bicor(exprvec, cogn_global_lv, use="pairwise.complete.obs")
  pmir=bicor(exprvec, pmi, use="pairwise.complete.obs")
  lastr=bicor(exprvec, last_tactivity_d_avg_new, use="pairwise.complete.obs")
  tactivityr=bicor(exprvec, tactivity_acth_avgnew, use="pairwise.complete.obs")
  geneSigsActical[,i]=c(ager, sexr, cogdxr, cognr, pmir, lastr, tactivityr)
  cat('Done for gene...',i,'\n')
}

#set colors
for (i in 1:nrow(geneSigsActical)){
  geneSigsActical[i,] <- numbers2colors(as.numeric(geneSigsActical[i,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
}

rownames(geneSigsActical) <- c("Age","Sex","cogdx", "cogn_global_lv", "pmi", "last_tactivity_d_avgnew", "tactivity_acth_avgnew")

#########################################################################
# Calculate modules
#########################################################################
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
mColorh1=cbind(mColorh,geneSigsActical[1,],geneSigsActical[2,],geneSigsActical[3,],geneSigsActical[4,],geneSigsActical[5,])
rownames_geneSigs = c(rownames(geneSigsActical))
mLabelh1=c(mLabelh,rownames_geneSigs)

pdf("ConsensusTOM.pdf",height=25,width=20)
plotDendroAndColors(consTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed consensus network with power = 20"));
dev.off()

# Final dendrogram
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
MEs <- merged$newMEs[[1]]$data

# Module color associated with each gene
moduleColor.cons <- labels2colors(merged$colors)

mColorh <- cbind(labels2colors(merged$colors))
mLabelh <- c("Merged Colors")
mColorh1 <- cbind(mColorh,geneSigsActical[1,],geneSigsActical[2,],geneSigsActical[3,],geneSigsActical[4,],geneSigsActical[5,], geneSigsActical[6,], geneSigsActical[7,])
rownames_geneSigs <- c(rownames(geneSigsActical))
mLabelh1 <- c(mLabelh,rownames_geneSigs)

pdf("ConsensusTOM_Final.pdf",height=10,width=16)
plotDendroAndColors(consTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power = 20, mms=",mms,"ds=",ds,"dthresh=",dthresh,"cquant=0.2"));
dev.off()

# Convert numerical lables to colors for labeling of modules
MEColors <- labels2colors(as.numeric(substring(names(MEs), 3)));
MEColorNames <- paste("ME", MEColors, sep="");
colnames(MEs)=c(MEColorNames)

#########################################################################
#  Get annotated gene list:
#########################################################################

ensembl <- read.csv("/home/vivek/FTD_Seeley/Analysis_Nov2017/ENSG85_Human.csv.gz") # Convert Ensembl gene ID to gene names
MEs$Ensembl.Gene.ID <- paste(rownames(MEs))

##################################
nSamples = nrow(datExpr.actical);
nGenes = ncol(datExpr.actical);

cogdx <- as.numeric(actical_meta$cogdx)
age <- as.numeric(actical_meta$age_death)
gender <- as.numeric(factor(actical_meta$msex))
pmi <- as.numeric(actical_meta$pmi)
last_tactivity_d_avg_new <- as.numeric(actical_meta$last_tactivity_d_avgnew)
tactivity_acth_avgnew <- as.numeric(actical_meta$tactivity_acth_avgnew)
cogn_global_lv <- as.numeric(actical_meta$cogn_global_lv)

factors1_actical <- cbind(cogdx, age, gender, pmi, last_tactivity_d_avg_new, tactivity_acth_avgnew, cogn_global_lv)
PCvalues <- MEs[,-ncol(MEs)] #exclude grey

moduleTraitCor = cor(PCvalues, factors1_actical, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
colnames(moduleTraitPvalue) = paste("p.value.", colnames(moduleTraitCor), sep="");

## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g.
txtMat <-  signif(moduleTraitPvalue, 2)
txtMat[txtMat>=0.05] <- ""
txtMat[txtMat <0.05&txtMat >0.01] <- "*"
txtMat[txtMat <0.01&txtMat >0.005] <- "**"
#txtMat[txtMat <0.005&txtMat >0] <- "***"

txtMat1 <- signif( moduleTraitCor,2)
#we only want to look at pearson correlations in certain range
txtMat1[txtMat1 > -0.3&txtMat1<0.2] <- ""
textMatrix1 = paste( txtMat1, '\n', '(',txtMat ,')', sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( moduleTraitPvalue),nrow=nrow(moduleTraitPvalue))


#Plot heatmap
pdf(paste('NetworkPlot.pdf'),width=16,height=30)
par( mar = c(8, 12, 3, 3) );
labeledHeatmap(Matrix = moduleTraitCor,
                xLabels = colnames(factors1_actical),
                yLabels = rownames(moduleTraitPvalue),
                ySymbols = rownames(moduleTraitPvalue),
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
plotEigengeneNetworks(MEs, "Eigengene Network", marHeatmap = c(3,4,2,2),
                      marDendro = c(0,4,1,2),cex.adjacency = 0.3,plotDendrograms = TRUE,
                      xLabelsAngle = 90,heatmapColors=blueWhiteRed(100)[51:100])

#Plot boxplots, scatterplots
# toplot=t(MEs)
# cols=substring(colnames(MEs),3,20)
# par(mfrow=c(4,4))
# par(mar=c(5,6,4,2))
# for (i in 1:nrow(toplot)) {
#   boxplot(toplot[i,]~factor(as.vector(as.factor(targets.Ref.MAYO$Diagnosis)),c('CONTROL','AD')),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
#   verboseScatterplot(x=as.numeric(targets.Ref.MAYO$AgeAtDeath),y=toplot[i,],xlab="Age",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
#   boxplot(toplot[i,]~factor(targets.Ref.MAYO$Gender),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
# }

dev.off()

#########################################################################
#  GO Analysis
#########################################################################

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
