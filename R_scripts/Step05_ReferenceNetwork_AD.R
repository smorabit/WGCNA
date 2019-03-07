setwd("~/WGCNA_2")
getwd()
library(WGCNA)
library(flashClust)
library(gplots)
library(cluster)
library(igraph); #for part 8
library(RColorBrewer); #for part 8
options(stringsAsFactors=FALSE)
enableWGCNAThreads()

#===============================================================================
#
#  Part 1: Load data
#
#===============================================================================

load("~/data/MAYO_TCX_ADcontrolOnly_Expression_metaData.rda")

datExpr=as.data.frame(t(normExpr.TCX))
normExpr=as.data.frame(t(datExpr))
save(targets.TCX,datExpr,normExpr,file="DiscoverySet_AD.rda")

#################################################

#===============================================================================
#
#  Part 2: Choose soft thresholding power
#
#===============================================================================

pdf("1.1_power.pdf", height=10, width=18)
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType="signed",corFnc="bicor")
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
abline(h=0.80,col="blue")
abline(h=0.70,col="orange")
abline(h=0.60,col="green")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#=====================================================================================
#
#  Part 3: Network Construction
#
#=====================================================================================

softPower=12 ## for AD discovery set

adjacency = adjacency(datExpr, power = softPower, type = "signed",corFnc="bicor",corOptions="use='pairwise.complete.obs'"); #this line takes a while to run

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = flashClust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
save(TOM,dissTOM,geneTree,adjacency,softPower,file="TOM.rda")

mColorh <- mLabelh <- colorLabels <- NULL
  for (minModSize in c(40,100,160)) {
    for (dthresh in c(0.1,0.2,0.25)) {
      for (ds in c(2,4)) {
        print("Trying parameters:")
        print(c(minModSize,dthresh,ds))
        tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE,
          minClusterSize = minModSize, cutHeight = 0.9999,
          deepSplit = ds, distM = as.matrix(dissTOM))

        merged <- mergeCloseModules(exprData = datExpr,colors = tree$labels,
                                    cutHeight = dthresh)
        mColorh <- cbind(mColorh,labels2colors(merged$colors))
        mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
      }
    }
  }

### Relating dendrogram with traits
#datTraits<- targets.TCX[,c(5,7,6,4,41,42)]
datTraits<- targets.TCX[,c(4,6,7,8,9,11,16,17)]

traitmat=as.data.frame(cbind(as.factor(datTraits[,1]),
                              as.factor(datTraits[,2]),
                              as.numeric(datTraits[,3]),
                              as.numeric(datTraits[,4]),
                              as.numeric(datTraits[,5]),
                              as.numeric(datTraits[,6]),
                              as.numeric(as.factor(datTraits[,7])),
                              as.numeric(as.factor(datTraits[,8]))))


rownames(traitmat)=rownames(datTraits)
colnames(traitmat)=c("Sex","BrainRegion.Diagnosis","RIN","RIN2","Age", "PCT_PF_READS_ALIGNED","Region", "Diagnosis")

geneSigs=matrix(NA,nrow=8,ncol=ncol(datExpr)) # create a vector to hold the data

for(i in 1:ncol(geneSigs)) {
	exprvec=as.numeric(datExpr[,i]) # get the expression vector for ith gene
	sexr=bicor(exprvec, traitmat[,1],use="pairwise.complete.obs")
	braindiagr=bicor(traitmat[,2],exprvec,use="pairwise.complete.obs")
	#rinr=sqrt(max(summary(lm(exprvec~as.factor(traitmat[,3])))$adj.r.squared,0)) # calculate adjusted R^2s square-root for categorical variables
  rinr=bicor(traitmat[,3],exprvec,use="pairwise.complete.obs")
  rin2r=bicor(traitmat[,4],exprvec,use="pairwise.complete.obs")# calculate r correlation value for numeric variables
  ager=bicor(traitmat[,5],exprvec,use="pairwise.complete.obs")# calculate r correlation value for numeric variables
  pct_pf_reads_alignedr=bicor(traitmat[,6],exprvec,use="pairwise.complete.obs")
	regionr=bicor(traitmat[,7],exprvec,use="pairwise.complete.obs")
  diagnosisr=bicor(traitmat[,8],exprvec,use="pairwise.complete.obs")

	geneSigs[,i]=c(sexr, braindiagr, rinr, rin2r,ager,pct_pf_reads_alignedr, regionr, diagnosisr)

	cat('Done for gene...',i,'\n')
} #there were a ton of warnings...

for (i in 1:nrow(geneSigs)){
  geneSigs[i,] <- numbers2colors(as.numeric(geneSigs[i,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
}

rownames(geneSigs)=c("Sex","BrainRegion.Diagnosis","RIN","RIN2","Age", "PCT_PF_READS_ALIGNED","Region", "Diagnosis")

rm(TOM,adjacency,dissTOM)
save(list=ls(),file="WGCNA.rda")
load("TOM.rda")

############# Calculate modules for each set of parameters #####################
mColorh <- mLabelh <- colorLabels <- NULL
for (minModSize in c(40,100,160)) {
  for (dthresh in c(0.1,0.2,0.25)) {
    for (ds in c(2,4)) {
      print("Trying parameters:")
      print(c(minModSize,dthresh,ds))
      tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE,
                          minClusterSize = minModSize, cutHeight = 0.99999999,
                          deepSplit = ds, distM = as.matrix(dissTOM))

      merged <- mergeCloseModules(exprData = datExpr,colors = tree$labels,
                                  cutHeight = dthresh)
      mColorh <- cbind(mColorh,labels2colors(merged$colors))
      mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
    }
  }
}

# Plotting modules for each set of params and traits
mColorh1 <- cbind(mColorh, geneSigs[1,], geneSigs[2,], geneSigs[3,], geneSigs[4,],
                  geneSigs[5,], geneSigs[6,], geneSigs[7,], geneSigs[8,])


rownames_geneSigs = c(rownames(geneSigs))
mLabelh1=c(mLabelh,rownames_geneSigs)

pdf("MAYO_TOM_MultiDendro_DSAD_02.07.19.pdf",height=25,width=20)
plotDendroAndColors(geneTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed network (MAYO) with power = 8"));
dev.off()

### selecting parameters for the final dendrogram based on observations ###
#select dendrogram parameters
mms <- 100
ds <- 4
dthresh <- 0.2

tree <- cutreeHybrid(dendro = geneTree, pamStage=FALSE,
                    minClusterSize = mms, cutHeight = 0.99999999,
                    deepSplit = ds, distM = as.matrix(dissTOM))

merged <- mergeCloseModules(exprData = datExpr,colors = tree$labels,
                            cutHeight = dthresh)

# Eigengenes of the new merged modules:
MEs <- merged$newMEs

# Module color associated with each gene
moduleColor.cons <- labels2colors(merged$colors)

mColorh <- cbind(labels2colors(merged$colors))
mLabelh <- c("colors")

mColorh1 <- cbind(mColorh, geneSigs[1,], geneSigs[2,], geneSigs[3,], geneSigs[4,],
                  geneSigs[5,], geneSigs[6,], geneSigs[7,], geneSigs[8,])

rownames_geneSigs = c(rownames(geneSigs))
mLabelh1=c(mLabelh,rownames_geneSigs)

pdf("MAYO_TOM_softpower16_02.20.19.pdf",height=10,width=16)
plotDendroAndColors(geneTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network (MAYO) with power = 16, mms=",mms,"ds=",ds,"dthresh=",dthresh,"cquant=0.2"));
dev.off()


#### stuff for later ####

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

save(list=ls(),file="MAYO_WGCNA.rda")
