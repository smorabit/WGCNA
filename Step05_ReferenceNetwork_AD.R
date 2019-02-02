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

load("~/data/MAYO_TCX_ADcontrolOnly_Expression_metaData.rda")

##########################

datExpr=as.data.frame(t(normExpr.TCX))
normExpr=as.data.frame(t(datExpr))
save(targets.TCX,datExpr,normExpr,file="DiscoverySet_AD.rda")

#################################################

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

##################################################

softPower=8 ## for AD discovery set

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
##stopped working here

 ### Relating dendrogram with traits
datTraits<- targets[,c(5,7,6,4,41,42)]

traitmat=as.data.frame(cbind(as.numeric(factor(datTraits[,1],c("Control","AD")))-1,
                              as.numeric(datTraits[,2]),
                              as.factor(datTraits[,3]),
                              as.numeric(datTraits[,4]),
                              as.numeric(datTraits[,5]),
                              as.numeric(datTraits[,6])))# convert categorical variables in factor and numeric as numeric


rownames(traitmat)=rownames(datTraits)
colnames(traitmat)=c("Diagnosis","Age","Sex","RIN","SeqPC1","SeqPC2")


### Stratifying by Age


geneSigs=matrix(NA,nrow=6,ncol=ncol(datExpr)) # create a vector to hold the data

for(i in 1:ncol(geneSigs)) {

	exprvec=as.numeric(datExpr[,i]) # get the expression vector for ith gene
	conditionr=bicor(exprvec, traitmat[,1],use="pairwise.complete.obs")
	ager=bicor(traitmat[,2],exprvec,use="pairwise.complete.obs")
	sexr=sqrt(max(summary(lm(exprvec~as.factor(traitmat[,3])))$adj.r.squared,0)) # calculate adjusted R^2s square-root for categorical variables
	rinr=bicor(traitmat[,4],exprvec,use="pairwise.complete.obs")# calculate r correlation value for numeric variables
	seq1r=bicor(traitmat[,5],exprvec,use="pairwise.complete.obs")
	seq2r=bicor(traitmat[,6],exprvec,use="pairwise.complete.obs")


	geneSigs[,i]=c(conditionr, ager, sexr, rinr,seq1r,seq2r)

	cat('Done for gene...',i,'\n')
}


geneSigs[1,] =numbers2colors(as.numeric(geneSigs[1,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[2,] =numbers2colors(as.numeric(geneSigs[2,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[3,] =numbers2colors(as.numeric(geneSigs[3,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1)) # For categorical variables like strain or wt_tg we do not want values, thus lim=c(0,1), and signed and centered=F

geneSigs[4,] =numbers2colors(as.numeric(geneSigs[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))

geneSigs[5,] =numbers2colors(as.numeric(geneSigs[5,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigs[6,] =numbers2colors(as.numeric(geneSigs[6,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))

rownames(geneSigs)=c("Diagnosis","Age","Sex","RIN","Seq.PC1","Seq.PC2")



 rm(TOM,adjacency,dissTOM)
 save(list=ls(),file="WGCNA.rda")
