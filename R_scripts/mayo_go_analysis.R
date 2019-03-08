library(WGCNA)
library(flashClust)
library(gplots)
library(cluster)
library(igraph); #for part 8
library(RColorBrewer); #for part 8
library(RegFacEnc)

options(stringsAsFactors=FALSE)
enableWGCNAThreads()


load("/home/vivek/AMP_AD/Mayo/Analysis/Step04_WGCNA/Discovery_Set/AD/rWGCNA_Mayo_ForPreservation.rda")

#=====================================================================================
#  Part 1: Module-trait relationship
#=====================================================================================
nSamples = nrow(datExpr.Ref)
nGenes = ncol(datExpr.Ref)

# get lists of relevant traits from metadata object, then collapse into one table
Diagnosis <- as.numeric(relevel(as.factor(targets.Ref$Diagnosis),'Control'))
Age <- as.numeric(targets.Ref$AgeAtDeath)
RIN <- as.numeric(targets.Ref$RIN)
Gender <- as.numeric(as.factor(targets.Ref$Gender))
reads_aligned <- as.numeric(targets.Ref$PCT_PF_READS_ALIGNED)
intergenic_bases <- as.numeric(targets.Ref$PCT_INTERGENIC_BASES)
intronic_bases <- as.numeric(targets.Ref$PCT_INTRONIC_BASES)
coding_bases <- as.numeric(targets.Ref$PCT_CODING_BASES)
ribosomal_bases <- as.numeric(targets.Ref$PCT_RIBOSOMAL_BASES)
pc1 <- as.numeric(targets.Ref$Seq.PC1)
pc2 <- as.numeric(targets.Ref$Seq.PC2)
traits_table <- cbind(Diagnosis, Age, RIN, Gender, reads_aligned, intergenic_bases,
                      intronic_bases, coding_bases, ribosomal_bases, pc1, pc2)

# correlate PCs to traits
PCvalues <- MEs.Mayo[,-ncol(MEs.Mayo)]
moduleTraitCor_MAYO <- cor(PCvalues, traits_table, use='p')
moduleTraitPvalue_MAYO <- corPvalueStudent(moduleTraitCor_MAYO, nSamples)
colnames(moduleTraitCor_MAYO) <- paste("p.value", colnames(moduleTraitCor_MAYO), sep="")

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
pdf(paste('module_trait_correlation_MAYO_03.07.19.pdf'),width=16,height=30)
par( mar = c(8, 12, 3, 3) );
labeledHeatmap(Matrix = moduleTraitCor_MAYO,
                xLabels = colnames(traits_table),
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
plotEigengeneNetworks(MEs.Mayo, "Eigengene Network", marHeatmap = c(3,4,2,2),
                      marDendro = c(0,4,1,2),cex.adjacency = 0.3,plotDendrograms = TRUE,
                      xLabelsAngle = 90,heatmapColors=blueWhiteRed(100)[51:100])

#Plot boxplots, scatterplots
toplot=t(MEs.Mayo)
cols=substring(colnames(MEs.Mayo),3,20)
par(mfrow=c(4,4))
par(mar=c(5,6,4,2))
for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(as.vector(as.factor(targets.Ref$Diagnosis)),c('CONTROL','AD')),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
  verboseScatterplot(x=as.numeric(targets.Ref$AgeAtDeath),y=toplot[i,],xlab="Age",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
  boxplot(toplot[i,]~factor(targets.Ref$Gender),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
}

dev.off()

#=====================================================================================
#  Part 2: GO analysis
#=====================================================================================

geneInfo.mayo <- read.csv("geneInfo.MAYO.csv")

dir.create("./geneInfo2")
dir.create("./geneInfo2/background/")
dir.create("./geneInfo2/input/")
dir.create("./geneInfo2/output/")

geneInfo.mayo$SystemCode <- rep("En",length=nrow(geneInfo.mayo))
background <- geneInfo.mayo[,"Ensembl.Gene.ID"]
background <- as.data.frame(background)

## Output files for GO elite
background <- cbind(background,rep("En",length=length(background)))
colnames(background) <- c("Source Identifier","SystemCode")
write.table(background,"./geneInfo2/background/denominator.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

uniquemodcolors <- unique(moduleColors.Mayo)
uniquemodcolors <- uniquemodcolors[uniquemodcolors!='grey'] #do I need to exclude grey60?

# i = Number of modules
for(i in 1:length(uniquemodcolors)){
  thismod= uniquemodcolors[i]
  ind=which(colnames(geneInfo.mayo)==paste("kME",thismod,sep=""))
  thisInfo=geneInfo.mayo[geneInfo.mayo$Initially.Assigned.Module.Color==thismod, c(1, dim(geneInfo.mayo)[2], ind)]
  colnames(thisInfo) <- c("Source Identifier","SystemCode","kME")
  write.table(thisInfo,file=paste("./geneInfo2/input/",thismod,"_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}

# Run GO elite as nohupped shell script:
codedir <- "/home/vivek/bin/GO-Elite_v.1.2.5-Py"
pathname <- "~/mayo_WGCNA/geneInfo2"
nperm=10000
system(paste("nohup python ",codedir,"/GO_Elite.py --species Hs --mod Ensembl --permutations ",
             nperm,"  --method \"z-score\" --zscore 1.96 --pval 0.01 --num 5 --input ",pathname,
             "/input --denom ",pathname,"/background --output ",pathname,"/output &",sep=""))


# Plotting the GO Output
pathname <- "~/mayo_WGCNA/geneInfo/output/GO-Elite_results/CompleteResults"

#had to get rid of cyan, yellow, lightgreen, blue, salmon
uniquemodcolors <- uniquemodcolors[uniquemodcolors!='salmon']

pdf("GOElite_plot_Modules.pdf",height=8,width=12)
for(i in 1:length(uniquemodcolors)){
  thismod = uniquemodcolors[i]
  cat('Starting ...',thismod,'\n')
  tmp=read.csv(file=paste(pathname,"/ORA_pruned/",thismod,"_Module-GO_z-score_elite.txt",sep=""),sep="\t")
  tmp=subset(tmp,Ontology.Type!='cellular_component')
  tmp=tmp[,c(2,9)] ## Select GO-terms and Z-score
  tmp=tmp[order(tmp$Z.Score,decreasing=T),] #
  if (nrow(tmp)<10){
    tmp1=tmp ## Take top 10 Z-score
    tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
    par(mar=c(5,40,5,2))
    #barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
    barplot(tmp1$Z.Score,horiz=T,col=thismod,names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
    abline(v=2,col="black")
  } else {
    tmp1=tmp[c(1:10),] ## Take top 10 Z-score
    tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
    par(mar=c(5,40,5,2))
    barplot(tmp1$Z.Score,horiz=T,col=thismod,names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
    abline(v=2,col="black")
    }

  cat('Done ...',thismod,'\n')
}

dev.off()

#=====================================================================================
#  Part 3: Cell type enrichment
#=====================================================================================
geneInfo <- read.csv('geneInfo.MAYO.csv')

# Get a list of genes to test for enrichment, e.g. genes with modules defined
datKME <- geneInfo[,c("Ensembl.Gene.ID","Initially.Assigned.Module.Color")]
testbackground <- as.character(geneInfo$Ensembl.Gene.ID) # background list
datKME <- subset(datKME,Initially.Assigned.Module.Color!="grey")
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
pdf("CellTypeEnrich_mayo_WGCNAMods.pdf", width=6,height=10)
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
#  Part 4: TOM network plot
#  Plots all modules as a network with hub genes in the center surrounded
#  by all other genes in the module.
#=====================================================================================

load("consensusTOM_final.rda")
geneInfo.mayo <- read.csv("geneInfo.MAYO.csv")

#Get the top connected genes in the module
uniquemodcolors = unique(moduleColors.Mayo);
uniquemodcolors <- uniquemodcolors[!uniquemodcolors %in% "grey"]
TOM.matrix = as.matrix(consensusTOM_final);

pdf("mayo_ModuleNetworks.pdf",height=9,width=10);
for (mod in uniquemodcolors)  {
  numgenesingraph = 100;
  numconnections2keep = 1500;
  cat('module:',mod,'\n');
  geneInfo=geneInfo.mayo[geneInfo.mayo$GeneSymbol!="NA",]
  colind = which(colnames(geneInfo.mayo)==paste("kME",mod, sep=""));
  rowind = which(geneInfo.mayo[,4]==mod);
  cat(' ',length(rowind),'probes in module\n');
  submatrix = geneInfo.mayo[rowind,];
  orderind = order(submatrix[,colind],decreasing=TRUE);
  if (length(rowind) < numgenesingraph) {
    numgenesingraph = length(rowind);
    numconnections2keep = numgenesingraph * (numgenesingraph - 1);
  }
  cat('Making network graphs, using top',numgenesingraph,'probes and',numconnections2keep,'connections of TOM\n');
  submatrix = submatrix[orderind[1:numgenesingraph],];

  #Identify the columns in the TOM that correspond to these hub probes
  matchind = match(submatrix$Ensembl.Gene.ID,colnames(datExpr.Ref));
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
