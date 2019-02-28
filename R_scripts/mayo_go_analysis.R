library(WGCNA)
library(flashClust)
library(gplots)
library(cluster)
library(igraph); #for part 8
library(RColorBrewer); #for part 8
options(stringsAsFactors=FALSE)
enableWGCNAThreads()

#=====================================================================================
#  Part 0: Create geneInfo.cons dataframe from ME and datExpr objects
#=====================================================================================

# load relevant .rda files:
load("/pub/smorabit/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda")
load("/pub/smorabit/Mayo/rWGCNA/rWGCNA.rda")

# compute signed eigengene based connectivity (also known as module membership)
kme <- signedKME(datExpr.Ref, MEs.cons_Ref)

# put relevant info into a dataframe, then save it!
geneInfo.cons <- as.data.frame(cbind(colnames(datExpr.Ref), colnames(datExpr.Ref), moduleColors.cons, kme))
colnames(geneInfo.cons)[1]= "Ensembl.Gene.ID"
colnames(geneInfo.cons)[2]= "GeneSymbol"
colnames(geneInfo.cons)[3]= "Initially.Assigned.Module.Color"
write.csv(geneInfo.cons,'geneInfo.cons.MAYO.csv') #Final annotated geneInfo file is input to the rest of the analysis

#=====================================================================================
#  Part 1: GO analysis
#=====================================================================================

dir.create("./geneInfo")
dir.create("./geneInfo/background/")
dir.create("./geneInfo/input/")
dir.create("./geneInfo/output/")

geneInfo.cons$SystemCode <- rep("En",length=nrow(geneInfo.cons))
background <- geneInfo.cons[,"Ensembl.Gene.ID"]
background <- as.data.frame(background)

## Output files for GO elite
background <- cbind(background,rep("En",length=length(background)))
colnames(background) <- c("Source Identifier","SystemCode")
write.table(background,"./geneInfo/background/denominator.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

uniquemodcolors <- unique(moduleColors.Mayo)
uniquemodcolors <- uniquemodcolors[uniquemodcolors!='grey'] #do I need to exclude grey60?

# i = Number of modules
for(i in 1:length(uniquemodcolors)){
  thismod= uniquemodcolors[i]
  ind=which(colnames(geneInfo.cons)==paste("kME",thismod,sep=""))
  thisInfo=geneInfo.cons[geneInfo.cons$Initially.Assigned.Module.Color==thismod, c(1, dim(geneInfo.cons)[2], ind)]
  colnames(thisInfo) <- c("Source Identifier","SystemCode","kME")
  write.table(thisInfo,file=paste("./geneInfo/input/",thismod,"_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}

# Run GO elite as nohupped shell script:
codedir <- "~/bin/GO-Elite_v.1.2.5-Py"
pathname <- "/pub/smorabit/WGCNA/geneInfo"
nperm=10000
system(paste("nohup python ",codedir,"/GO_Elite.py --species Hs --mod Ensembl --permutations ",
             nperm,"  --method \"z-score\" --zscore 1.96 --pval 0.01 --num 5 --input ",pathname,
             "/input --denom ",pathname,"/background --output ",pathname,"/output &",sep=""))
#
# #########################################
# # stopped here last time
# #########################################
#
# # Plotting the GO Output
# pathname <- "~/WGCNA/geneInfo/output/GO-Elite_results/CompleteResults"
#
# uniquemodcolors=uniquemodcolors[-c(2, 14)] # For some reason sometimes modules are not run correctly, therefore won't be able to be plotted so they are excluded
#
# #manually set uniquemodcolors:
# uniquemodcolors = c("black", "brown", "darkgreen", "darkred", "darkturquoise", "greenyellow",
#                     "lightcyan", "lightgreen", "lightyellow", "magenta", "midnightblue", "orange",
#                     "pink", "purple", "royalblue", "tan", "grey60")
#
# #for some reason the royalblue module gives errors so I removed that for now
# uniquemodcolors = c("black", "brown", "darkgreen", "darkred", "darkturquoise", "greenyellow",
#                     "lightcyan", "lightgreen", "lightyellow", "magenta", "midnightblue", "orange",
#                     "pink", "purple", "tan", "grey60")
#
#
# pdf("GOElite_plot_Modules_01.23.19.pdf",height=8,width=12)
# for(i in 1:length(uniquemodcolors)){
#   thismod = uniquemodcolors[i]
#   tmp=read.csv(file=paste(pathname,"/ORA_pruned/",thismod,"_Module-GO_z-score_elite.txt",sep=""),sep="\t")
#   tmp=subset(tmp,Ontology.Type!='cellular_component')
#   tmp=tmp[,c(2,9)] ## Select GO-terms and Z-score
#   tmp=tmp[order(tmp$Z.Score,decreasing=T),] #
#   if (nrow(tmp)<10){
#     tmp1=tmp ## Take top 10 Z-score
#     tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
#     par(mar=c(5,40,5,2))
#     barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
#     abline(v=2,col="red")
#   } else {
#     tmp1=tmp[c(1:10),] ## Take top 10 Z-score
#     tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
#     par(mar=c(5,40,5,2))
#     barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
#     abline(v=2,col="red")
#     }
#
#   cat('Done ...',thismod,'\n')
# }
#
# dev.off()
