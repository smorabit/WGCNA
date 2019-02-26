library(WGCNA)
setwd("~/modulePreservation")

##############################################################################
# rTG4510 mouse [x]
##############################################################################

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda') #this is where the file will be (no file now)
load('mp_rTG4510mouse_Mayo_AD.rda')

setLabels=c("Mayo.AD","rTG4510mouse")
shortLabels=setLabels
nSets=2

MEs = moduleEigengenes(datExpr.rTG4510, moduleColors.Mayo)$eigengenes
MEs=orderMEs(MEs)
uniquemodcolors=unique(moduleColors.Mayo)
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']

Tg=as.character(targets_rTG4510$Genotype)
Tg[Tg=="-"] <- "Wt"
Tg[Tg=="+"] <- "Tg"
Age <- as.numeric(targets_rTG4510$Age_months)

ID=factor(paste(Tg,Age,sep="."))
ID=factor(ID,c('Wt.2.5','Tg.2.5','Wt.4.5','Tg.4.5','Wt.6','Tg.6'))

pdf("ME_Condition_rTG4510mouse_MayoRef.pdf",height=4,width=4.5)
for (j in 1:length(uniquemodcolors)){
	   thismod= uniquemodcolors[j]
	thisME <- MEs[,paste("ME",thismod,sep="")]
	# lm1 <- summary(lm(thisME~ID))
	boxplot(thisME~ID,col=c("blue","red","blue","red","blue","red"),main=paste(thismod, "Module ME ",sep=""), ylab="Module Eigengene Value")

}
dev.off()


nSamples = nrow(datExpr.rTG4510);
nGenes = ncol(datExpr.rTG4510);

Genotype <- as.numeric(as.factor(targets_rTG4510$Genotype))
Age <- as.numeric(targets_rTG4510$Age_months)
RIN <- as.numeric(targets_rTG4510$RIN)
Sex <- as.numeric(factor(targets_rTG4510$Sex))

factors1_rTG4510=cbind(Genotype, Age, RIN, Sex)
PCvalues <-MEs [,-ncol(MEs)] #exclude grey

moduleTraitCor_rTG4510= cor(PCvalues, factors1_rTG4510, use = "p");
moduleTraitPvalue_rTG4510 = corPvalueStudent(moduleTraitCor_rTG4510, nSamples);
colnames(moduleTraitPvalue_rTG4510) = paste("p.value.", colnames(moduleTraitCor_rTG4510), sep="");

## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g.
txtMat <-  signif(moduleTraitPvalue_rTG4510,2)
txtMat[txtMat>=0.05] <- ""
txtMat[txtMat <0.05&txtMat >0.01] <- "*"
txtMat[txtMat <0.01&txtMat >0.005] <- "**"
#txtMat[txtMat <0.005&txtMat >0] <- "***"
txtMat1 <- signif( moduleTraitCor_rTG4510,2)

#we only want to look at pearson correlations in certain range
txtMat1[txtMat1> -0.3&txtMat1<0.2] <- ""
textMatrix1 = paste( txtMat1, '\n', '(',txtMat ,')', sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( moduleTraitPvalue_rTG4510),nrow=nrow(moduleTraitPvalue_rTG4510))

#Plot heatmap
pdf(paste('NetworkPlot_rTG4510_consensus_02.20.19.pdf'),width=16,height=30)
par( mar = c(8, 12, 3, 3) );
labeledHeatmap( Matrix = moduleTraitCor_rTG4510,
                xLabels = colnames(factors1_rTG4510),
                yLabels = rownames(moduleTraitPvalue_rTG4510),
                ySymbols = rownames(moduleTraitPvalue_rTG4510),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix1,
                setStdMargins = FALSE,
                cex.text = 1.5,
                zlim = c(-1, 1),
                cex.lab.x = 1.2,
                main = paste("Module-trait relationships rTG4510")
);

#Plot eigengene heatmap
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene Network", marHeatmap = c(3,4,2,2),
                      marDendro = c(0,4,1,2),cex.adjacency = 0.3,plotDendrograms = TRUE,
                      xLabelsAngle = 90,heatmapColors=blueWhiteRed(100)[51:100])

#Plot boxplots, scatterplots
toplot=t(MEs)
cols=substring(colnames(MEs),3,20)
par(mfrow=c(4,4))
par(mar=c(5,6,4,2))
for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(as.vector(as.factor(targets_rTG4510$Genotype))),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
  verboseScatterplot(x=as.numeric(targets_rTG4510$Age_months),y=toplot[i,],xlab="Age (months)",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
  boxplot(toplot[i,]~factor(targets_rTG4510$Sex),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
}

dev.off()

##############################################################################
########## APP_Oligomers -- Dentate Gyrus [x]
##############################################################################
rm(list=ls())
load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda') #this is where the file will be (no file now)
load('mp_DG_Oligomer_mouse_Mayo_AD.rda')

setLabels=c("Mayo.AD","Oligomer_DG_mouse")
shortLabels=setLabels
nSets=2

####Trajectory

MEs=moduleEigengenes(datExpr.DG, colors = moduleColors.Mayo, nPC=1)$eigengenes
MEs=orderMEs(MEs)
uniquemodcolors=unique(moduleColors.Mayo)
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']

group=factor(targets$strain)
group=factor(group,c('WT','oligomerogenic','fibrillogenic'))

pdf("ME_Condition_OligomerMouse_DG_MayoRef.pdf",height=4,width=4.5)
for (j in 1:length(uniquemodcolors)){
	   thismod= uniquemodcolors[j]
	thisME <- MEs[,paste("ME",thismod,sep="")]
	boxplot(thisME~group,col=c("blue","red","green"),main=paste(thismod, "Module ME ",sep=""), ylab="Module Eigengene Value")

}
dev.off()

nSamples = nrow(datExpr.DG);
nGenes = ncol(datExpr.DG);

Genotype <- as.numeric(as.factor(targets_rTG4510$Genotype))
Age <- as.numeric(targets_rTG4510$Age_months)
RIN <- as.numeric(targets_rTG4510$RIN)
Sex <- as.numeric(factor(targets_rTG4510$Sex))

factors1_rTG4510=cbind(Genotype, Age, RIN, Sex)
PCvalues <-MEs [,-ncol(MEs)] #exclude grey

moduleTraitCor_rTG4510= cor(PCvalues, factors1_rTG4510, use = "p");
moduleTraitPvalue_rTG4510 = corPvalueStudent(moduleTraitCor_rTG4510, nSamples);
colnames(moduleTraitPvalue_rTG4510) = paste("p.value.", colnames(moduleTraitCor_rTG4510), sep="");

## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g.
txtMat <-  signif(moduleTraitPvalue_rTG4510,2)
txtMat[txtMat>=0.05] <- ""
txtMat[txtMat <0.05&txtMat >0.01] <- "*"
txtMat[txtMat <0.01&txtMat >0.005] <- "**"
#txtMat[txtMat <0.005&txtMat >0] <- "***"
txtMat1 <- signif( moduleTraitCor_rTG4510,2)

#we only want to look at pearson correlations in certain range
txtMat1[txtMat1> -0.3&txtMat1<0.2] <- ""
textMatrix1 = paste( txtMat1, '\n', '(',txtMat ,')', sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( moduleTraitPvalue_rTG4510),nrow=nrow(moduleTraitPvalue_rTG4510))

#Plot heatmap
pdf(paste('NetworkPlot_rTG4510_consensus_02.20.19.pdf'),width=16,height=30)
par( mar = c(8, 12, 3, 3) );
labeledHeatmap( Matrix = moduleTraitCor_rTG4510,
                xLabels = colnames(factors1_rTG4510),
                yLabels = rownames(moduleTraitPvalue_rTG4510),
                ySymbols = rownames(moduleTraitPvalue_rTG4510),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix1,
                setStdMargins = FALSE,
                cex.text = 1.5,
                zlim = c(-1, 1),
                cex.lab.x = 1.2,
                main = paste("Module-trait relationships rTG4510")
);

#Plot eigengene heatmap
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene Network", marHeatmap = c(3,4,2,2),
                      marDendro = c(0,4,1,2),cex.adjacency = 0.3,plotDendrograms = TRUE,
                      xLabelsAngle = 90,heatmapColors=blueWhiteRed(100)[51:100])

#Plot boxplots, scatterplots
toplot=t(MEs)
cols=substring(colnames(MEs),3,20)
par(mfrow=c(4,4))
par(mar=c(5,6,4,2))
for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(as.vector(as.factor(targets_rTG4510$Genotype))),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
  verboseScatterplot(x=as.numeric(targets_rTG4510$Age_months),y=toplot[i,],xlab="Age (months)",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
  boxplot(toplot[i,]~factor(targets_rTG4510$Sex),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
}

dev.off()

###############################################################################
########## APP_Oligomers -- Entorhinal Cortex [x]
###############################################################################
rm(list=ls())
load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda') #this is where the file will be (no file now)
load('mp_EC_Oligomer_mouse_Mayo_AD.rda')

setLabels=c("Mayo.AD","Oligomer_DG_mouse")
shortLabels=setLabels
nSets=2

MEs=moduleEigengenes(datExpr.EC, colors = moduleColors.Mayo, nPC=1)$eigengenes
MEs=orderMEs(MEs)
uniquemodcolors=unique(moduleColors.Mayo)
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']

group=factor(targets$strain)
group=factor(group,c('WT','oligomerogenic','fibrillogenic'))

pdf("ME_Condition_OligomerMouse_eCTX_MayoRef.pdf",height=4,width=4.5)
for (j in 1:length(uniquemodcolors)){
	   thismod= uniquemodcolors[j]
	thisME <- MEs[,paste("ME",thismod,sep="")]
	boxplot(thisME~group,col=c("blue","red","green"),main=paste(thismod, "Module ME ",sep=""), ylab="Module Eigengene Value")

}
dev.off()
