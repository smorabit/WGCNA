library(WGCNA)
setwd("~/modulePreservation")

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda') #this is where the file will be (no file now)
load('mp_rTG4510mouse_Mayo_AD.rda')

setLabels=c("Mayo.AD","rTG4510mouse")
shortLabels=setLabels
nSets=2

MEs = moduleEigengenes(datExpr.rTG4510, moduleColors.Mayo)$eigengenes

##############################################################################
# change below code from mt sinai to whatever other dataset you are working on
##############################################################################

nSamples = nrow(datExpr.rTG4510);
nGenes = ncol(datExpr.rTG4510);

##### next time pick labels to correlate!!

Diagnosis <- as.numeric(relevel(factor(as.character(targets.Ref.MSSM$Diagnosis)), 'CONTROL'))
Age <- as.numeric(targets.Ref.MSSM$AOD)
RIN <- as.numeric(targets.Ref.MSSM$RIN)
Gender <- as.numeric(factor(targets.Ref.MSSM$SEX))PCT_PF_READS_ALIGNED <- as.numeric(targets.Ref.MSSM$PCT_PF_READS_ALIGNED)
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
