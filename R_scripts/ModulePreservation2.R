
library(WGCNA)
load('mp_ROSMAP_Mayo_AD.rda') ### load corresponding file for each module preservation data
mp=modulePreservation(multiExpr,multiColors,referenceNetworks=1,nPermutations=1000,networkType = "signed", corFnc="bicor",randomSeed=1,quickCor=0,verbose=3)
save(list=ls(),file=filename)

ref = 1 ## is the one whose moduleColors is provided
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1],mp$preservation$log.pBonf[[ref]][[test]][,-1])
zstats=signif(statsZ[, c("Zsummary.pres", "log.p.Bonfsummary.pres")], 2)

fle=paste("Zstats_",setLabels[1],"_vs_",setLabels[2],".csv",sep="")
write.csv(zstats,file=fle)

##Plotting
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
text1=text

# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
#sizeGrWindow(10, 5);

# Start the plot
#pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
pdf(paste("ModulePreserv1_",setLabels[1],"_vs_",setLabels[2],".pdf",sep=""),height=5,width=10)
text1=text


par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))

for (p in 1:2)
{
min = min(plotData[, p], na.rm = TRUE);
max = max(plotData[, p], na.rm = TRUE);
# Adjust ploting ranges appropriately
if (p==2)
{
if (min > -max/10) min = -max/10
ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
} else
ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
# set.seed(12312)
plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,main = mains[p], cex = 2.4, ylab = mains[p], xlab = "Module size", log = "x", ylim = ylim, xlim = c(40, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
 # text(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offset = 0.7,pos=3)
 text(moduleSizes[plotMods], plotData[plotMods, p], text1, cex = 1, offset = 0.7,pos=3)

# plot(jitter(moduleSizes[plotMods],12), plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,main = mains[p], cex = 2.4, ylab = mains[p], xlab = "Module size", log = "x", ylim = ylim, xlim = c(90,120), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
#labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
# For Zsummary, add threshold lines
if (p==2)
{
abline(h=0)
abline(h=2, col = "blue", lty = 2)
abline(h=10, col = "darkgreen", lty = 2)
}
}
dev.off()
