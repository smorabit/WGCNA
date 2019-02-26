



### APP Mouse
####Trajectory
load('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/CRND8/Normalized_APP.rda')
datExpr=as.data.frame(t(normExpr_APP))

ensembl=read.csv('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/CRND8/GRCm38_hg38_homolog.csv')
ensembl=subset(ensembl,Homology.Type=='ortholog_one2one')

gnS=intersect(ensembl$Ensembl.Gene.ID,colnames(datExpr))
ensembl=ensembl[match(gnS,ensembl$Ensembl.Gene.ID),]
datExpr=datExpr[,match(gnS,colnames(datExpr))]
colnames(datExpr) <- ensembl$Human.Ensembl.Gene.ID
datExpr.APP=datExpr


library(WGCNA)

goodSamplesGenes(datExpr.APP) -> tmp ## one gene has no variance; removing the gene
rm.Gene=colnames(datExpr.APP)[which(tmp$goodGenes==FALSE)] ##remove this gene
a=setdiff(colnames(datExpr.APP),rm.Gene)
datExpr.APP=datExpr.APP[,match(a,colnames(datExpr.APP))]



load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')
gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.APP))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.APP=datExpr.APP[,match(gnS,colnames(datExpr.APP))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


MEs=moduleEigengenes(datExpr.APP, colors = moduleColors.Mayo, nPC=1)$eigengenes
MEs=orderMEs(MEs)
uniquemodcolors=unique(moduleColors.Mayo)
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']

Tg=as.character(targets_APP$Genotype)
Tg[Tg=="-"] <- "Wt"
 Tg[Tg=="+"] <- "Tg"
Age <- as.numeric(targets_APP$Age_months)

ID=factor(paste(Tg,Age,sep="."))
ID=factor(ID,c('Wt.9','Tg.9','Wt.12','Tg.12','Wt.20','Tg.20'))

pdf("ME_Condition_APPmouse_MayoRef.pdf",height=4,width=4.5)
for (j in 1:length(uniquemodcolors)){
	   thismod= uniquemodcolors[j]
	thisME <- MEs[,paste("ME",thismod,sep="")]
	# lm1 <- summary(lm(thisME~ID))
	boxplot(thisME~ID,col=c("blue","red","blue","red","blue","red"),main=paste(thismod, "Module ME ",sep=""), ylab="Module Eigengene Value")

}
dev.off()






####Trajectory CRND8 mouse

load('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/CRND8/Normalized_CRND8.rda')
datExpr=as.data.frame(t(normExpr_CRND8))

ensembl=read.csv('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/CRND8/GRCm38_hg38_homolog.csv')
ensembl=subset(ensembl,Homology.Type=='ortholog_one2one')

gnS=intersect(ensembl$Ensembl.Gene.ID,colnames(datExpr))
ensembl=ensembl[match(gnS,ensembl$Ensembl.Gene.ID),]
datExpr=datExpr[,match(gnS,colnames(datExpr))]
colnames(datExpr) <- ensembl$Human.Ensembl.Gene.ID
datExpr.CRND8=datExpr


library(WGCNA)

goodSamplesGenes(datExpr.CRND8) -> tmp ## one gene has no variance; removing the gene
rm.Gene=colnames(datExpr.CRND8)[which(tmp$goodGenes==FALSE)] ##remove this gene
a=setdiff(colnames(datExpr.CRND8),rm.Gene)
datExpr.CRND8=datExpr.CRND8[,match(a,colnames(datExpr.CRND8))]

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')
gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.CRND8))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.CRND8=datExpr.CRND8[,match(gnS,colnames(datExpr.CRND8))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


MEs=moduleEigengenes(datExpr.CRND8, colors = moduleColors.Mayo, nPC=1)$eigengenes
MEs=orderMEs(MEs)
uniquemodcolors=unique(moduleColors.Mayo)
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']

Tg=as.character(targets_CRND8$Genotype)
Tg[Tg=="-"] <- "Wt"
 Tg[Tg=="+"] <- "Tg"
Age <- as.numeric(targets_CRND8$Age_months)

ID=factor(paste(Tg,Age,sep="."))
ID=factor(ID,c('Wt.3','Tg.3','Wt.6','Tg.6','Wt.12','Tg.12','Wt.20','Tg.20'))

pdf("ME_Condition_CRND8mouse_MayoRef.pdf",height=4,width=4.5)
for (j in 1:length(uniquemodcolors)){
	   thismod= uniquemodcolors[j]
	thisME <- MEs[,paste("ME",thismod,sep="")]
	# lm1 <- summary(lm(thisME~ID))
	boxplot(thisME~ID,col=c("blue","red","blue","red","blue","red","blue","red"),main=paste(thismod, "Module ME ",sep=""), ylab="Module Eigengene Value")

}
dev.off()





####Trajectory

########## APP_Oligomers -- Entorhinal Cortex


load('/home/vivek//AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/APP_Oligomers/Regressed_EC_data.rda')
datExpr=as.data.frame(t(normExpr.reg))

ensembl=read.csv('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/CRND8/GRCm38_hg38_homolog.csv')
ensembl=subset(ensembl,Homology.Type=='ortholog_one2one')

gnS=intersect(ensembl$Ensembl.Gene.ID,colnames(datExpr))
ensembl=ensembl[match(gnS,ensembl$Ensembl.Gene.ID),]
datExpr=datExpr[,match(gnS,colnames(datExpr))]
colnames(datExpr) <- ensembl$Human.Ensembl.Gene.ID
datExpr.EC=datExpr


library(WGCNA)

goodSamplesGenes(datExpr.EC) -> tmp ## one gene has no variance; removing the gene
rm.Gene=colnames(datExpr.EC)[which(tmp$goodGenes==FALSE)] ##remove this gene
a=setdiff(colnames(datExpr.EC),rm.Gene)
datExpr.EC=datExpr.EC[,match(a,colnames(datExpr.EC))]



load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')
gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.EC))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.EC=datExpr.EC[,match(gnS,colnames(datExpr.EC))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


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





########## APP_Oligomers -- Dendate Gyrus

rm(list=ls())
load('/home/vivek//AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/APP_Oligomers/CombatCorrected_DG.rda')
datExpr=as.data.frame(t(datExpr.Combat))

ensembl=read.csv('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/CRND8/GRCm38_hg38_homolog.csv')
ensembl=subset(ensembl,Homology.Type=='ortholog_one2one')

gnS=intersect(ensembl$Ensembl.Gene.ID,colnames(datExpr))
ensembl=ensembl[match(gnS,ensembl$Ensembl.Gene.ID),]
datExpr=datExpr[,match(gnS,colnames(datExpr))]
colnames(datExpr) <- ensembl$Human.Ensembl.Gene.ID
datExpr.DG=datExpr


library(WGCNA)

goodSamplesGenes(datExpr.DG) -> tmp ## one gene has no variance; removing the gene
rm.Gene=colnames(datExpr.DG)[which(tmp$goodGenes==FALSE)] ##remove this gene
a=setdiff(colnames(datExpr.DG),rm.Gene)
datExpr.DG=datExpr.DG[,match(a,colnames(datExpr.DG))]



load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')
gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.DG))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.DG=datExpr.DG[,match(gnS,colnames(datExpr.DG))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]

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





###### Mouse rTG4510 and rTG4510

rm(list=ls())
load('/home/vivek//AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/Tg4510/Normalized_Expressed_genes.rda')
normExpr=normExpr[,targets$outlier==FALSE]
targets=subset(targets,outlier==FALSE)

normExpr_rTG4510=normExpr[,targets$Experiment=='MAPT_rTG4510']
targets_rTG4510=subset(targets,Experiment=='MAPT_rTG4510')

datExpr=as.data.frame(t(normExpr_rTG4510))

ensembl=read.csv('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/CRND8/GRCm38_hg38_homolog.csv')
ensembl=subset(ensembl,Homology.Type=='ortholog_one2one')

gnS=intersect(ensembl$Ensembl.Gene.ID,colnames(datExpr))
ensembl=ensembl[match(gnS,ensembl$Ensembl.Gene.ID),]
datExpr=datExpr[,match(gnS,colnames(datExpr))]
colnames(datExpr) <- ensembl$Human.Ensembl.Gene.ID
datExpr.rTG4510=datExpr


library(WGCNA)

goodSamplesGenes(datExpr.rTG4510) -> tmp ## one gene has no variance; removing the gene
rm.Gene=colnames(datExpr.rTG4510)[which(tmp$goodGenes==FALSE)] ##remove this gene
a=setdiff(colnames(datExpr.rTG4510),rm.Gene)
datExpr.rTG4510=datExpr.rTG4510[,match(a,colnames(datExpr.rTG4510))]



load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')
gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.rTG4510))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.rTG4510=datExpr.rTG4510[,match(gnS,colnames(datExpr.rTG4510))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]

####Trajectory

MEs=moduleEigengenes(datExpr.rTG4510, colors = moduleColors.Mayo, nPC=1)$eigengenes
MEs=orderMEs(MEs)
uniquemodcolors=unique(moduleColors.Mayo)
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']

Tg=as.character(targets_rTG4510$Genotype)
Tg[Tg=="-"] <- "Wt"
 Tg[Tg=="+"] <- "Tg"
Age <- as.numeric(targets_rTG4510$Age_months)

ID=factor(paste(Tg,Age,sep="."))
ID=factor(ID,c('Wt.2','Tg.2','Wt.6','Tg.6','Wt.12','Tg.12'))

pdf("ME_Condition_rTG4510mouse_MayoRef.pdf",height=4,width=4.5)
for (j in 1:length(uniquemodcolors)){
	   thismod= uniquemodcolors[j]
	thisME <- MEs[,paste("ME",thismod,sep="")]
	# lm1 <- summary(lm(thisME~ID))
	boxplot(thisME~ID,col=c("blue","red","blue","red","blue","red"),main=paste(thismod, "Module ME ",sep=""), ylab="Module Eigengene Value")

}
dev.off()




####rTg4510

rm(list=ls())
load('/home/vivek//AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/Tg4510/Normalized_Expressed_genes.rda')
normExpr=normExpr[,targets$outlier==FALSE]
targets=subset(targets,outlier==FALSE)

normExpr_rTG4510=normExpr[,targets$Experiment=='rTG4510']
targets_rTG4510=subset(targets,Experiment=='rTG4510')

datExpr=as.data.frame(t(normExpr_rTG4510))

ensembl=read.csv('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/CRND8/GRCm38_hg38_homolog.csv')
ensembl=subset(ensembl,Homology.Type=='ortholog_one2one')

gnS=intersect(ensembl$Ensembl.Gene.ID,colnames(datExpr))
ensembl=ensembl[match(gnS,ensembl$Ensembl.Gene.ID),]
datExpr=datExpr[,match(gnS,colnames(datExpr))]
colnames(datExpr) <- ensembl$Human.Ensembl.Gene.ID
datExpr.rTG4510=datExpr


library(WGCNA)

goodSamplesGenes(datExpr.rTG4510) -> tmp ## one gene has no variance; removing the gene
rm.Gene=colnames(datExpr.rTG4510)[which(tmp$goodGenes==FALSE)] ##remove this gene
a=setdiff(colnames(datExpr.rTG4510),rm.Gene)
datExpr.rTG4510=datExpr.rTG4510[,match(a,colnames(datExpr.rTG4510))]





load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')
gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.rTG4510))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.rTG4510=datExpr.rTG4510[,match(gnS,colnames(datExpr.rTG4510))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]

####Trajectory

MEs=moduleEigengenes(datExpr.rTG4510, colors = moduleColors.Mayo, nPC=1)$eigengenes
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
