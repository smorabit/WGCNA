


#########Zhang et al HBRTC Prefrontal cortex

###
library(WGCNA)
load('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/Zhang_HBRTC_PFC/Zhang_ENSG_PFC.rda')
load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda') #this is where the file will be (no file now)

setLabels=c("Mayo.AD","ZhangPFC")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.Zhang)) #intersect genes between zhang and mayo
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.Zhang=datExpr.Zhang[,match(gnS,colnames(datExpr.Zhang))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]

multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), ZhangPFC =list(data= datExpr.Zhang))
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_Zhang_Mayo_AD.rda"
save(list=ls(),file=filename)


################Mt Sinai STG region
library(WGCNA)

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/MSSM/MSSM_FP_STG_PHG_IFG_Expression_metaData.rda')

targets.MSSM.STG=targets.MSSM[targets.MSSM$Region=="STG",]
normExpr.MSSM.STG=normExpr.MSSM[,targets.MSSM$Region=="STG"]
datExpr.MSSM.STG=as.data.frame(t(normExpr.MSSM.STG))

##Similarly separate other regions and perform separate module preservation analysis (not together)

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')


setLabels=c("Mayo.AD","MtSinai")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.MSSM.STG))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.MSSM.STG=datExpr.MSSM.STG[,match(gnS,colnames(datExpr.MSSM.STG))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), MtSinai =list(data= datExpr.MSSM.STG))
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_MSSM.STG_Mayo_AD.rda"
save(list=ls(),file=filename)



################Mt Sinai FP region
library(WGCNA)

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/MSSM/MSSM_FP_STG_PHG_IFG_Expression_metaData.rda')

targets.MSSM.FP=targets.MSSM[targets.MSSM$Region=="FP",]
normExpr.MSSM.FP=normExpr.MSSM[,targets.MSSM$Region=="FP"]
datExpr.MSSM.FP=as.data.frame(t(normExpr.MSSM.FP))


load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')


setLabels=c("Mayo.AD","MtSinai")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.MSSM.FP))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.MSSM.FP=datExpr.MSSM.FP[,match(gnS,colnames(datExpr.MSSM.FP))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), MtSinai =list(data= datExpr.MSSM.FP))
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_MSSM.FP_Mayo_AD.rda"
save(list=ls(),file=filename)




################Mt Sinai IFG region
library(WGCNA)

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/MSSM/MSSM_FP_STG_PHG_IFG_Expression_metaData.rda')

targets.MSSM.IFG=targets.MSSM[targets.MSSM$Region=="IFG",]
normExpr.MSSM.IFG=normExpr.MSSM[,targets.MSSM$Region=="IFG"]
datExpr.MSSM.IFG=as.data.frame(t(normExpr.MSSM.IFG))

##Similarly separate other regions and perform separate module preservation analysis (not together)

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')


setLabels=c("Mayo.AD","MtSinai")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.MSSM.IFG))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.MSSM.IFG=datExpr.MSSM.IFG[,match(gnS,colnames(datExpr.MSSM.IFG))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), MtSinai =list(data= datExpr.MSSM.IFG))
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_MSSM.IFG_Mayo_AD.rda"
save(list=ls(),file=filename)



################Mt Sinai PHG region
library(WGCNA)

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/MSSM/MSSM_FP_STG_PHG_IFG_Expression_metaData.rda')

targets.MSSM.PHG=targets.MSSM[targets.MSSM$Region=="PHG",]
normExpr.MSSM.PHG=normExpr.MSSM[,targets.MSSM$Region=="PHG"]
datExpr.MSSM.PHG=as.data.frame(t(normExpr.MSSM.STG))

##Similarly separate other regions and perform separate module preservation analysis (not together)

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')


setLabels=c("Mayo.AD","MtSinai")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.MSSM.PHG))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.MSSM.PHG=datExpr.MSSM.PHG[,match(gnS,colnames(datExpr.MSSM.PHG))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), MtSinai =list(data= datExpr.MSSM.PHG))
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_MSSM.PHG_Mayo_AD.rda"
save(list=ls(),file=filename)




#####################
#########  Myers et al AD Frontal Cortex
library(WGCNA)

load('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/Myers_AD/Myers_frontal_ENSG.rda')
load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')

setLabels=c("Mayo.AD","Myers_FrontalCTX")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.F))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.F=datExpr.F[,match(gnS,colnames(datExpr.F))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), Myers_FrontalCTX =list(data= datExpr.F))
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_MyersFC_Mayo_AD.rda"
save(list=ls(),file=filename)

####

##### Myers et al AD Temporal Cortex

library(WGCNA)

load('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/Myers_AD/Myers_Temporal_ENSG.rda')
load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')

setLabels=c("Mayo.AD","Myers_TemporalCTX")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.T))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.T=datExpr.T[,match(gnS,colnames(datExpr.T))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), Myers_TemporalCTX =list(data= datExpr.T))
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_MyersTC_Mayo_AD.rda"
save(list=ls(),file=filename)

################
###########################


#############Berchtold data

library(WGCNA)

load('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/Berchtold_AD/Berchtoldlog2QN_RMA_outlierRemoved_combat_unfiltered_collapseRows_22647probes_techRegressed.RData')

targets.F=datPheno[datPheno$Region=='superior frontal gyrus',]
targets.T=datPheno[datPheno$Region=='postcentral gyrus',]

datExpr=as.data.frame(t(datExpr))

datExpr.F=datExpr[datPheno$Region=='superior frontal gyrus',]
datExpr.T=datExpr[datPheno$Region=='postcentral gyrus',]


########Berchtold data Frontal Cortex
library(WGCNA)

load('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/Berchtold_AD/Berchtold_FrontalCTX.rda')
load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')

setLabels=c("Mayo.AD","Berchtold_FrontalCTX")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.F))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.F=datExpr.F[,match(gnS,colnames(datExpr.F))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), Berchtold_FrontalCTX =list(data= datExpr.F))
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_BerchtoldFC_Mayo_AD.rda"
save(list=ls(),file=filename)

####Berchtold data Temporal Cortex

library(WGCNA)

load('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/Berchtold_AD/Berchtold_TemporalCTX.rda')
load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')

setLabels=c("Mayo.AD","Berchtold_TemporalCTX")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.T))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.T=datExpr.T[,match(gnS,colnames(datExpr.T))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), Berchtold_TemporalCTX =list(data= datExpr.T))
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_BerchtoldTC_Mayo_AD.rda"
save(list=ls(),file=filename)


#######
#####################################
#####BLSA Emory Proteomics Dataset
library(WGCNA)

rm(list=ls())
load('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/Emory_BLSA_Proteomics/BLSA_ProteomicsData_Final.rda')

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')

ensembl=read.csv('/home/vivek/AMP_AD/Mayo/Analysis/hg38_77_GeneSymbol.csv')
ensembl.subset=ensembl[match(colnames(datExpr.Ref),ensembl$Ensembl.Gene.ID),]



setLabels=c("Mayo.AD","BLSA_Proteomics")
shortLabels=setLabels
nSets=2

gnS=intersect(ensembl.subset$Associated.Gene.Name,geneId)## geneId from BLSA data, ensembl.subset from Discovery Set

datExpr.Ref1=datExpr.Ref[,match(gnS,ensembl.subset$Associated.Gene.Name)]
datExpr.BLSA.subset=datExpr.BLSA[,match(gnS,geneId)]

moduleColors.Mayo=moduleColors.Mayo[match(gnS,ensembl.subset$Associated.Gene.Name)]
colnames(datExpr.BLSA.subset)=colnames(datExpr.Ref1)


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), BLSA_Proteomics =list(data= datExpr.BLSA.subset))
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_BLSAEmory_Mayo_AD.rda"
save(list=ls(),file=filename)


######### Reference EMORY Proteomics
library(WGCNA)

load('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/Emory_BLSA_Proteomics/BLSA_ProteomicsData_Final.rda')

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')

ensembl=read.csv('~/AMP_AD/Mayo/Analysis/hg38_77_GeneSymbol.csv')
ensembl.subset=ensembl[match(colnames(datExpr.Ref),ensembl$Ensembl.Gene.ID),]



setLabels=c("BLSA_Proteomics","Mayo.AD")
shortLabels=setLabels
nSets=2

gnS=intersect(ensembl.subset$Associated.Gene.Name,geneId)## geneId from BLSA data, ensembl.subset from Discovery Set

datExpr.Ref1=datExpr.Ref[,match(gnS,ensembl.subset$Associated.Gene.Name)]
datExpr.BLSA.subset=datExpr.BLSA[,match(gnS,geneId)]
colnames(datExpr.BLSA.subset)=colnames(datExpr.Ref1)

geneInfo.BLSA=geneInfo.BLSA[match(gnS,geneId),]
moduleColors.BLSA=geneInfo.BLSA$net.colors


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(BLSA_Proteomics =list(data= datExpr.BLSA.subset), Mayo.AD =list(data= datExpr.Ref1))
checkSets(multiExpr) # check data size
multiColors=list(BLSA_Proteomics = moduleColors.BLSA)
filename="mp_Mayo_AD_BLSAEmoryRef.rda"
save(list=ls(),file=filename)


############ROSMAP Dataset
library(WGCNA)

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/ROSMAP/ROSMAP_DLPFC_expression_metaData.rda')

datExpr.ROSMAP=as.data.frame(t(normExpr.ROSMAP))

load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')


setLabels=c("Mayo.AD","ROSMAP")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.ROSMAP))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.ROSMAP=datExpr.ROSMAP[,match(gnS,colnames(datExpr.ROSMAP))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), ROSMAP =list(data= datExpr.ROSMAP))
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_ROSMAP_Mayo_AD.rda"
save(list=ls(),file=filename)








#####Mouse preservation -- APP

library(WGCNA)

load('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/CRND8/Normalized_APP.rda')
datExpr=as.data.frame(t(normExpr_APP))

ensembl=read.csv('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/CRND8/GRCm38_hg38_homolog.csv')
ensembl=subset(ensembl,Homology.Type=='ortholog_one2one')

gnS=intersect(ensembl$Ensembl.Gene.ID,colnames(datExpr))
ensembl=ensembl[match(gnS,ensembl$Ensembl.Gene.ID),]
datExpr=datExpr[,match(gnS,colnames(datExpr))]
colnames(datExpr) <- ensembl$Human.Ensembl.Gene.ID
datExpr.APP=datExpr



goodSamplesGenes(datExpr.APP) -> tmp ## one gene has no variance; removing the gene
rm.Gene=colnames(datExpr.APP)[which(tmp$goodGenes==FALSE)] ##remove this gene
a=setdiff(colnames(datExpr.APP),rm.Gene)
datExpr.APP=datExpr.APP[,match(a,colnames(datExpr.APP))]


load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')

######
setLabels=c("Mayo.AD","APPmouse")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.APP))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]


datExpr.APP=datExpr.APP[,match(gnS,colnames(datExpr.APP))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]



multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), APPmouse =list(data= datExpr.APP))
library(WGCNA)
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_APPmouse_Mayo_AD.rda"
save(list=ls(),file=filename)





#### CRND8 mouse
library(WGCNA)

rm(list=ls())

load('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/CRND8/Normalized_CRND8.rda')
datExpr=as.data.frame(t(normExpr_CRND8))

ensembl=read.csv('/home/vivek/AMP_AD/Mayo/Analysis/Step05_ModulePreserv/AD/MouseModels/CRND8/GRCm38_hg38_homolog.csv')
ensembl=subset(ensembl,Homology.Type=='ortholog_one2one')

gnS=intersect(ensembl$Ensembl.Gene.ID,colnames(datExpr))
ensembl=ensembl[match(gnS,ensembl$Ensembl.Gene.ID),]
datExpr=datExpr[,match(gnS,colnames(datExpr))]
colnames(datExpr) <- ensembl$Human.Ensembl.Gene.ID
datExpr.CRND8=datExpr



goodSamplesGenes(datExpr.CRND8) -> tmp ## one gene has no variance; removing the gene
rm.Gene=colnames(datExpr.CRND8)[which(tmp$goodGenes==FALSE)] ##remove this gene
a=setdiff(colnames(datExpr.CRND8),rm.Gene)
datExpr.CRND8=datExpr.CRND8[,match(a,colnames(datExpr.CRND8))]



load('/home/vivek/AMP_AD/AMP_AD_NormalizedData_rWGCNAs/Expression_Jan2019/Mayo/rWGCNA/rWGCNA_Mayo_ForPreservation.rda')

######
setLabels=c("Mayo.AD","CRND8mouse")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.CRND8))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]
datExpr.CRND8=datExpr.CRND8[,match(gnS,colnames(datExpr.CRND8))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]


multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), CRND8mouse =list(data= datExpr.CRND8))
library(WGCNA)
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_CRND8mouse_Mayo_AD.rda"
save(list=ls(),file=filename)



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

######
setLabels=c("Mayo.AD","Oligomer_EC_mouse")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.EC))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]


datExpr.EC=datExpr.EC[,match(gnS,colnames(datExpr.EC))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]



multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), Oligomer_EC_mouse =list(data= datExpr.EC))
library(WGCNA)
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_EC_Oligomer_mouse_Mayo_AD.rda"
save(list=ls(),file=filename)



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

######
setLabels=c("Mayo.AD","Oligomer_DG_mouse")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.DG))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]


datExpr.DG=datExpr.DG[,match(gnS,colnames(datExpr.DG))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]



multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), Oligomer_DG_mouse =list(data= datExpr.DG))
library(WGCNA)
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_DG_Oligomer_mouse_Mayo_AD.rda"
save(list=ls(),file=filename)






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

######
setLabels=c("Mayo.AD","rTG4510mouse")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.rTG4510))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]


datExpr.rTG4510=datExpr.rTG4510[,match(gnS,colnames(datExpr.rTG4510))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]



multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), rTG4510mouse =list(data= datExpr.rTG4510))
library(WGCNA)
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_rTG4510mouse_Mayo_AD.rda"
save(list=ls(),file=filename)




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

######
setLabels=c("Mayo.AD","rTG4510mouse")
shortLabels=setLabels
nSets=2

gnS=intersect(colnames(datExpr.Ref),colnames(datExpr.rTG4510))
datExpr.Ref1=datExpr.Ref[,match(gnS,colnames(datExpr.Ref))]


datExpr.rTG4510=datExpr.rTG4510[,match(gnS,colnames(datExpr.rTG4510))]
moduleColors.Mayo=moduleColors.Mayo[match(gnS,colnames(datExpr.Ref))]



multiExpr =vector(mode="list",length=nSets)
multiExpr =list(Mayo.AD =list(data= datExpr.Ref1), rTG4510mouse =list(data= datExpr.rTG4510))
library(WGCNA)
checkSets(multiExpr) # check data size
multiColors=list(Mayo.AD = moduleColors.Mayo)
filename="mp_rTG4510mouse_Mayo_AD.rda"
save(list=ls(),file=filename)
