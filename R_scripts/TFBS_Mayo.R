
uniquemodcolors <- as.character(unique(geneInfo.Mayo$Initially.Assigned.Module.Color))

for(i in 1:length(uniquemodcolors)){
	thismod= uniquemodcolors[i]
	thisInfo=geneInfo.Mayo[geneInfo.Mayo$Initially.Assigned.Module.Color==thismod,c(2,1)] ##1  =Ensembl.ID, 2= geneName
  thisInfo$Sno=seq(1:nrow(thisInfo))
  thisInfo=thisInfo[,c(3,1,2)]
	colnames(thisInfo) <- c("Sno","Genename","EnsemblID")
	write.table(thisInfo,file=paste("./",thismod,"_TFBS_user.file.csv",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep=",")
}
write.csv(uniquemodcolors,'uniquemodcolors.csv')


my.modules <- read.csv('uniquemodcolors.csv')
my.modules <- as.character(my.modules[,2])

#needed to install these packages in order to get this code to work:
install.packages("seqinr")
BiocManager::install("pcaMethods")
install.packages("MEET")
install.packages("Cairo")

for(i in 1:length(my.modules)){
  mymods <- my.modules[i]
  TFBSenrich (user.file=paste(mymods,"_TFBS_user.file.csv",sep=""),
      TF.db=system.file("data/JASPAR_CLOVER",package="RegFacEnc"),
      TF.nome=system.file("data/JASPAR_NOMENCLATURE_TABLE",package="RegFacEnc"),
      db.seq=system.file("data/UP1000_Protien_Coding_HUMANS_unique.fasta",package="RegFacEnc"),
      cpg.seq=system.file("data/HUMAN_CpG.fa",package="RegFacEnc"),
      chro.seq=system.file("data/HUMAN_chr20.fa",package="RegFacEnc"),option= "-t",
      pval = "0.05", species = "Human", TF_motifs="JASPAR", BF_Type="Protein Coding")
}


#TFBSenrich (user.file="lightgreen_TFBS_user.file.csv", TF.db=system.file("data/JASPAR_CLOVER",package="RegFacEnc"), TF.nome=system.file("data/JASPAR_NOMENCLATURE_TABLE",package="RegFacEnc"),                      db.seq=system.file("data/UP1000_Protien_Coding_HUMANS_unique.fasta",package="RegFacEnc"), cpg.seq=system.file("data/HUMAN_CpG.fa",package="RegFacEnc"),
#            chro.seq=system.file("data/HUMAN_chr20.fa",package="RegFacEnc"),option= "-t",
#            pval = "0.05", species = "Human", TF_motifs="JASPAR", BF_Type="Protein Coding")
