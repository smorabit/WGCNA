#!/bin/bash
#$ -N mayo_go
#$ -q free64
#$ -ckpt restart

# module load R/3.4.1
# Rscript ~/bin/mayo_go_analysis.R
#
# codedir <- "~/bin/GO-Elite_v.1.2.5-Py"
# pathname <- "/pub/smorabit/WGCNA/geneInfo"
# nperm=10000
# system(paste("nohup python ",codedir,"/GO_Elite.py --species Hs --mod Ensembl --permutations ",
#              nperm,"  --method \"z-score\" --zscore 1.96 --pval 0.01 --num 5 --input ",pathname,
#              "/input --denom ",pathname,"/background --output ",pathname,"/output &",sep=""))

python ~/bin/GO-Elite_v.1.2.5-Py/GO_Elite.py --species Hs --mod Ensembl --permutations 10000 --method z-score --zscore 1.96 --pval 0.01 --num 5 --input /pub/smorabit/WGCNA/geneInfo/input --denom /pub/smorabit/WGCNA/geneInfo/background --output /pub/smorabit/WGCNA/geneInfo/output
