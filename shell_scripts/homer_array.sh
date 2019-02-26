#!/bin/bash
#$ -N motif_enrichment
#$ -q free64
#$ -t 1-138

module load homer/4.7

findMotifsGenome.pl peaks_control/*.$SGE_TASK_ID /pub/smorabit/resources/hg19.fa homer/*.$SGE_TASK_ID
