#!/bin/bash

#load appropriate module files:
module load samtools/1.9
module load bedtools/2.25.0

# navigate to appropriate directory on HPC:
cd /pub/smorabit/ATAC

# get fasta sequences within ATAC peaks:
for line in $(cat control_ids.txt)
do
  name=$(basename peaks/*$line*)
  echo $name
  # zcat peaks/*$line* | bedtools getfasta -fi /pub/smorabit/resources/hg19.fa -bed - -fo peak_sequences/$name.fa

  # if [[ name =~ *.narrowPeak.gz]]
  # then
  #   echo $line
  # fi

done

findMotifsGenome.pl test_peak.bed ../resources/hg19.fa test_homer/ -size 200
