#!/bin/bash
#$ -N transport_data
#$ -q free64
#$ -m beas
#$ -ckpt restart

rsync -r smorabit@128.97.126.96:~/ ~/swaruplab/smorabit/nessie
