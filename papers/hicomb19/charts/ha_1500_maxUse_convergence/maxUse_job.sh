#!/bin/bash
#PBS -j oe
#PBS -l "nodes=1:ppn=1"
#PBS -l "mem=4gb"
#PBS -N "maxUse_0.10"

maxUse=0.10

cd $PBS_O_WORKDIR

/usr/bin/time -v ~/research/peace/src/peace --print-progress --fastaFile /shared/raodm/sutharzan/datasets/HA_1500_rnd.fa --clusterMaker mst --maxUse ${maxUse}

