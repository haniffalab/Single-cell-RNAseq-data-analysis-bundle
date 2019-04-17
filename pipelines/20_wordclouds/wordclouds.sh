#!/bin/bash

#$ -cwd
#$ -N wordclouds
#$ -V
#$ -l h_rt=23:59:59
#$ -l h_vmem=100G

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

Rscript wordclouds.R $1

echo "End on `date`"