#!/bin/bash -l

#$ -S /bin/bash

#$ -P camplab

#$ -j y

#$ -l h_rt=20:00:00

#$ -l mem_per_core=8G

#$ -pe omp 15

#$ -o module_heatmap.out

#$ -e module_heatmap.err

#$ -N report

echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="

module load pandoc
module load R/4.1.2
Rscript supplementary_module_heatmap.R

echo "Job Finished"
echo "Finish date : $(date)"
