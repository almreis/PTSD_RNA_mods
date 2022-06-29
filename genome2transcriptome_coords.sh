#!/bin/bash
#$ -N genome2transcriptome_coords 
#$ -V
#$ -cwd
#$ -pe smp 16
#$ -j y
#$ -l mem_requested=8G,tmp_requested=120G
#$ -e /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files
#$ -o /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files

# ARG1: prefix of bed file to convert

Rscript genome2transcriptome_coords.R ${1}
