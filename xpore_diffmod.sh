#!/bin/bash
#$ -cwd
#$ -V
#$ -N xpore_diffmod 
#$ -q long.q
#$ -e /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files
#$ -o /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files
#$ -S /bin/bash
#$ -j y
#$ -pe smp 8
#$ -l mem_requested=24G
#$ -l tmp_requested=200G

# ARG1: male or female

CONFIG=/directflow/KCCGGenometechTemp/projects/PTSD_mods/rna_modifications/xpore/${1}/config.yml

xpore diffmod --config ${CONFIG}

