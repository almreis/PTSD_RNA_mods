#!/bin/bash
#$ -cwd
#$ -V
#$ -N xpore 
#$ -q long.q
#$ -e /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files
#$ -o /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files
#$ -S /bin/bash
#$ -j y
#$ -pe smp 8
#$ -l mem_requested=24G
#$ -l tmp_requested=200G

# ARG1: eventalign file
# ARG2: male or female
# ARG3: sample prefix

dir=/directflow/KCCGGenometechTemp/projects/PTSD_mods
input=${dir}/rna_modifications/nanocompare/alignments
EVENTALIGN=${input}/${1}
OUT_DIR=${dir}/rna_modifications/xpore/${2}/${3}

xpore dataprep --eventalign ${EVENTALIGN} --out_dir ${OUT_DIR}

