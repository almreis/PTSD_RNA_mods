#!/bin/bash
#$ -cwd
#$ -V
#$ -N f5c-eventalign
#$ -e /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files
#$ -o /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files
#$ -S /bin/bash
#$ -j y
#$ -pe smp 32
#$ -l mem_requested=4G
#$ -l h_vmem=4G
#$ -l tmp_requested=150G
#$ -l nvgpu=1

# ARG1: Sample prefix

source $DF/software/nanocompare/bin/activate

dir=/directflow/KCCGGenometechTemp/projects/PTSD_mods
input=${dir}/rna_modifications/nanocompare/alignments/${1}.eventalign.tsv
output=${dir}/rna_modifications/nanocompare/alignments/${1}.eventalign.collapsed

nanocompore eventalign_collapse -t 32 -i ${input} -o ${output}

