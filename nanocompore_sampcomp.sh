#!/bin/bash
#$ -cwd
#$ -V
#$ -N nanocompare_sampcomp 
#$ -q long.q
#$ -e /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files
#$ -o /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files
#$ -S /bin/bash
#$ -j y
#$ -pe smp 8
#$ -l mem_requested=24G
#$ -l tmp_requested=200G

# ARG1: male or female

source $DF/software/nanocompare/bin/activate

dir=/directflow/KCCGGenometechTemp/projects/PTSD_mods
input=${dir}/rna_modifications/nanocompare/alignments
reference=${dir}/reference_files/gencode.vM25.transcripts.renamed.fa
c1=${input}/control.${1}.1.eventalign.collapsed/out_eventalign_collapse.tsv
c2=${input}/control.${1}.2.eventalign.collapsed/out_eventalign_collapse.tsv
c3=${input}/control.${1}.5.eventalign.collapsed/out_eventalign_collapse.tsv
s1=${input}/stress.${1}.1.eventalign.collapsed/out_eventalign_collapse.tsv
s2=${input}/stress.${1}.4.eventalign.collapsed/out_eventalign_collapse.tsv
s3=${input}/stress.${1}.5.eventalign.collapsed/out_eventalign_collapse.tsv
output=${dir}/rna_modifications/nanocompare/sampcomp

nanocompore sampcomp --min_coverage 30 --downsample_high_coverage 1000 --overwrite --nthreads 8 --file_list1 ${c1},${c2},${c3} --file_list2 ${s1},${s2},${s3} --label1 control --label2 stress --fasta ${reference} --outpath ${output}

