#!/bin/bash
#$ -N map2transcriptome 
#$ -V
#$ -cwd
#$ -pe smp 32
#$ -j y
#$ -l mem_requested=16G,tmp_requested=120G
#$ -e /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files
#$ -o /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files

# ARG1: sample name

reference=/directflow/KCCGGenometechTemp/projects/PTSD_mods/reference_files/gencode.vM25.transcripts.renamed.fa
dir=/directflow/KCCGGenometechTemp/projects/PTSD_mods/rna_modifications/nanocompare/alignments
fastq=/directflow/KCCGGenometechTemp/projects/PTSD_mods/raw_data/fastqs/${1}.fastq
bam=${dir}/${1}.bam
ncpus=32

minimap2 -ax map-ont -t ${ncpus} -L ${reference} ${fastq} | samtools view -bh -F 2324 -q 10 | samtools sort  -@ ${ncpus} -O bam > ${bam}
samtools index ${bam}
