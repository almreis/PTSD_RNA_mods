#!/bin/bash
#$ -N f5c_index
#$ -e /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files
#$ -o /directflow/KCCGGenometechTemp/projects/PTSD_mods/log_files
#$ -q short.q
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 32 
#$ -j y
#$ -l mem_requested=8G,tmp_requested=60G

# ARG1: sequencing RunID
# ARG2: sample name

module load centos6.10/hasgam/f5c-cpu 

dir=/directflow/KCCGGenometechTemp/projects/PTSD_mods
fast5=${dir}/raw_data/${1}/${1}_reads.tar

location=`tar -tf ${fast5} | grep -m 1 fast5_pass | rev | cut -d/ -f2- | rev`

fast5=${dir}/raw_data/${1}/${location}
fastq=/directflow/KCCGGenometechTemp/projects/PTSD_mods/raw_data/fastqs/${2}.fastq

echo ${fast5}
echo ${fastq}

f5c index --iop 64 -t 32 -d ${fast5} ${fastq} 
