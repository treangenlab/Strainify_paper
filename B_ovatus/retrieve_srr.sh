#!/bin/bash
# Retrieve and download FASTQ files for SAMN range

for i in $(seq 11950081 11950286); do
  esearch -db biosample -query SAMN${i} | \
    elink -target sra | \
    efetch -format runinfo | \
    cut -d',' -f1 | grep ^SRR
done > srr_list.txt

# Prefetch and download
for srr in $(cat srr_list.txt); do
  prefetch $srr
  fasterq-dump $srr --split-files --threads 4
done

