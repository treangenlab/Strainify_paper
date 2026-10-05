#!/bin/bash

# 4-strain mock dataset (not weighted)
./strainify \
    --genome_folder ../Strainify_paper/4_strain_mock_community/genomes \
    --fastq_folder ../Strainify_paper/4_strain_mock_community/fastq \
    --outdir ../Strainify_paper/4_strain_mock_community/strainify/results \
    --max_cpus 12