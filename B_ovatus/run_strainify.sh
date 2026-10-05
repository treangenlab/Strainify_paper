#!/bin/bash

./strainify \
    --genome_folder ../Strainify_paper/B_ovatus/all_fastas \
    --fastq_folder ../Strainify_paper/B_ovatus/fastqs \
    --outdir ../Strainify_paper/B_ovatus/results \
    --weight_by_entropy \
    --parsnp_flags " --skip-ani-filter" \
    --max_cpus 12
