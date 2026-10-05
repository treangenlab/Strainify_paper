# Please change the paths to your local setup
./strainify \
    --genome_folder /home/Users/rl152/Strainify/UMB/UMB18_genomes \
    --fastq_folder /dodo/rl152/UMB/PRJNA400628_UMB18/UMB18_selected_fastqs \
    --outdir /home/Users/rl152/Strainify/UMB/UMB18_results_weighted \
    --threads 12 \
    --weight_by_entropy