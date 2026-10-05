# !/bin/bash

# 4-strain simulated dataset (not weighted)
./strainify \
    --genome_folder ../Strainify_paper/simulated_exp/4_strain_simulated/genomes \
    --fastq_folder ../Strainify_paper/simulated_exp/4_strain_simulated/fastq \
    --outdir ../Strainify_paper/simulated_exp/strainify/results \
    --max_cpus 12

# 30-strain simulated dataset (not weighted)
./strainify \
    --genome_folder ../Strainify_paper/simulated_exp/30_strain_simulated/cacnes/cacnes_downloads \
    --fastq_folder ../Strainify_paper/simulated_exp/30_strain_simulated/cacnes/fastq \
    --outdir ../Strainify_paper/simulated_exp/strainify/30_strain_simulated/cacnes/not_weighted \
    --max_cpus 12
./strainify \
    --genome_folder ../Strainify_paper/simulated_exp/30_strain_simulated/cdiff/cdiff_downloads \
    --fastq_folder ../Strainify_paper/simulated_exp/30_strain_simulated/cdiff/fastq \
    --outdir ../Strainify_paper/simulated_exp/strainify/30_strain_simulated/cdiff/not_weighted \
    --max_cpus 12
./strainify \
    --genome_folder ../Strainify_paper/simulated_exp/30_strain_simulated/ecoli/ecoli_downloads \
    --fastq_folder ../Strainify_paper/simulated_exp/30_strain_simulated/ecoli/fastq \
    --outdir ../Strainify_paper/simulated_exp/strainify/30_strain_simulated/ecoli/not_weighted \
    --max_cpus 12
./strainify \
    --genome_folder ../Strainify_paper/simulated_exp/30_strain_simulated/mtuberculosis/mtuberculosis_downloads \
    --fastq_folder ../Strainify_paper/simulated_exp/30_strain_simulated/mtuberculosis/fastq \
    --outdir ../Strainify_paper/simulated_exp/strainify/30_strain_simulated/mtuberculosis/not_weighted \
    --max_cpus 12
./strainify \
    --genome_folder ../Strainify_paper/simulated_exp/30_strain_simulated/sepidermidis/sepidermidis_downloads \
    --fastq_folder ../Strainify_paper/simulated_exp/30_strain_simulated/sepidermidis/fastq \
    --outdir ../Strainify_paper/simulated_exp/strainify/30_strain_simulated/sepidermidis/not_weighted \
    --max_cpus 12

# 30-strain simulated dataset (weighted)
./strainify \
    --genome_folder ../Strainify_paper/simulated_exp/30_strain_simulated/cacnes/cacnes_downloads \
    --fastq_folder ../Strainify_paper/simulated_exp/30_strain_simulated/cacnes/fastq \
    --outdir ../Strainify_paper/simulated_exp/strainify/30_strain_simulated/cacnes/weighted \
    --weight_by_entropy \
    --max_cpus 12
./strainify \
    --genome_folder ../Strainify_paper/simulated_exp/30_strain_simulated/cdiff/cdiff_downloads \
    --fastq_folder ../Strainify_paper/simulated_exp/30_strain_simulated/cdiff/fastq \
    --outdir ../Strainify_paper/simulated_exp/strainify/30_strain_simulated/cdiff/weighted \
    --weight_by_entropy \
    --max_cpus 12
./strainify \
    --genome_folder ../Strainify_paper/simulated_exp/30_strain_simulated/ecoli/ecoli_downloads \
    --fastq_folder ../Strainify_paper/simulated_exp/30_strain_simulated/ecoli/fastq \
    --outdir ../Strainify_paper/simulated_exp/strainify/30_strain_simulated/ecoli/weighted \
    --weight_by_entropy \
    --max_cpus 12
./strainify \
    --genome_folder ../Strainify_paper/simulated_exp/30_strain_simulated/mtuberculosis/mtuberculosis_downloads \
    --fastq_folder ../Strainify_paper/simulated_exp/30_strain_simulated/mtuberculosis/fastq \
    --outdir ../Strainify_paper/simulated_exp/strainify/30_strain_simulated/mtuberculosis/weighted \
    --weight_by_entropy \
    --max_cpus 12
./strainify \
    --genome_folder ../Strainify_paper/simulated_exp/30_strain_simulated/sepidermidis/sepidermidis_downloads \
    --fastq_folder ../Strainify_paper/simulated_exp/30_strain_simulated/sepidermidis/fastq \
    --outdir ../Strainify_paper/simulated_exp/strainify/30_strain_simulated/sepidermidis/weighted \
    --weight_by_entropy \
    --max_cpus 12