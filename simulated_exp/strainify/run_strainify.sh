# !/bin/bash

# 4-strain simulated dataset (not weighted)
./strainify run --cores 12 --configfile ../Strainify_paper/simulated_exp/strainify/config/4_strain_ecoli_simulated.yaml

# 30-strain simulated dataset (not weighted)
./strainify run --cores 12 --configfile ../Strainify_paper/simulated_exp/strainify/config/30_strain_simulated_cdiff.yaml
./strainify run --cores 12 --configfile ../Strainify_paper/simulated_exp/strainify/config/30_strain_simulated_ecoli.yaml
./strainify run --cores 12 --configfile ../Strainify_paper/simulated_exp/strainify/config/30_strain_simulated_mtuberculosis.yaml
./strainify run --cores 12 --configfile ../Strainify_paper/simulated_exp/strainify/config/30_strain_simulated_sepidermidis.yaml
./strainify run --cores 12 --configfile ../Strainify_paper/simulated_exp/strainify/config/30_strain_simulated_cacnes.yaml

# 30-strain simulated dataset (weighted)
./strainify run --cores 12 --configfile ../Strainify_paper/simulated_exp/strainify/config/30_strain_simulated_cdiff.yaml --config weight_by_entropy=true output_dir=../Strainify_paper/simulated_exp/strainify/30_strain_simulated/cdiff/weighted
./strainify run --cores 12 --configfile ../Strainify_paper/simulated_exp/strainify/config/30_strain_simulated_ecoli.yaml --config weight_by_entropy=true output_dir=../Strainify_paper/simulated_exp/strainify/30_strain_simulated/ecoli/weighted
./strainify run --cores 12 --configfile ../Strainify_paper/simulated_exp/strainify/config/30_strain_simulated_mtuberculosis.yaml --config weight_by_entropy=true output_dir=../Strainify_paper/simulated_exp/strainify/30_strain_simulated/mtuberculosis/weighted
./strainify run --cores 12 --configfile ../Strainify_paper/simulated_exp/strainify/config/30_strain_simulated_sepidermidis.yaml --config weight_by_entropy=true output_dir=../Strainify_paper/simulated_exp/strainify/30_strain_simulated/sepidermidis/weighted
./strainify run --cores 12 --configfile ../Strainify_paper/simulated_exp/strainify/config/30_strain_simulated_cacnes.yaml --config weight_by_entropy=true output_dir=../Strainify_paper/simulated_exp/strainify/30_strain_simulated/cacnes/weighted