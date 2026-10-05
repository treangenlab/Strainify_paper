# CAMI II dataset

We followed the ChronoStrain paper's method for the CAMI II analyses. The only change we made was using the core genome only for all reference strains. 

## 1. Download dataset
Please follow the ChronoStrain paper's instructions to download the CAMI II strain-madness dataset (https://github.com/gibsonlab/chronostrain_cami)

## 2. Running Strainify

Run Strainify with the original reference genomes from the CAMI dataset (not the core genome FASTAs). Strainify will first run Parsnp to extract the core genomes. This will produce the same MAF files for the next step. Strainify will then infer relative abundances based on variant frequencies in the core genome only. See the `Strainify` folder for the relevant commands. 

## 3. Extracting core genomes
Run the `extract_core.py` script with each of the `parsnp.maf` file in the Strainify output (inside the `parsnp_results` subdirectory). These files are the Parsnp alignment results and there is a separate alignment for each of the five species. The resulting FASTA files containing the core genomes can also be found in the `fastas` folder. These will be used as reference genomes for ChronoStrain and StrainGE.

## 4. Running ChronoStrain and StrainGE

### ChronoStrain
Update the `GOLD_STANDARD_INDEX` in the original `settings_global.sh` script (ChronoStrain paper's version) to use the core genome FASTAs. See the examples in the `ChronoStrain` folder here. Make sure to change the paths in the `core_genome_index.tsv` file to match your local setup. Follow the `chronostrain_cami` repository for instructions on database construction and inference. For database construction, use the `gold_standard_only.sh` script.


### StrainGE
1. Run the `straingst_build_db.sh` script to build the reference databases (one for each species)
2. Run the `strainge_inference_parallelized.sh` script to infer abundances 


## 4. Evaluation and plotting
All metrics calculation and visualization are contained in the `Strainify_paper/plotting_scripts/Figure_4.py` script. The same metrics used in the ChronoStrain paper are applied here, with one extra metric (`corr^2`) reported. 






