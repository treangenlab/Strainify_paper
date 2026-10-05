# Strainify_paper

## Simulated datasets
[ART Illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) is required to create the simulated datasets
Create simulated datasets using the `simulate_reads.sh` script in each subfolder of the `/simulated_exp` folder. Please set the `art_bin` configuration in each script to the path of your ART Illumina binary. 

The simulated datasets used in the manuscript are also available on Zenodo (10.5281/zenodo.23107250)

## Running Strainify on simulated and mock datasets
Make sure the `/Strainify_paper` directory and the `/Strainify` directory are located in the same parent directory. Run Strainify from the `/Strainify` directory using the following command:
```bash
../Strainify_paper/simulated_exp/strainify/run_strainify.sh
../Strainify_paper/4_strain_mock_community/strainify/run_strainify.sh
```
Note: these scripts assume that you installed Strainify via git. 

## Running StrainScan, StrainR2, PHLAME, PanTax and StrainGE on benchmarking datasets
Please find each tool's own folder in this repository and run the scripts inside. Please change the paths in all scripts to match your local setup. 

## Running Strainify on <i>B. ovatus</i> dataset
All scripts required for this dataset are in the `B_ovatus` directory.
1. Download the fastq files by running the `retrieve_srr.sh` script
2. The reference genomes for this dataset are in the `all_fastas` folder (the script to download the reference genomes is `retrieve_fastas.sh` and it will download the assemblies in `PRJNA544527_AssemblyDetails.txt`)
3. Run Strainify via `run_strainify.sh`. Make sure to run it from the Strainify directory. Please modify the paths to match your local setup.
4. To analyze the data, run the following scripts: `reorder_srr.py`, `extract_genome_names.py` and `stack_plot.py`. Please change the paths in these scripts to match your local setup. You will need the `srr_list.txt` file for the `reorder_srr.py` script and the `assembly_names.tree` files for the `extract_genome_names.py` script. 

## CAMI II dataset
All relevant scripts are in the `CAMI` folder. Detailed instructions to run the analyses are here: [Strainify_paper/CAMI/CAMI.md](Strainify_paper/CAMI/CAMI.md).

## UMB dataset
All relevant scripts are in the `UMB` folder.
### Download dataset 
Run `download_dataset.sh` to download the UMB18 samples, and then run `copy_fastqs.sh` to copy the stool samples to another folder (only stool samples are included in our analyses). 

Run `get_genome_input_for_strainify.sh` to download reference genomes for this dataset.

### Run Strainify
Run the `run_strainify.sh` script. Make sure to update the paths based on your local setup. 

