# Strainify_paper

## Simulated datasets
[ART Illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) is required to create the simulated datasets
Create simulated datasets using the `simulate_reads.sh` script in each subfolder of the `/simulated_exp` folder. Please set the `art_bin` configuration in each script to the path of your ART Illumina binary. 

Note: the simulated datasets used in the manuscript will be available on Zenodo soon. 

## Running Strainify on simulated and mock datasets
Make sure the `/Strainify_paper` directory and the `/Strainify` directory are located in the same parent directory. Run Strainify from the `/Strainify` directory using the following command:
```bash
../Strainify_paper/simulated_exp/strainify/run_strainify.sh
../Strainify_paper/4_strain_mock_community/strainify/run_strainify.sh
```
Note: these scripts assume that you installed Strainify via git. For runtime and memory benchmarking, please uncomment the lines related to the `benchmark` feature in the Snakefile (lines 47, 78, 88, 104, 117, 133, 143, 154, 166, 176, 190, 209, 224, 243, 262, 283). The runtime and memory benchmarking results analysis scripts are in the `strainify_timings` directory.

## Running StrainScan, StrainR2 and PHLAME on benchmarking datasets
Please find each tool's own folder in this repository and run the scripts inside. Please change the paths in all scripts to match your local setup. 


