# Please change these paths to your local setup.

# Note that Strainify was run with the original reference genomes as it runs Parsnp to align the core genomes. To operate under the same strain definition, the same core genomes are then provided to StrainGE and ChronoStrain. 

./strainify \
    --genome_folder /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/gold_standard_genomes/ecoli \
    --fastq_folder /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/reads/paired_fastqs \
    --outdir /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/strainify_cami/force_include_all/ecoli \
    --weight_by_entropy \
    --parsnp_flags " -c --skip-ani-filter"

./strainify \
    --genome_folder /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/gold_standard_genomes/efaecium \
    --fastq_folder /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/reads/paired_fastqs \
    --outdir /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/strainify_cami/force_include_all/efaecium \
    --weight_by_entropy \
    --parsnp_flags " -c --skip-ani-filter"

./strainify \
    --genome_folder /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/gold_standard_genomes/kpneumoniae \
    --fastq_folder /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/reads/paired_fastqs \
    --outdir /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/strainify_cami/force_include_all/kpneumoniae \
    --weight_by_entropy \
    --parsnp_flags " -c --skip-ani-filter"

./strainify \
    --genome_folder /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/gold_standard_genomes/saureus \
    --fastq_folder /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/reads/paired_fastqs \
    --outdir /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/strainify_cami/force_include_all/saureus \
    --weight_by_entropy \
    --parsnp_flags " -c --skip-ani-filter"

./strainify \
    --genome_folder /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/gold_standard_genomes/spneumoniae \
    --fastq_folder /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/reads/paired_fastqs \
    --outdir /dodo/rl152/ChronoStrain/chronostrain_cami/cami_dataset/strainify_cami/force_include_all/spneumoniae \
    --weight_by_entropy \
    --parsnp_flags " -c --skip-ani-filter"