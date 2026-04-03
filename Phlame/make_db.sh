#!/bin/sh

# Please change all paths in this script to match your local setup
mkdir -p makedb/db_fastq
mkdir -p makedb/db_sam
mkdir -p makedb/db_bam
mkdir -p makedb/counts
wgsim -e 0.0 -d 500 -N 1000000 -1 150 -2 150 -r 0.0 -R 0.0 -X 0.0 genomes/E24377A.fna makedb/db_fastq/E24377A_r1.fastq makedb/db_fastq/E24377A_r2.fastq
wgsim -e 0.0 -d 500 -N 1000000 -1 150 -2 150 -r 0.0 -R 0.0 -X 0.0 genomes/H10407.fna makedb/db_fastq/H10407_r1.fastq makedb/db_fastq/H10407_r2.fastq
wgsim -e 0.0 -d 500 -N 1000000 -1 150 -2 150 -r 0.0 -R 0.0 -X 0.0 genomes/Sakai.fna makedb/db_fastq/Sakai_r1.fastq makedb/db_fastq/Sakai_r2.fastq
wgsim -e 0.0 -d 500 -N 1000000 -1 150 -2 150 -r 0.0 -R 0.0 -X 0.0 genomes/UTI89.fna makedb/db_fastq/UTI89_r1.fastq makedb/db_fastq/UTI89_r2.fastq
gzip makedb/db_fastq/E24377A_r1.fastq makedb/db_fastq/E24377A_r2.fastq
gzip makedb/db_fastq/H10407_r1.fastq makedb/db_fastq/H10407_r2.fastq
gzip makedb/db_fastq/Sakai_r1.fastq makedb/db_fastq/Sakai_r2.fastq
gzip makedb/db_fastq/UTI89_r1.fastq makedb/db_fastq/UTI89_r2.fastq
bowtie2-build -q genomes/Sakai.fna reference_genome_idx/Sakai_idx

bowtie2 -X 2000 --no-mixed --dovetail -x reference_genome_idx/Sakai_idx -1 makedb/db_fastq/E24377A_r1.fastq.gz -2 makedb/db_fastq/E24377A_r2.fastq.gz -S makedb/db_sam/E24377A.sam
bowtie2 -X 2000 --no-mixed --dovetail -x reference_genome_idx/Sakai_idx -1 makedb/db_fastq/UTI89_r1.fastq.gz -2 makedb/db_fastq/UTI89_r2.fastq.gz -S makedb/db_sam/UTI89.sam
bowtie2 -X 2000 --no-mixed --dovetail -x reference_genome_idx/Sakai_idx -1 makedb/db_fastq/H10407_r1.fastq.gz -2 makedb/db_fastq/H10407_r2.fastq.gz -S makedb/db_sam/H10407.sam
bowtie2 -X 2000 --no-mixed --dovetail -x reference_genome_idx/Sakai_idx -1 makedb/db_fastq/Sakai_r1.fastq.gz -2 makedb/db_fastq/Sakai_r2.fastq.gz -S makedb/db_sam/Sakai.sam

samtools view -bS makedb/db_sam/E24377A.sam | samtools sort - -o makedb/db_bam/E24377A.bam
samtools view -bS makedb/db_sam/H10407.sam | samtools sort - -o makedb/db_bam/H10407.bam
samtools view -bS makedb/db_sam/UTI89.sam | samtools sort - -o makedb/db_bam/UTI89.bam
samtools view -bS makedb/db_sam/Sakai.sam | samtools sort - -o makedb/db_bam/Sakai.bam

samtools index makedb/db_bam/E24377A.bam
samtools index makedb/db_bam/H10407.bam
samtools index makedb/db_bam/UTI89.bam
samtools index makedb/db_bam/Sakai.bam


phlame counts -i makedb/db_bam/E24377A.bam -r genomes/Sakai.fna -o makedb/counts/E24377A.counts
phlame counts -i makedb/db_bam/H10407.bam -r genomes/Sakai.fna -o makedb/counts/H10407.counts
phlame counts -i makedb/db_bam/Sakai.bam -r genomes/Sakai.fna -o makedb/counts/Sakai.counts
phlame counts -i makedb/db_bam/UTI89.bam -r genomes/Sakai.fna -o makedb/counts/UTI89.counts

phlame cmt -i makedb/counts_files.txt -s makedb/sample_names.txt -r genomes/Sakai.fna -o makedb/ecoli.pickle.gz

phlame tree -i makedb/ecoli.pickle.gz -p makedb/ecoli.phylip -r makedb/ecoli_phylip2names.txt -o makedb/ecoli.tre --rescale

# Run the following command instead for low cov mode of phlame 
#phlame tree -i makedb/ecoli.pickle.gz -p makedb/ecoli_2.phylip -r makedb/ecoli_phylip2names_2.txt -o makedb/ecoli_2.tre --rescale --min_cov 1

phlame makedb -i makedb/ecoli.pickle.gz -t makedb/rescaled_ecoli.tre -o makedb/ecoli_db.classifier -p makedb/ecoli_cladeIDs.txt --min_branchlen 0 --min_leaves 1 --min_snps 1 -c makedb/ecoli_manual_clade_ID.txt

# Run the following command instead for low cov mode of phlame
#phlame makedb -i makedb/ecoli.pickle.gz -t makedb/rescaled_ecoli_2.tre -o makedb/ecoli_db_2.classifier -p makedb/ecoli_cladeIDs_2.txt --min_branchlen 0 --min_leaves 1 --min_snps 1 -c makedb/ecoli_manual_clade_ID.txt