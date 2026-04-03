# 1) Download the SRA file
prefetch SRR13355226

# 2) Convert to FASTQ
fasterq-dump SRR13355226 --split-files -O .

# 3) Move and rename the FASTQ files
mkdir -p fastq
mv SRR13355226_1.fastq fastq/SRR13355226_r1.fq
mv SRR13355226_2.fastq fastq/SRR13355226_r2.fq