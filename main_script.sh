#!/usr/bin/env bash
set -euo pipefail

SRR=SRR2589044
REF_ACC=NC_000913.3

mkdir -p input_data work output_data

# 1. Reference genome (for comparison)
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${REF_ACC}&rettype=fasta&retmode=text" > input_data/ref.fa

# Download all paired-end reads for the run, split into forward/reverse fastq files
fastq-dump --split-files -O input_data ${SRR}

# 3. Assembly (de novo)
skesa --fastq input_data/${SRR}_1.fastq,input_data/${SRR}_2.fastq --cores 1 --memory 8 > output_data/assembly.fasta

# 4. Assembly stats — side-by-side comparison (assembly vs reference)
seqkit stats -a output_data/assembly.fasta input_data/ref.fa > output_data/comparison_stats.tsv

# 5. Annotation
prokka --outdir work/prokka --prefix genome --cpus 1 --force output_data/assembly.fasta
cp work/prokka/genome.txt output_data/prokka_summary.txt
cp work/prokka/genome.gff output_data/prokka_annotation.gff
cp work/prokka/genome.faa output_data/prokka_proteins.faa

# 6. Checksums — data file(s) used + every output file, explicitly listed
sha256sum \
  input_data/ref.fa \
  input_data/${SRR}_1.fastq \
  input_data/${SRR}_2.fastq \
  output_data/assembly.fasta \
  output_data/comparison_stats.tsv \
  output_data/prokka_summary.txt \
  output_data/prokka_annotation.gff \
  output_data/prokka_proteins.faa \
  > checksum_values.txt

# 7. Clean up Prokka's raw intermediate/duplicate files
# (only genome.txt / genome.gff / genome.faa were needed, and those are
# already copied into output_data/ above)
rm -f \
  work/prokka/genome.err \
  work/prokka/genome.ffn \
  work/prokka/genome.fna \
  work/prokka/genome.fsa \
  work/prokka/genome.gbk \
  work/prokka/genome.log \
  work/prokka/genome.sqn \
  work/prokka/genome.tbl