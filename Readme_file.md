# Bacterial Genome Assembly and Annotation

## Background
Whole-genome sequencing produces raw reads that need to be assembled into
contigs and then annotated before they are biologically useful. This
assignment builds a small, reproducible pipeline that takes public
paired-end Illumina reads from *E. coli*, assembles them de novo with
SKESA, compares the assembly against the reference genome, and annotates
the assembly with Prokka.

## Aim
To assemble a bacterial genome from raw public sequencing reads and
annotate its gene content, then check how the assembly compares to the
known reference genome in terms of size and contiguity.

## Datasets Used
- Reference genome: *Escherichia coli* K-12 MG1655, NCBI accession
  **NC_000913.3**. Downloaded directly from NCBI E-utilities and used
  only for the stats comparison, not for the assembly itself.
- Sequencing reads: SRA run **SRR2589044** (paired-end Illumina WGS
  reads, *E. coli*). The full run is used (1,107,090 read pairs,
  approximately 70x coverage).

No data files are included in this repository. `main_script.sh`
downloads everything it needs directly from NCBI/SRA when it runs.

## Requirements
- SRA Toolkit (fastq-dump) 3.1.1
- SKESA 2.4.0
- SeqKit 2.8.2
- Prokka 1.14.6
- curl

All versions are pinned in `environment_info.yml`.

The script is set to run single-threaded (`--cores 1`, `--cpus 1`) with
8 GB of RAM as a safe default. Before running, check your own machine's
resources with `nproc` and `free -h`, and adjust the `--cores`,
`--cpus`, and `--memory` values in `run_analysis.sh` if you have more
available. Using more cores will speed up both SKESA and Prokka, but may
change the exact ordering of contigs SKESA produces (see
Reproducibility Notes below).

## Setting Up the Environment
Create the Conda environment:
```
conda env create -f environment_info.yml
```

Activate it:
```
conda activate assembly-annotation
```

## Execution Steps
Run the full pipeline with:
```
bash main_script.sh
```

The script performs the following steps in order:
1. Downloads the *E. coli* reference genome (NC_000913.3) into `input_data/`.
2. Downloads the full SRR2589044 read set from SRA into `input_data/` using `fastq-dump`.
3. Assembles the reads de novo into contigs using SKESA.
4. Generates a side-by-side comparison of assembly vs. reference
   statistics using SeqKit.
5. Annotates the assembly using Prokka (raw output written to `work/prokka/`).
6. Computes SHA256 checksums of every downloaded input and generated
   output file, written to `checksum_values.txt`.
7. Deletes Prokka's raw intermediate/duplicate files from `work/prokka/`,
   keeping only the three annotation files already copied into
   `output_data/`.

On the machine used for this submission (single-threaded, `--cores 1` /
`--cpus 1`), the full run took approximately **28 minutes**: SKESA
assembly took about 21 minutes across its 20 progressive k-mer rounds,
and Prokka annotation took 7.10 minutes (as reported directly in
Prokka's own log). Runtime will be shorter on a machine with more cores
allocated.

## Generated Files
```
output_data/
├── assembly.fasta            # de novo genome assembly from SKESA
├── comparison_stats.tsv      # assembly vs. reference stats, side by side
├── prokka_summary.txt        # Prokka gene-count summary
├── prokka_annotation.gff     # full gene annotation
└── prokka_proteins.faa       # predicted protein sequences
```

Results from the completed run:
- Assembly: assembly recovered 4,532,199 bp out of the reference's 4,641,652 bp — about 97.6% of the genome's total length. 82 contigs
- Annotation: 4,188 predicted coding sequences (CDS), 78 tRNAs, 3 rRNAs,
  and 2 CRISPR arrays.

## Pipeline Details
SKESA builds and extends the assembly graph across a series of
progressively longer k-mer values (19 up to 511), then resolves repeat
regions in a separate mate-pair-connection phase that iterates through
k-mer lengths again. This makes it fast and deterministic compared to
assemblers that rely on randomized graph traversal.

Prokka annotates the assembly by first predicting tRNAs (Aragorn), rRNAs
(Barrnap), and CRISPR arrays (Minced), then predicting coding sequences
with Prodigal. Each predicted CDS is searched in turn against curated
databases (IS elements, then AMR genes, then Swiss-Prot) and, for
anything still unmatched, against the HAMAP database using HMMER.
Anything left over is labeled "hypothetical protein."

## Checking Results
After running the analysis, verify file integrity with:
```
sha256sum -c checksum_values.txt
```
Every listed file should report `OK`.

## Folder Layout
```
.
├── Readme_file.md
├── environment_info.yml
├── main_script.sh
├── checksum_values.txt
└── output_data/
    ├── assembly.fasta
    ├── comparison_stats.tsv
    ├── prokka_summary.txt
    ├── prokka_annotation.gff
    └── prokka_proteins.faa
```
