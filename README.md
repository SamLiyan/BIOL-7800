# Bacterial Genome Assembly and Annotation (SKESA + Prokka)

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

All versions are pinned in `environment_info.yml`. `curl` is left
unpinned deliberately — it's only used to fetch the reference FASTA over
HTTP and has no effect on reproducibility, and pinning it to a newer
build can pull in system libraries that conflict with SKESA's older
dependency chain during environment solving.

The script is set to run with `--cores 8`, `--cpus 8`, and 10 GB of RAM,
matching a machine with 10 CPUs and 11 GB available (leaving a small
margin for the OS). Before running, check your own machine's resources
with `nproc` and `free -h`, and adjust the `--cores`, `--cpus`, and
`--memory` values in `main_script.sh` to match what you actually have
available. Using more or fewer cores will change runtime, and may also
change the exact ordering of contigs SKESA produces (see Reproducibility
Notes below).

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
3. Assembles the reads de novo into contigs using SKESA (`--cores 8`, `--memory 10`).
4. Generates a side-by-side comparison of assembly vs. reference
   statistics using SeqKit.
5. Annotates the assembly using Prokka (`--cpus 8`, raw output written
   to `work/prokka/`), then copies the three files needed downstream
   (`genome.txt`, `genome.gff`, `genome.faa`) into `output_data/`.
6. Computes SHA256 checksums of every downloaded input and generated
   output file, written to `checksum_values.txt`.
7. Deletes Prokka's raw intermediate/duplicate files from `work/prokka/`,
   keeping only `genome.txt`, `genome.gff`, `genome.faa`, and
   `genome.tsv`.

On the machine used for this submission (`--cores 8`, `--cpus 8`), the
full run took approximately **[UPDATE WITH YOUR MEASURED TIME]**.
Runtime will vary depending on how many cores/RAM you allocate.

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
- Assembly: 82 contigs, 4,532,199 bp total length (reference genome is
  4,641,652 bp — about 97.6% recovered)
- Annotation: 4,188 predicted coding sequences (CDS), 77 tRNAs, and 3 rRNAs.

## Pipeline Details
SKESA builds and extends the assembly graph across a series of
progressively longer k-mer values, then resolves repeat
regions in a separate mate-pair-connection phase that iterates through
k-mer lengths again. This makes it fast and largely deterministic
compared to assemblers that rely on randomized graph traversal.

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
The `output_data/` entries should report `OK`. The `input_data/` entries
(`ref.fa`, `SRR2589044_1.fastq`, `SRR2589044_2.fastq`) will report
"No such file or directory" instead, since step 8 deletes that folder
after the checksums are recorded — those three lines exist for
record-keeping of exactly what was downloaded, not for later
verification.

## Folder Layout
```
.
├── README.md
├── environment_info.yml
├── main_script.sh
├── checksum_values.txt
├── work/
│   └── prokka/
│       ├── genome.txt
│       ├── genome.gff
│       ├── genome.faa
│       └── genome.tsv
└── output_data/
    ├── assembly.fasta
    ├── comparison_stats.tsv
    ├── prokka_summary.txt
    ├── prokka_annotation.gff
    └── prokka_proteins.faa
```
