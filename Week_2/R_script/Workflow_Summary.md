# Workflow Summary - Blue Catfish Testis DGE Analysis

Compared testis gene expression between 2 year old (Age2) and 7 year old (Age7) blue catfish, n=3 per group.

Full code: `R_markdown_main RNA-seq analysis` (main analysis, steps below), `Reproducibility and Accuracy Test of the Pipeline.html` (checks further down).

### Step 1: Differential expression (edgeR)

Read raw gene counts (`Counts.csv`). Annotated gene IDs with gene names from `bluecatfish_annotation.gtf`. Set Age2/Age7 groups. Built a DGEList, filtered low-count genes, TMM normalized. Fit an edgeR quasi-likelihood model. Tested the Age7 vs Age2 contrast. Called significant genes at FDR < 0.01, >= 3-fold change. Saved results, up, and down gene lists as CSVs.

### Step 2: PCA and sample clustering

Computed log-CPM values. Ran PCA, plotted samples colored by age, labeled by sample name. Computed sample correlation matrix, plotted as a heatmap with age annotation, to check samples cluster by group.

### Step 3: DEG visualization

Built a volcano plot (logFC vs -log10 FDR), colored by up/down/not significant, with up/down counts on the plot. Pulled the top 20 significant genes by FDR (dropping unnamed "LOC" genes), plotted their expression as a heatmap.

### Step 4: GO enrichment

Took gene symbols for up and down gene sets (dropping unnamed genes). Ran GO enrichment against zebrafish annotation (g:Profiler), since gene symbols were assigned by homology to zebrafish. Saved GO results, plotted top 10 terms per direction as horizontal bar plots, colored by GO category (BP/CC/MF).

## Testing Mindset Controls

### 1. Known answer

Method 1: housekeeping genes should not be significant.

```r
housekeeping_genes <- c("actb", "gapdh", "tubb", "rpl13a", "b2m")
hk_check <- results[results$GeneName %in% housekeeping_genes,
                    c("GeneName", "gene_id", "logFC", "FDR")]
sum(hk_check$FDR < 0.01, na.rm = TRUE)  # expect 0
```

Result: pending re-run (earlier attempt had an indexing bug).

Method 2: reversing the contrast must flip logFC sign, keep FDR same.

```r
contrast_reversed <- makeContrasts(Age2 - Age7, levels = design)
qlf_reversed <- glmQLFTest(fit, contrast = contrast_reversed)
results_reversed <- topTags(qlf_reversed, n = Inf, sort.by = "none")$table
results_original <- topTags(qlf, n = Inf, sort.by = "none")$table

isTRUE(all.equal(results_original$logFC, -results_reversed$logFC))
isTRUE(all.equal(results_original$FDR, results_reversed$FDR))
```

Result: logFC reversed - TRUE. FDR unchanged - TRUE. Passed.

### 2. Invariant

Up and down gene sets can't overlap.

```r
overlap_genes <- intersect(up$gene_id, down$gene_id)
stopifnot(length(overlap_genes) == 0)
```

Result: 2146 up, 2013 down, overlap = 0. Passed.

### 3. Positive control

Planted a 3-fold change in one non-significant gene, checked it's detected.

```r
spike_gene <- setdiff(rownames(dge), c(up$gene_id, down$gene_id))[1]
age7_cols <- colnames(count_data)[4:6]
spiked <- count_data
spiked[spike_gene, age7_cols] <- round(spiked[spike_gene, age7_cols] * 3)

s_dge <- DGEList(counts = spiked, group = Condition)
s_dge <- s_dge[filterByExpr(s_dge), , keep.lib.sizes = FALSE]
s_dge <- calcNormFactors(s_dge)
s_dge <- estimateDisp(s_dge, design)
s_fit <- glmQLFit(s_dge, design)
s_qlf <- glmQLFTest(s_fit, contrast = contrast)
s_res <- topTags(s_qlf, n = Inf)$table
```

Result: gene XM_053610102.1. Planted logFC 1.58, observed 2.17. FDR 5.0e-05. Detected. Passed.

### 4. Order of magnitude

Expected range set before running: 15-35% of genes tested. Reason: Age2 vs Age7 spans juvenile to mature testis, a big developmental shift, not a subtle aging change.

```r
real_deg_count <- nrow(up) + nrow(down)
pct_deg <- 100 * real_deg_count / nrow(dge)
```

Result: 16188 genes tested, 4159 DEGs (25.7%). Within expected range. Passed.

### 5. Determinism

Ran full analysis twice from a fresh R session, compared output checksums.

```r
md5_run1 <- tools::md5sum("AllResults_Age7_vs_Age2.csv")
# restart R, re-run full script, then:
md5_run2 <- tools::md5sum("AllResults_Age7_vs_Age2.csv")
identical(md5_run1, md5_run2)
```

Result: both runs gave checksum `4082f8c881e696eecd3b45f77949c202`. Identical. Passed.

## Environment and Software Versions

| Package | Version |
|---|---|
| edgeR | 4.8.2 |
| rtracklayer | 1.70.1 |
| ggplot2 | 4.0.3 |
| pheatmap | 1.0.13 |
| dplyr | 1.2.1 |
| gprofiler2 | 0.2.4 |

## Data Provenance

| File | MD5 |
|---|---|---|
| Counts.csv | f0545873c72926d13a89ab7f2ec07159 |
| bluecatfish_annotation.gtf | cdfaad498c02494354b8698cb0669d48 |

## GitHub
GitHub repo: https://github.com/SamLiyan/BIOL-7800-Computational-Biology-Colloquium/tree/main/Assignment_2
