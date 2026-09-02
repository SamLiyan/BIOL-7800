# Differential Gene Expression Analysis - Blue Catfish Testis Tissue Samples

# Reproducibility and Accuracy Test of the Pipeline

## 1. Known answer

Two methods use under this category.

**Method 1: housekeeping genes.** Known, from outside this pipeline. These genes should stay stable with age. Should not show up as significant.

```r
housekeeping_genes <- c("actb", "gapdh", "tubb", "rpl13a", "b2m")
hk_check <- results[results$GeneName %in% housekeeping_genes,
                    c("GeneName", "gene_id", "logFC", "FDR")]
sum(hk_check$FDR < 0.01, na.rm = TRUE)  # expect 0
```

Result: 0 of 5 found.

**Method 2: reverse the contrast.** logFC must flip sign, FDR must stay same.

```r
contrast_reversed <- makeContrasts(Age2 - Age7, levels = design)
qlf_reversed <- glmQLFTest(fit, contrast = contrast_reversed)
results_reversed <- topTags(qlf_reversed, n = Inf, sort.by = "none")$table
results_original <- topTags(qlf, n = Inf, sort.by = "none")$table

isTRUE(all.equal(results_original$logFC, -results_reversed$logFC))
isTRUE(all.equal(results_original$FDR, results_reversed$FDR))
```

Result: logFC signs reversed - TRUE. FDR unchanged - TRUE. Passed.

## 2. Invariant

A gene can't be both up and down at once. Must hold no matter what.

```r
overlap_genes <- intersect(up$gene_id, down$gene_id)
stopifnot("Up and down regulated sets must not overlap" = length(overlap_genes) == 0)
```

Result: 2146 up, 2013 down. Overlap = 0. Holds. Passed.

## 3. Positive control

Planted a 3-fold change in one gene that was not significant. Checked if pipeline finds it back.

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

Result: gene XM_053610102.1. Planted logFC = 1.58, observed = 2.17. FDR = 5.0e-05. Detected. Passed.

## 4. Order of magnitude

Expected range written before running: 15-35% of genes tested. Reason: Age2 vs Age7 is immature vs mature testis, not small aging shift. Cell type makeup changes a lot at this stage.

```r
real_deg_count <- nrow(up) + nrow(down)
pct_deg <- 100 * real_deg_count / nrow(dge)
```

Result: 16188 genes tested. 4159 DEGs (2146 up, 2013 down). 25.7% of genes. Falls in expected range. Passed.

## 5. Determinism

Ran full analysis twice, from a fresh R session each time. Compared output file checksums.

```r
md5_run1 <- tools::md5sum("AllResults_Age7_vs_Age2.csv")
# restart R, re-run full script, then:
md5_run2 <- tools::md5sum("AllResults_Age7_vs_Age2.csv")
identical(md5_run1, md5_run2)
```

Result: both runs gave checksum `4082f8c881e696eecd3b45f77949c202`. Identical. Passed.
