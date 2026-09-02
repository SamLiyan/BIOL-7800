# Differential expression: 2 year old vs 7 year old fish (testis gene expression study)

setwd("C:/Users/samit/Desktop/class_assignment_2")

#############STEP 1####################################
# edgeR quasi-likelihood (QL) test

# Load edgeR
library(edgeR)

# Read count data
count_data <- read.csv("Counts.csv", header = TRUE, row.names = 1)

# Annotate the gene Ids with their common gene names
library(rtracklayer)

annotation <- import("bluecatfish_annotation.gtf")
annotation_df <- as.data.frame(annotation)

# Build a gene_id which is a geneName lookup table from the annotation
gene_lookup <- unique(annotation_df[!is.na(annotation_df$transcript_id),
                                    c("transcript_id", "gene")])
names(gene_lookup) <- c("gene_id", "GeneName")

# Define the treatments / columns 1-3 = Age2, columns 4-6 = Age7
Condition <- factor(rep(c("Age2", "Age7"), each = 3), levels = c("Age2", "Age7"))

# Build DGEList
dge <- DGEList(counts = count_data, group = Condition)

# Filter low count genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalization / normalize the gene counts
dge <- calcNormFactors(dge, method = "TMM")

# Design matrix
design <- model.matrix(~0 + Condition)
colnames(design) <- levels(Condition)

# Estimate dispersion and fit QL model
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# Comparison - Age7 vs Age2
contrast <- makeContrasts(Age7 - Age2, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)


# Get results, filter by FDR < 0.01
results <- topTags(qlf, n = Inf)$table
results$gene_id <- rownames(results)
results <- merge(gene_lookup, results, by = "gene_id", all.y = TRUE)

sig <- results[results$FDR < 0.01, ]

# Split into up/down regulated (>= 2-fold change)
fold_change_cutoff <- 3
up <- sig[sig$logFC > log2(fold_change_cutoff), ]
down <- sig[sig$logFC < -log2(fold_change_cutoff), ]

# to get the gene count
gene_counts <- data.frame(Category = c("Upregulated", "Downregulated"),
                          Count = c(nrow(up), nrow(down)))
print(gene_counts)

# Save results
write.csv(results, "AllResults_Age7_vs_Age2.csv")
write.csv(up, "Upregulated_Age7_vs_Age2.csv")
write.csv(down, "Downregulated_Age7_vs_Age2.csv")

####################################STEP 2############################################

# Creating the PCA plot

library(ggplot2)
library(pheatmap)

logCPM <- cpm(dge, log = TRUE)
pca <- prcomp(t(logCPM), scale. = TRUE)
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Condition = Condition, Sample = colnames(logCPM))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("Age2" = "orangered", "Age7" = "navy")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

ggsave(filename = "PCA_Plot_Age2_vs_Age7.png", plot = p_pca, width = 6, height = 5, dpi = 600)

#Sample clustering heat map generation
library(pheatmap)

sample_cor <- cor(logCPM)
annotation_col <- data.frame(Condition = Condition, row.names = colnames(logCPM))

pheatmap(sample_cor,annotation_col = annotation_col,filename = "Heatmap_SampleCorrelation_Age2_vs_Age7.png")

########################################STEP 3###################################
#DEG Visualization

# Creating the volcano plot
library(ggplot2)

# Get results and label genes
results <- topTags(qlf, n = Inf)$table
results$Category <- "Not Significant"
results$Category[results$FDR < 0.01 & results$logFC > log2(fold_change_cutoff)] <- "Upregulated"
results$Category[results$FDR < 0.01 & results$logFC < -log2(fold_change_cutoff)] <- "Downregulated"

up_count <- sum(results$Category == "Upregulated")
down_count <- sum(results$Category == "Downregulated")

# Volcano plot
ggplot(results, aes(logFC, -log10(FDR), color = Category)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c(Upregulated = "orangered", Downregulated = "navy", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-log2(fold_change_cutoff), log2(fold_change_cutoff)), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  annotate("text", x = Inf, y = Inf, label = paste("Up:", up_count), hjust = 1.1, vjust = 2, color = "orangered", fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, label = paste("Down:", down_count), hjust = -0.1, vjust = 2, color = "navy", fontface = "bold") +
  theme_minimal()

ggsave("Volcano_Age7_vs_Age2.png", width = 6, height = 6, dpi = 300)


# Significant DEG heat map
library(dplyr)

# Get the top 20 most significant genes by FDR, excluding unnamed "LOC..." genes
top20_genes <- sig %>%
  filter(!startsWith(GeneName, "LOC")) %>%
  arrange(FDR) %>%
  pull(gene_id) %>%
  head(20)

top20_expr <- logCPM[rownames(logCPM) %in% top20_genes, ]

# Use gene names as row labels where available, gene_id otherwise
gene_labels <- gene_lookup$GeneName[match(rownames(top20_expr), gene_lookup$gene_id)]
gene_labels[is.na(gene_labels)] <- rownames(top20_expr)[is.na(gene_labels)]
rownames(top20_expr) <- gene_labels

pheatmap(top20_expr,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         filename = "Top20_SigGenes_Heatmap_Age7_vs_Age2.png")

###################################STEP 4#######################################

# Gene ontology analysis
library(gprofiler2)

# Use gene symbols  against zebrafish annotation,
# Drop LOC-prefixed / missing names first
up_names <- up$GeneName[!is.na(up$GeneName) & !startsWith(up$GeneName, "LOC")]
down_names <- down$GeneName[!is.na(down$GeneName) & !startsWith(down$GeneName, "LOC")]

# GO enrichment
go_up <- gost(query = up_names, organism = "drerio", sources = "GO")
go_down <- gost(query = down_names, organism = "drerio", sources = "GO")

go_up$result
go_down$result

# Save results
write.csv(go_up$result[, 1:13], "GO_Upregulated_Age7_vs_Age2.csv", row.names = FALSE)
write.csv(go_down$result[, 1:13], "GO_Downregulated_Age7_vs_Age2.csv", row.names = FALSE)


#Draw the GO bar plot for upregulated GO term at age 7
library(ggplot2)
library(dplyr)

go_colors <- c("GO:BP" = "steelblue", "GO:CC" = "darkorange", "GO:MF" = "forestgreen")

# Top 10 GO terms by significance - Upregulated
top_go_up <- go_up$result %>% arrange(p_value) %>% head(10)

ggplot(top_go_up, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value), fill = source)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = go_colors) +
  labs(x = NULL, y = "-log10(p-value)", fill = "GO Category", title = "Top GO Terms - Upregulated") +
  theme_minimal()

ggsave("GO_Upregulated_BarPlot.png", width = 8, height = 6, dpi = 300)

#Draw the GO bar plot for downregulated GO term at age 7
go_colors <- c("GO:BP" = "steelblue", "GO:CC" = "darkorange", "GO:MF" = "forestgreen")
# Top 10 GO terms by significance - Downregulated
top_go_down <- go_down$result %>% arrange(p_value) %>% head(10)
ggplot(top_go_down, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value), fill = source)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = go_colors) +
  labs(x = NULL, y = "-log10(p-value)", fill = "GO Category", title = "Top GO Terms - Downregulated") +
  theme_minimal()

ggsave("GO_Downregulated_BarPlot.png", width = 8, height = 6, dpi = 300)
