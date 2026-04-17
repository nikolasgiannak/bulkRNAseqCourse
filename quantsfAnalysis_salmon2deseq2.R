# =============================================================================
# DESeq2 from Salmon output — quant.sf (transcript-level)
# Gene symbol annotation via org.Mm.eg.db (no internet required)
# =============================================================================

# --- 1. Install & load packages ----------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("tximport", "DESeq2", "GenomicFeatures", "apeglm",
                       "clusterProfiler", "org.Mm.eg.db", "EnhancedVolcano"),
                     update = FALSE, ask = FALSE)
install.packages(c("ggplot2", "pheatmap"))

library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(apeglm)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Mm.eg.db)
library(EnhancedVolcano)

# --- 2. Sample table ---------------------------------------------------------

sampleinfo <- data.frame(
  sample    = c("WT_1", "WT_2", "WT_3", "TKO_1", "TKO_2", "TKO_3"),
  condition = factor(c("WT", "WT", "WT", "KO", "KO", "KO")),
  row.names = c("WT_1", "WT_2", "WT_3", "TKO_1", "TKO_2", "TKO_3")
)

# --- 3. Build tx2gene from GTF -----------------------------------------------

txdb <- makeTxDbFromGFF(
  "/Users/nikolasgiannak/Desktop/Teaching/RNAseqCourse/Data/Mus_musculus.GRCm38.correct_chrom_names.102.chr.gtf.gz",
  format = "gtf"
)
k       <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")

# --- 4. Point to quant.sf files ----------------------------------------------

files <- file.path(
  "/Users/nikolasgiannak/Desktop/Teaching/RNAseqCourse/Data",
  sampleinfo$sample,
  "quant.sf"
)
names(files) <- sampleinfo$sample
stopifnot(all(file.exists(files)))

# --- 5. Import with tximport -------------------------------------------------

txi <- tximport(
  files,
  type            = "salmon",
  tx2gene         = tx2gene,
  ignoreTxVersion = TRUE
)

# --- 6. Create DESeq2 object -------------------------------------------------

dds <- DESeqDataSetFromTximport(txi, colData = sampleinfo, design = ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]

# --- 7. Run DESeq2 -----------------------------------------------------------

dds <- DESeq(dds)
resultsNames(dds)

# --- 8. Extract & shrink results ---------------------------------------------

res        <- results(dds, contrast = c("condition", "KO", "WT"))
res_shrunk <- lfcShrink(dds, coef = "condition_WT_vs_KO", type = "apeglm")
summary(res_shrunk, alpha = 0.05)

# --- 9. Add gene symbols via org.Mm.eg.db ------------------------------------

ensembl_ids <- sub("\\..*", "", rownames(res_shrunk))

gene_map <- select(
  org.Mm.eg.db,
  keys    = ensembl_ids,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL", "GENENAME")
)

# Remove duplicate mappings — keep first occurrence
gene_map <- gene_map[!duplicated(gene_map$ENSEMBL), ]

# Merge into results
res_df         <- as.data.frame(res_shrunk)
res_df$ENSEMBL <- sub("\\..*", "", rownames(res_df))

res_df <- merge(res_df, gene_map, by = "ENSEMBL", all.x = TRUE)
res_df <- res_df[order(res_df$padj), ]

write.csv(res_df, "WT_vs_KO_DESeq2_quantsf_results.csv", row.names = FALSE)

# --- 10. Filter significant genes --------------------------------------------

sig <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
cat("Upregulated:  ", sum(sig$log2FoldChange >  1, na.rm = TRUE), "\n")
cat("Downregulated:", sum(sig$log2FoldChange < -1, na.rm = TRUE), "\n")

# --- 11. QC plots ------------------------------------------------------------

vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "condition")
plotMA(res_shrunk, ylim = c(-5, 5))
plotDispEsts(dds)

# --- 12. EnhancedVolcano -----------------------------------------------------

# Remove NA symbols, deduplicate keeping lowest padj per symbol
res_plot <- res_df[!is.na(res_df$SYMBOL), ]
res_plot <- res_plot[order(res_plot$padj), ]
res_plot <- res_plot[!duplicated(res_plot$SYMBOL), ]
rownames(res_plot) <- res_plot$SYMBOL

EnhancedVolcano(res_plot,
    lab             = rownames(res_plot),
    x               = 'log2FoldChange',
    y               = 'pvalue',
    selectLab       = c('Dnmt1', 'Scml2', 'Dnmt3b',
                        'Utrn', 'Rbm44', 'Rrm2b'),
    xlab            = bquote(~Log[2]~ 'fold change'),
    pCutoff      = 10e-14,
    FCcutoff     = 2.0,
    pointSize    = 3.0,
    labSize      = 6.0,
    labCol       = 'black',
    labFace      = 'bold',
    boxedLabels  = TRUE,
    parseLabels = TRUE,
    col = c('gray', 'lightblue3', 'maroon4', 'brown3'),
    colAlpha     = 5/5,
    legendPosition  = 'right',
    legendLabSize   = 14,
    legendIconSize  = 4.0,
    drawConnectors  = TRUE,
    widthConnectors = 0.5,
    colConnectors   = 'black')

# --- 13. Heatmap (ComplexHeatmap) --------------------------------------------
# --- Heatmap with gene symbols only ------------------------------------------

my_col_fun <- colorRamp2(c(-2, 0, 2), c("gray", "white", "coral3"))

up   <- res_df$ENSEMBL[res_df$log2FoldChange >  1 & res_df$padj < 0.05 & !is.na(res_df$padj)]
down <- res_df$ENSEMBL[res_df$log2FoldChange < -1 & res_df$padj < 0.05 & !is.na(res_df$padj)]

# Match back to rownames in dds
dds_ids  <- sub("\\..*", "", rownames(dds))
up_idx   <- rownames(dds)[dds_ids %in% up]
down_idx <- rownames(dds)[dds_ids %in% down]
top_genes <- c(up_idx, down_idx)

mat        <- assay(vst(dds))[top_genes, ]
mat_scaled <- t(scale(t(mat)))

# Get gene symbols for each row
clean_ids   <- sub("\\..*", "", rownames(mat_scaled))
row_symbols <- gene_map$SYMBOL[match(clean_ids, gene_map$ENSEMBL)]

# *** Keep only rows that have a gene symbol ***
has_symbol  <- !is.na(row_symbols) & row_symbols != ""
mat_scaled  <- mat_scaled[has_symbol, ]
row_symbols <- row_symbols[has_symbol]

# Update split vector to match filtered matrix
up_clean   <- sub("\\..*", "", up_idx)
down_clean <- sub("\\..*", "", down_idx)
row_clean  <- sub("\\..*", "", rownames(mat_scaled))

split <- factor(
  ifelse(row_clean %in% up_clean, "Upregulated", "Downregulated"),
  levels = c("Upregulated", "Downregulated")
)

# Set gene symbols as rownames
rownames(mat_scaled) <- row_symbols

col_ann <- HeatmapAnnotation(
  condition = sampleinfo$condition,
  col = list(condition = c("WT" = "firebrick2", "KO" = "cyan4"))
)

Heatmap(
  mat_scaled,
  top_annotation    = col_ann,
  name              = "Z-score",
  col               = my_col_fun,
  row_split         = split,
  row_title_gp      = gpar(fontsize = 12, fontface = "bold"),
  show_row_names    = TRUE,    # now safe to show — all rows have symbols
  show_column_names = TRUE,
  row_names_gp      = gpar(fontsize = 7)   # reduce size if many genes
)
# --- 14. GO enrichment -------------------------------------------------------

sig_symbols <- sig$SYMBOL[!is.na(sig$SYMBOL) & sig$SYMBOL != ""]

ego <- enrichGO(
  gene          = sig_symbols,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)
barplot(ego, showCategory = 20)
dotplot(ego,  showCategory = 20)

# --- 15. GSEA ----------------------------------------------------------------

ranked        <- res_df$log2FoldChange
names(ranked) <- res_df$SYMBOL
ranked        <- ranked[!is.na(names(ranked)) & names(ranked) != ""]
ranked        <- ranked[!duplicated(names(ranked))]
ranked        <- sort(ranked, decreasing = TRUE)

gsea_res <- gseGO(
  geneList = ranked,
  OrgDb    = org.Mm.eg.db,
  keyType  = "SYMBOL",
  ont      = "BP"
)
dotplot(gsea_res, showCategory = 20)
