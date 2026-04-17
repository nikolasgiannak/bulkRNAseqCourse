
#   DESeq2 from Salmon output — full walkthrough
# The standard approach uses tximport to load Salmon counts into R, 
# then DESeq2 for differential expression.



# 1.	Install required packages (R)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# This updates all Bioconductor packages to match your R version
BiocManager::install() 

#   Validate your installation:
#  If the error persists, use BiocManager's built-in validation tool to identify specific broken dependencies:
BiocManager::valid()

# install out of date packages
BiocManager::install(c(
  "arm", "devtools", "ggsci", "glue", "leidenAlg", "mvtnorm", "openssl", "pak", "parallelly", "Rcpp", "renv", "rJava", "S7", "scater", "V8"
), update = TRUE, ask = FALSE, force = TRUE)

BiocManager::install("S4Vectors", type = "source")



BiocManager::install(c("tximport", "DESeq2", "tximeta", "AnnotationDbi"), force =TRUE)
install.packages("readr")
BiocManager::install("GenomicFeatures", force =TRUE)
# 2.	Set up your sample table
library(tximport)
library(DESeq2)

# Define sample info
sampleinfo <- data.frame(
  sample = c("WT_1", "WT_2", "WT_3", "TKO_1", "TKO_2", "TKO_3"),
  condition = factor(c("WT", "WT", "WT", "KO", "KO", "KO")),
  row.names = c("WT1", "WT2", "WT2", "KO1", "KO2", "KO3")
)

# Point to quant.sf files (transcript-level)
files <- file.path("/Users/nikolasgiannak/Desktop/Teaching/RNAseqCourse/Data", sampleinfo$sample, "quant.sf")
names(files) <- sampleinfo$sample

# 3. Build tx2gene map from your GTF
# This maps transcript IDs → gene IDs, required by tximport.
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("/Users/nikolasgiannak/Desktop/Teaching/RNAseqCourse/Data/Mus_musculus.GRCm38.correct_chrom_names.102.chr.gtf.gz", format = "gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
# tx2gene has 2 cols: TXNAME, GENEID

#4. Import with tximport
#You have two options depending on which file you use:
#  Option A — use quant.sf (transcript-level → summarised to gene)

txi <- tximport(files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE  # strips version suffixes like ENST00000123.4 → ENST00000123)
  
#  Option B — use quant.genes.sf directly (already gene-level)
  
  gene_files <- file.path("/Users/nikolasgiannak/Desktop/Teaching/RNAseqCourse/Data", sampleinfo$sample, "quant.genes.sf")
  names(gene_files) <- sampleinfo$sample
  
  txi <- tximport(
    gene_files,
    type = "salmon",
    txOut = FALSE,       # already gene-level, no summarisation needed
    tx2gene = NULL)
  
#  Recommendation: Option A (quant.sf) is preferred — it uses the full transcript-level uncertainty during summarisation, which is more statistically sound.
  
  #  5. Run DESeq2
  dds <- DESeqDataSetFromTximport(
    txi,
    colData = sampleinfo,
    design = ~ condition)
  
  # Optional: pre-filter low-count genes (speeds things up)
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  #  6. Extract results
  # KO vs WT
  res <- results(dds, contrast = c("condition", "KO", "WT"))
  
  # Shrink LFC (recommended for ranking/visualization)
  library(apeglm)
  res_shrunk <- lfcShrink(dds, coef = "condition_KO_vs_WT", type = "apeglm")
  
  # Summary
  summary(res_shrunk, alpha = 0.05)
  
  # Save results
  res_df <- as.data.frame(res_shrunk)
  res_df <- res_df[order(res_df$padj), ]
  write.csv(res_df, "KO_vs_WT_DESeq2_results.csv")
  
  7. Quick QC plots
  # PCA
  vsd <- vst(dds, blind = TRUE)
  plotPCA(vsd, intgroup = "condition")
  
  # MA plot
  plotMA(res_shrunk, ylim = c(-5, 5))
  
  # Dispersion estimates
  plotDispEsts(dds)
  
  
  
  ##  Key decisions to make
  #  Question	Recommendation
  #  quant.sf vs quant.genes.sf?	Use quant.sf with tx2gene
  #  LFC shrinkage method?	apeglm (default best practice)
  #  Pre-filtering threshold?	≥10 counts in ≥3 samples
  # Significance cutoff?	padj < 0.05, |LFC| > 1
  
  
  
  #   Here's the simplest possible DESeq2 script for your setup:
  # Install if needed
  BiocManager::install(c("tximport", "DESeq2", "GenomicFeatures"))
  
  library(tximport)
  library(DESeq2)
  library(GenomicFeatures)
  
  # 1. tx2gene from your GTF
  txdb <- makeTxDbFromGFF("fgtf.gz", format = "gtf")
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
  
  # 2. Sample table — edit paths to match your folder structure
  sampleinfo <- data.frame(
    condition = factor(c("WT","WT","WT","KO","KO","KO")),
    row.names  = c("WT1","WT2","WT3","TKO1","TKO2","TKO3"))
  
  files <- c(
    WT1 = "WT_1/quant.sf",
    WT2 = "WT_2/quant.sf",
    WT3 = "WT_3/quant.sf",
    TKO1 = "TKO_1/quant.sf",
    TKO2 = "TKO_2/quant.sf",
    TKO3 = "TKO_3/quant.sf")
  
  # 3. Import
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
  
  # 4. DESeq2
  dds <- DESeqDataSetFromTximport(txi, colData = sampleinfo, design = ~ condition)
  dds <- DESeq(dds)
  
  # 5. Results
  res <- results(dds, contrast = c("condition", "KO", "WT"))
  write.csv(as.data.frame(res), "results.csv")
  
  
  ##  PREVIOUS SCRIPT STEP BY STEP
  # Step 1 — Load libraries
  
  library(tximport)
  library(DESeq2)
  library(GenomicFeatures)
  
  # Step 2 — Build tx2gene from GTF
  
#  txdb <- makeTxDbFromGFF("fgtf.gz", format = "gtf")
#  k <- keys(txdb, keytype = "TXNAME")
#  tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
  
  
  txdb <- makeTxDbFromGFF("/Users/nikolasgiannak/Desktop/Teaching/RNAseqCourse/Data/Mus_musculus.GRCm38.correct_chrom_names.102.chr.gtf.gz", format = "gtf")
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
  
  # Step 3 — Sample table
  
  sampleinfo <- data.frame(
    condition = factor(c("WT","WT","WT","KO","KO","KO")),
    row.names  = c("WT1","WT2","WT3","TKO1","TKO2","TKO3")
  )
  
  #   Step 4 — File paths
  
  files <- c(
    WT1 = "WT_1/quant.sf",
    WT2 = "WT_2/quant.sf",
    WT3 = "WT_3/quant.sf",
    TKO1 = "TKO_1/quant.sf",
    TKO2 = "TKO_2/quant.sf",
    TKO3 = "TKO_3/quant.sf"
    
    # Step 5 — Import Salmon data
    
    txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
    
    # Step 6 — Create DESeq2 object
    
    dds <- DESeqDataSetFromTximport(txi, colData = sampleinfo, design = ~ condition)
    
    # Step 7 — Run DESeq2
    
    dds <- DESeq(dds)
    
    # Step 8 — Get results and save
    
    res <- results(dds, contrast = c("condition", "KO", "WT"))
    write.csv(as.data.frame(res), "results.csv")
    
    -----------------------------------------------------------------------------------------------------------------------------------
      DIFFERENCE OF quant.sf and quant.genes.sf
    quant.sf — transcript level
    •	Each row is a transcript (e.g. ENST00000123)
    •	One gene can have multiple rows (one per isoform)
    •	More detailed, but needs to be summarised to gene level before DESeq2 (tximport does this for you via tx2gene)
    quant.genes.sf — gene level
    •	Each row is a gene (e.g. ENSG00000123)
    •	Already summarised, one row per gene
    •	Simpler, ready to use directly
    
    Which to use?
      •	Use quant.sf (transcript level). It's better because tximport summarises it to gene level while accounting for transcript-level uncertainty, which gives more accurate counts for DESeq2.
•	quant.genes.sf is fine too and simpler, but slightly less statistically rigorous.
•	For your simple analysis either works, but the script I gave you already uses quant.sf which is the recommended way.

The GTF file is needed to build the tx2gene table — a simple map that links every transcript to its gene:
  ENST00000123  →  ENSG00000001
ENST00000124  →  ENSG00000001
ENST00000125  →  ENSG00000002
...
Because quant.sf is at the transcript level, tximport needs this map to know which transcripts belong to the same gene so it can add them up into a single gene count.
If you use quant.genes.sf instead, you don't need the GTF at all, because the summarisation to gene level is already done by Salmon. That's one advantage of using the genes file.
-----------------------------------------------------------------------------------------------------------------------------------
  Here are the most common things people do with DESeq2 results:
  Filter significant genes
r
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
Volcano plot — visualise all genes at once
r
library(ggplot2)
df <- as.data.frame(res)
ggplot(df, aes(log2FoldChange, -log10(padj))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_vline(xintercept = c(-1, 1), color = "blue")
Heatmap — show expression patterns of top genes
r
library(pheatmap)
top_genes <- head(order(res$padj), 50)
pheatmap(assay(vst(dds))[top_genes, ])
PCA — check if samples cluster as expected
r
plotPCA(vst(dds), intgroup = "condition")
Filter significant genes
r
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
Volcano plot — visualise all genes at once
r
library(ggplot2)
df <- as.data.frame(res)
ggplot(df, aes(log2FoldChange, -log10(padj))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_vline(xintercept = c(-1, 1), color = "blue")
Heatmap — show expression patterns of top genes
r
library(pheatmap)
top_genes <- head(order(res$padj), 50)
pheatmap(assay(vst(dds))[top_genes, ])


OR

r
library(ComplexHeatmap)

# get top 50 significant genes
top_genes <- head(order(res$padj), 50)

# get normalised counts for those genes
mat <- assay(vst(dds))[top_genes, ]

# scale by row (so colours show relative expression)
mat_scaled <- t(scale(t(mat)))

# annotation bar on top showing WT vs KO
col_ann <- HeatmapAnnotation(
  condition = sampleinfo$condition,
  col = list(condition = c("WT" = "blue", "KO" = "red"))
)

# draw heatmap
Heatmap(
  mat_scaled,
  top_annotation = col_ann,
  name = "Z-score",
  show_row_names = TRUE,
  show_column_names = TRUE
)

It gives you:
  •	Rows = top 50 genes by adjusted p-value
•	Columns = your 6 samples
•	Colours = scaled expression (Z-score)
•	Top bar = WT vs KO annotation
Here is an example with split between up- and down-regulated genes

r
library(ComplexHeatmap)

# get significant genes
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# split into up and down
up   <- rownames(subset(sig, log2FoldChange > 1))
down <- rownames(subset(sig, log2FoldChange < -1))

# combine and get normalised counts
top_genes <- c(up, down)
mat <- assay(vst(dds))[top_genes, ]
mat_scaled <- t(scale(t(mat)))

# split vector — tells heatmap which genes are up vs down
split <- factor(c(rep("Upregulated", length(up)), 
                  rep("Downregulated", length(down))),
                levels = c("Upregulated", "Downregulated"))

# top annotation
col_ann <- HeatmapAnnotation(
  condition = sampleinfo$condition,
  col = list(condition = c("WT" = "blue", "KO" = "red"))
)

# draw heatmap
Heatmap(
  mat_scaled,
  top_annotation = col_ann,
  name = "Z-score",
  row_split = split,
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  show_row_names = FALSE,   # set TRUE if you want gene names (can get crowded)
  show_column_names = TRUE
)
The key addition is row_split = split which divides the heatmap into two blocks — upregulated genes on top, downregulated below, each with their own cluster.
Set show_row_names = TRUE if you want gene names on the side, but if you have many significant genes it gets crowded.

PCA — check if samples cluster as expected
r
plotPCA(vst(dds), intgroup = "condition")

Gene Ontology (GO) enrichment — what biological processes are affected?
  r
library(clusterProfiler)
ego <- enrichGO(gene = rownames(sig), OrgDb = org.Mm.eg.db, ont = "BP")
barplot(ego)
KEGG pathway analysis — which pathways are enriched?
  r
ekegg <- enrichKEGG(gene = rownames(sig))
dotplot(ekegg)
GSEA — gene set enrichment using the full ranked gene list (not just significant ones)
r
ranked <- res$log2FoldChange
names(ranked) <- rownames(res)
ranked <- sort(ranked, decreasing = TRUE)
gsea_res <- gseGO(ranked, OrgDb = org.Mm.eg.db, ont = "BP")
