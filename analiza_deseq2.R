# =============================================================================
# DESeq2 analiza scRNA-seq podatkov S. cerevisiae - GSE201386
# Pogoji: kontrola vs. hipotonično, stradanje AA, hipoglikemično
# =============================================================================

# --- 0. Namestitev paketov -------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("DESeq2"), ask = FALSE)
#install.packages(c("ggplot2", "pheatmap", "ggrepel", "dplyr", "tidyr"), repos = "https://cloud.r-project.org")

# --- 1. Knjižnice ------------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(dplyr)

# --- 2. Poti -----------------------------------------------------------------
DATA_DIR    <- file.path(dirname(rstudioapi::getSourceEditorContext()$path),
                         "podatki", "GSE201386_RAW")
RESULT_DIR  <- file.path(dirname(rstudioapi::getSourceEditorContext()$path),
                         "rezultati")
DESEQ_DIR   <- file.path(RESULT_DIR, "deseq2")
dir.create(RESULT_DIR,  showWarnings = FALSE)
dir.create(DESEQ_DIR,   showWarnings = FALSE)

QC_ROWS <- c("no_feature", "ambiguous", "too_low_aQual",
             "not_aligned", "alignment_not_unique")

# --- 3. Uvoz in sestava count matrike ----------------------------------------
files <- list.files(DATA_DIR, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
cat("Najdenih datotek:", length(files), "\n")

read_counts <- function(path) {
  df <- read.table(path, sep = "\t", header = FALSE, col.names = c("gene", "count"),
                   colClasses = c("character", "character"))
  df <- df[!df$gene %in% QC_ROWS, ]
  df$count <- suppressWarnings(as.integer(df$count))
  df[!is.na(df$count), ]
}

count_list <- lapply(files, read_counts)

# Vzorčna imena: del med GSM\d+_ in _QC
sample_names <- sub(".*_([^_]+)_QC.*", "\\1", basename(files))
names(count_list) <- sample_names

# Referenčni seznam genov iz prvega vzorca
all_genes <- count_list[[1]]$gene

count_matrix <- sapply(count_list, function(df) {
  counts <- df$count[match(all_genes, df$gene)]
  ifelse(is.na(counts), 0L, counts)
})
rownames(count_matrix) <- all_genes

cat("Dimenzije matrike:", nrow(count_matrix), "genov x", ncol(count_matrix), "vzorcev\n")

# --- 4. Metadata / colData ---------------------------------------------------
assign_condition <- function(name) {
  if (grepl("^[AEG]\\d|^Y7", name))        return("control")
  if (grepl("^[CF]",          name))        return("hypotonic")
  if (grepl("^D",             name))        return("aa_starvation")
  if (grepl("^Y[25]",         name))        return("hypoglycemic")
  return(NA_character_)
}

col_data <- data.frame(
  sample    = sample_names,
  condition = factor(sapply(sample_names, assign_condition),
                     levels = c("control", "hypotonic", "aa_starvation", "hypoglycemic")),
  row.names = sample_names
)

# Preverimo da ni nerazvrščenih vzorcev
unassigned <- col_data$sample[is.na(col_data$condition)]
if (length(unassigned) > 0)
  warning("Nerazvrščeni vzorci: ", paste(unassigned, collapse = ", "))

cat("Vzorci po pogojih:\n")
print(table(col_data$condition))

# --- 5. DESeq2 ---------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = col_data,
  design    = ~condition
)

# Osnoven prefilter: vsaj 10 countov v vsaj 3 vzorcih
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("Genov po filtriranju:", nrow(dds), "\n")

dds <- DESeq(dds)

# --- 6. Rezultati za vsako primerjavo ----------------------------------------
padj_cutoff <- 0.05
lfc_cutoff  <- 1

get_results <- function(dds, contrast_name, numerator) {
  res <- results(dds, contrast = c("condition", numerator, "control"),
                 alpha = padj_cutoff)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$comparison <- contrast_name
  res_df
}

res_hypotonic    <- get_results(dds, "hypotonic_vs_control",    "hypotonic")
res_aa           <- get_results(dds, "aa_starvation_vs_control", "aa_starvation")
res_hypoglycemic <- get_results(dds, "hypoglycemic_vs_control", "hypoglycemic")

# Shrani CSV
write.csv(res_hypotonic,    file.path(DESEQ_DIR, "DESeq2_hypotonic_vs_control.csv"),    row.names = FALSE)
write.csv(res_aa,           file.path(DESEQ_DIR, "DESeq2_aa_starvation_vs_control.csv"), row.names = FALSE)
write.csv(res_hypoglycemic, file.path(DESEQ_DIR, "DESeq2_hypoglycemic_vs_control.csv"), row.names = FALSE)

# --- 7. Skupni DEG-i (presek vseh 3 primerjav) --------------------------------
sig_genes <- function(res_df) {
  res_df %>%
    filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff) %>%
    pull(gene)
}

sig_hypo  <- sig_genes(res_hypotonic)
sig_aa    <- sig_genes(res_aa)
sig_glyc  <- sig_genes(res_hypoglycemic)

common_degs_all <- Reduce(intersect, list(sig_hypo, sig_aa, sig_glyc))
# Odstrani ERCC spike-ine — niso yeast geni
common_degs <- common_degs_all[!grepl("^ERCC-", common_degs_all)]
cat("Skupnih DEG-ov (vsi 3 pogoji, brez ERCC):", length(common_degs), "\n")

# Tabela skupnih DEG-ov z LFC iz vseh primerjav
common_table <- data.frame(gene = common_degs) %>%
  left_join(res_hypotonic    %>% select(gene, log2FoldChange, padj) %>% rename(lfc_hypotonic = log2FoldChange, padj_hypotonic = padj),    by = "gene") %>%
  left_join(res_aa           %>% select(gene, log2FoldChange, padj) %>% rename(lfc_aa = log2FoldChange, padj_aa = padj),                   by = "gene") %>%
  left_join(res_hypoglycemic %>% select(gene, log2FoldChange, padj) %>% rename(lfc_hypoglycemic = log2FoldChange, padj_hypoglycemic = padj), by = "gene")

write.csv(common_table, file.path(DESEQ_DIR, "skupni_DEGs_vsi_3_pogoji.csv"), row.names = FALSE)
cat("Tabela skupnih DEG-ov shranjena.\n")

# --- 8. Vizualizacije --------------------------------------------------------

# 8a. PCA plot
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  scale_color_manual(values = c(control = "#4CAF50", hypotonic = "#2196F3",
                                aa_starvation = "#FF9800", hypoglycemic = "#E91E63")) +
  theme_bw() +
  ggtitle("PCA — vzorci po pogojih")

ggsave(file.path(DESEQ_DIR, "PCA_pogoji.png"), p_pca, width = 7, height = 5, dpi = 150)

# 8b. Volcano ploti
make_volcano <- function(res_df, title) {
  df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      sig = case_when(
        padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "up",
        padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "down",
        TRUE ~ "ns"
      ),
      label = ifelse(sig != "ns" & abs(log2FoldChange) > 3 & padj < 1e-10, gene, NA)
    )

  counts_sig <- table(df$sig)

  ggplot(df, aes(log2FoldChange, -log10(padj), color = sig)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
    scale_color_manual(values = c(up = "#E53935", down = "#1E88E5", ns = "grey70"),
                       labels = c(up   = paste0("Up (", counts_sig["up"],   ")"),
                                  down = paste0("Down (", counts_sig["down"], ")"),
                                  ns   = "NS")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_cutoff),        linetype = "dashed", color = "grey40") +
    theme_bw() +
    labs(title = title, x = "log2 Fold Change", y = "-log10(padj)", color = NULL)
}

p_v1 <- make_volcano(res_hypotonic,    "Volcano: hipotonično vs. kontrola")
p_v2 <- make_volcano(res_aa,           "Volcano: stradanje AA vs. kontrola")
p_v3 <- make_volcano(res_hypoglycemic, "Volcano: hipoglikemično vs. kontrola")

ggsave(file.path(DESEQ_DIR, "volcano_hypotonic.png"),    p_v1, width = 7, height = 5, dpi = 150)
ggsave(file.path(DESEQ_DIR, "volcano_aa_starvation.png"), p_v2, width = 7, height = 5, dpi = 150)
ggsave(file.path(DESEQ_DIR, "volcano_hypoglycemic.png"),  p_v3, width = 7, height = 5, dpi = 150)

# 8c. Heatmap skupnih DEG-ov
if (length(common_degs) > 1) {
  mat <- assay(vsd)[common_degs, ]
  control_mean <- rowMeans(mat[, col_data$condition == "control"])
  mat <- mat - control_mean

  annotation_col <- data.frame(condition = col_data$condition, row.names = colnames(mat))

  ann_colors <- list(condition = c(control = "#4CAF50", hypotonic = "#2196F3",
                                   aa_starvation = "#FF9800", hypoglycemic = "#E91E63"))

  # Omeji na max 100 genov za preglednost
  if (nrow(mat) > 100) {
    top_genes <- common_table %>%
      mutate(mean_lfc = (abs(lfc_hypotonic) + abs(lfc_aa) + abs(lfc_hypoglycemic)) / 3) %>%
      arrange(desc(mean_lfc)) %>%
      slice(1:100) %>%
      pull(gene)
    mat <- mat[top_genes, ]
  }

  png(file.path(DESEQ_DIR, "heatmap_skupni_DEGs.png"),
      width = 1200, height = max(800, nrow(mat) * 12), res = 150)
  pheatmap(mat,
           annotation_col  = annotation_col,
           annotation_colors = ann_colors,
           show_rownames   = nrow(mat) <= 60,
           show_colnames   = FALSE,
           clustering_method = "ward.D2",
           color = colorRampPalette(c("#1E88E5", "white", "#E53935"))(100),
           main  = paste0("Skupni DEG-i (", nrow(mat), " genov)"))
  dev.off()
} else {
  cat("Premalo skupnih DEG-ov za heatmap.\n")
}

# --- 9. Povzetek -------------------------------------------------------------
cat("\n========== POVZETEK ==========\n")
cat("Signifikantni DEG-i (padj <", padj_cutoff, ", |LFC| >", lfc_cutoff, "):\n")
cat("  Hipotonično:        ", length(sig_hypo), "\n")
cat("  Stradanje AA:       ", length(sig_aa),   "\n")
cat("  Hipoglikemično:     ", length(sig_glyc), "\n")
cat("  SKUPNI (presek):    ", length(common_degs), "\n")
cat("===============================\n")
cat("Rezultati shranjeni v:", RESULT_DIR, "\n")
