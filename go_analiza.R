# =============================================================================
# GO enrichment analiza skupnih DEG-ov - S. cerevisiae GSE201386
# =============================================================================

# --- 0. Namestitev paketov (zakomentiraj po prvem zagonu) --------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("clusterProfiler", "org.Sc.sgd.db"), ask = FALSE)

# --- 1. Knjižnice ------------------------------------------------------------
library(clusterProfiler)
library(org.Sc.sgd.db)
library(ggplot2)
library(dplyr)

# --- 2. Poti -----------------------------------------------------------------
SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
RESULT_DIR <- file.path(SCRIPT_DIR, "rezultati")
GO_DIR     <- file.path(RESULT_DIR, "GO")
dir.create(GO_DIR, showWarnings = FALSE)

# --- 3. Uvoz skupnih DEG-ov --------------------------------------------------
common_table <- read.csv(file.path(RESULT_DIR, "skupni_DEGs_vsi_3_pogoji.csv"))

# Odstrani ERCC spike-ine
common_table <- common_table[!grepl("^ERCC-", common_table$gene), ]
cat("Skupnih yeast DEG-ov:", nrow(common_table), "\n")

# Celoten seznam genov iz ene od DESeq2 tabel (ozadje za enrichment)
background <- read.csv(file.path(RESULT_DIR, "DESeq2_aa_starvation_vs_control.csv"))
background <- background$gene[!grepl("^ERCC-", background$gene)]
cat("Ozadje (vsi testirani geni):", length(background), "\n")

# --- 4. GO enrichment --------------------------------------------------------
run_go <- function(gene_list, ont, label) {
  # Najprej brez cutoffa da dobimo vse rezultate
  res <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Sc.sgd.db,
    keyType       = "ORF",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1.0,
    qvalueCutoff  = 1.0,
    minGSSize     = 3,
    readable      = FALSE
  )
  df <- as.data.frame(res)
  if (nrow(df) == 0) {
    cat(label, "— ni rezultatov\n")
    return(res)
  }
  # Filtriraj na pvalue < 0.01 (nekorigirano, eksploratorni rezultat)
  df_sig <- df[df$pvalue < 0.01, ]
  cat(label, "— terminov z pvalue < 0.01:", nrow(df_sig),
      " | z p.adjust < 0.05:", sum(df$p.adjust < 0.05), "\n")
  # Vrni filtriran data frame kot enrichResult
  res@result <- df_sig
  res
}

deg_genes <- common_table$gene

# Debug: preveri mapping genov v SGD bazo
mapped <- bitr(deg_genes, fromType = "ORF", toType = c("GENENAME", "GO"), OrgDb = org.Sc.sgd.db)
cat("Genov ki se mapirajo v SGD:", length(unique(mapped$ORF)), "od", length(deg_genes), "\n")

# Debug: poglej surove rezultate brez cutoffov
go_bp_raw <- enrichGO(
  gene          = deg_genes,
  OrgDb         = org.Sc.sgd.db,
  keyType       = "ORF",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1.0,
  qvalueCutoff  = 1.0,
  minGSSize     = 3,
  readable      = FALSE
)
df_raw <- as.data.frame(go_bp_raw)
cat("BP terminov z pvalueCutoff=1:", nrow(df_raw), "\n")
if (nrow(df_raw) > 0) {
  cat("Najboljši p.adjust vrednosti:\n")
  print(head(df_raw[order(df_raw$p.adjust), c("Description", "pvalue", "p.adjust", "Count")], 10))
  cat("Minimalni p.adjust:", min(df_raw$p.adjust), "\n")
  cat("Minimalni pvalue:", min(df_raw$pvalue), "\n")
  cat("Število terminov z pvalue < 0.05:", sum(df_raw$pvalue < 0.05), "\n")
  cat("Število terminov z p.adjust < 0.05:", sum(df_raw$p.adjust < 0.05), "\n")
}

go_bp <- run_go(deg_genes, "BP", "Biological Process")
go_mf <- run_go(deg_genes, "MF", "Molecular Function")
go_cc <- run_go(deg_genes, "CC", "Cellular Component")

# Shrani CSV — uporabi nekorigirane če so korigirani prazni
save_go <- function(res, path) {
  df <- as.data.frame(res)
  if (nrow(df) > 0) write.csv(df, path, row.names = FALSE)
}
save_go(go_bp, file.path(GO_DIR, "GO_BP.csv"))
save_go(go_mf, file.path(GO_DIR, "GO_MF.csv"))
save_go(go_cc, file.path(GO_DIR, "GO_CC.csv"))

# --- 5. Vizualizacije --------------------------------------------------------
plot_go <- function(go_res, title, filename) {
  df <- as.data.frame(go_res)
  if (nrow(df) == 0) {
    cat("Ni enrichanih terminov za:", title, "\n")
    return(invisible(NULL))
  }

  df <- df %>%
    arrange(p.adjust) %>%
    slice(1:min(20, nrow(.))) %>%
    mutate(
      Description = factor(Description, levels = rev(Description)),
      GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text = x)))
    )

  p <- ggplot(df, aes(x = GeneRatio_num, y = Description, color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "#E53935", high = "#90CAF9", name = "p.adjust") +
    scale_size_continuous(name = "Število genov") +
    theme_bw() +
    labs(title = title, x = "Gene Ratio", y = NULL) +
    theme(axis.text.y = element_text(size = 9))

  ggsave(file.path(GO_DIR, filename), p,
         width = 9, height = max(4, nrow(df) * 0.35 + 2), dpi = 150)
}

plot_go(go_bp, "GO Biological Process — skupni DEG-i", "GO_BP_dotplot.png")
plot_go(go_mf, "GO Molecular Function — skupni DEG-i", "GO_MF_dotplot.png")
plot_go(go_cc, "GO Cellular Component — skupni DEG-i", "GO_CC_dotplot.png")

# --- 6. KEGG pathway analiza -------------------------------------------------
go_kegg <- enrichKEGG(
  gene          = deg_genes,
  organism      = "sce",
  universe      = background,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  minGSSize     = 3
)
cat("KEGG — enrichanih pathwayev:", nrow(as.data.frame(go_kegg)), "\n")
write.csv(as.data.frame(go_kegg), file.path(GO_DIR, "KEGG_pathways.csv"), row.names = FALSE)
plot_go(go_kegg, "KEGG Pathways — skupni DEG-i", "KEGG_dotplot.png")

# --- 7. Povzetek -------------------------------------------------------------
cat("\n========== POVZETEK GO ==========\n")
cat("Biological Process enrichanih terminov:", nrow(as.data.frame(go_bp)), "\n")
cat("Molecular Function enrichanih terminov:", nrow(as.data.frame(go_mf)), "\n")
cat("Cellular Component enrichanih terminov:", nrow(as.data.frame(go_cc)), "\n")
cat("KEGG pathways:                         ", nrow(as.data.frame(go_kegg)), "\n")
cat("Rezultati shranjeni v:", GO_DIR, "\n")
