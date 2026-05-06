# =============================================================================
# Vizualizacije DEG-ov - S. cerevisiae GSE201386
# =============================================================================

# --- 0. Namestitev paketov (zakomentiraj po prvem zagonu) --------------------
#install.packages(c("ggplot2", "dplyr", "ggrepel", "RColorBrewer"), repos = "https://cloud.r-project.org")
#install.packages("ggVennDiagram", repos = "https://cloud.r-project.org")

# --- 1. Knjižnice ------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggVennDiagram)

# --- 2. Poti -----------------------------------------------------------------
SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
RESULT_DIR <- file.path(SCRIPT_DIR, "rezultati")
DESEQ_DIR  <- file.path(RESULT_DIR, "deseq2")
VIZ_DIR    <- file.path(RESULT_DIR, "vizualizacije")
dir.create(VIZ_DIR, showWarnings = FALSE)

padj_cutoff <- 0.05
lfc_cutoff  <- 1

# --- 3. Uvoz podatkov --------------------------------------------------------
res_hypo  <- read.csv(file.path(DESEQ_DIR, "DESeq2_hypotonic_vs_control.csv"))
res_aa    <- read.csv(file.path(DESEQ_DIR, "DESeq2_aa_starvation_vs_control.csv"))
res_glyc  <- read.csv(file.path(DESEQ_DIR, "DESeq2_hypoglycemic_vs_control.csv"))
common    <- read.csv(file.path(DESEQ_DIR, "skupni_DEGs_vsi_3_pogoji.csv"))

# Odstrani ERCC
res_hypo  <- res_hypo[!grepl("^ERCC-", res_hypo$gene), ]
res_aa    <- res_aa[!grepl("^ERCC-", res_aa$gene), ]
res_glyc  <- res_glyc[!grepl("^ERCC-", res_glyc$gene), ]
common    <- common[!grepl("^ERCC-", common$gene), ]

# Signifikantni DEG-i po pogojih
sig <- function(df) df %>% filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
sig_hypo  <- sig(res_hypo)$gene
sig_aa    <- sig(res_aa)$gene
sig_glyc  <- sig(res_glyc)$gene

cat("DEG-i hipotonično:", length(sig_hypo), "\n")
cat("DEG-i AAS:        ", length(sig_aa), "\n")
cat("DEG-i hipoglik.:  ", length(sig_glyc), "\n")
cat("Skupni presek:    ", nrow(common), "\n")

# --- 4. Venn diagram ---------------------------------------------------------
venn_list <- list(
  Hipotonicno    = sig_hypo,
  AA_stradanje   = sig_aa,
  Hipoglikemicno = sig_glyc
)

p_venn <- ggVennDiagram(
  venn_list,
  label_alpha = 0,
  edge_size   = 1.2,
  label_size  = 5
) +
  scale_fill_gradient(low = "#EFF3FF", high = "#2171B5") +
  scale_color_manual(values = c("#E53935", "#FF9800", "#1E88E5")) +
  labs(title = "Overlap DEG-ov po pogojih",
       subtitle = "Filtrirano: padj < 0.05, |log2FC| > 1") +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle   = element_text(hjust = 0.5, size = 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(VIZ_DIR, "venn_diagram.png"), p_venn,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("Venn diagram shranjen.\n")

# --- 5. Bubble plot: top skupni DEG-i ----------------------------------------
# Za vsak skupni gen izračunaj povprečni -log10(padj) in povprečni |LFC|
bubble_df <- common %>%
  mutate(
    mean_padj    = (padj_hypotonic + padj_aa + padj_hypoglycemic) / 3,
    mean_lfc     = (abs(lfc_hypotonic) + abs(lfc_aa) + abs(lfc_hypoglycemic)) / 3,
    neg_log_padj = -log10(mean_padj),
    direction    = case_when(
      lfc_hypotonic > 0 & lfc_aa > 0 & lfc_hypoglycemic > 0 ~ "up v vseh",
      lfc_hypotonic < 0 & lfc_aa < 0 & lfc_hypoglycemic < 0 ~ "down v vseh",
      TRUE ~ "mešano"
    )
  ) %>%
  arrange(desc(neg_log_padj)) %>%
  slice(1:40)

p_bubble <- ggplot(bubble_df, aes(x = mean_lfc, y = neg_log_padj,
                                   size = mean_lfc, color = direction, label = gene)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(size = 2.8, max.overlaps = 30) +
  scale_color_manual(values = c("up v vseh" = "#E53935", "down v vseh" = "#1E88E5", "mešano" = "#FF9800")) +
  scale_size_continuous(range = c(2, 10)) +
  theme_bw() +
  labs(
    title    = "Top 40 skupnih DEG-ov čez vse 3 stresne pogoje",
    subtitle = "Velikost = povprečni |log2FC|, Y = statistična zanesljivost",
    x        = "Povprečni |log2FC|",
    y        = "-log10(povprečni padj)",
    color    = "Smer ekspresije",
    size     = "|log2FC|"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(VIZ_DIR, "bubble_skupni_DEGs.png"), p_bubble, width = 10, height = 7, dpi = 150)
cat("Bubble plot shranjen.\n")

# --- 6. Stripplot: LFC po pogojih za top skupne DEG-e ------------------------
# Prikaži kako se isti gen obnaša v vseh 3 pogojih
top20 <- common %>%
  mutate(mean_lfc = (abs(lfc_hypotonic) + abs(lfc_aa) + abs(lfc_hypoglycemic)) / 3) %>%
  arrange(desc(mean_lfc)) %>%
  slice(1:20) %>%
  pull(gene)

strip_df <- bind_rows(
  data.frame(gene = common$gene, lfc = common$lfc_hypotonic,  pogoj = "Hipotonicno"),
  data.frame(gene = common$gene, lfc = common$lfc_aa,          pogoj = "AA stradanje"),
  data.frame(gene = common$gene, lfc = common$lfc_hypoglycemic, pogoj = "Hipoglikemicno")
) %>% filter(gene %in% top20) %>%
  mutate(gene = factor(gene, levels = top20))

p_strip <- ggplot(strip_df, aes(x = lfc, y = gene, color = pogoj, shape = pogoj)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, alpha = 0.85, position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c(
    "Hipotonicno"   = "#2196F3",
    "AA stradanje"  = "#FF9800",
    "Hipoglikemicno" = "#E91E63"
  )) +
  theme_bw() +
  labs(
    title    = "log2FC top 20 skupnih DEG-ov po pogojih",
    subtitle = "Vsaka pika = vrednost v enem pogoju; geni sortirani po povprečnem |LFC|",
    x        = "log2 Fold Change",
    y        = NULL,
    color    = "Pogoj",
    shape    = "Pogoj"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(VIZ_DIR, "stripplot_top20_DEGs.png"), p_strip, width = 9, height = 7, dpi = 150)
cat("Stripplot shranjen.\n")

# --- 7. Heatmap LFC: top skupni DEG-i ----------------------------------------
library(pheatmap)

top50 <- common %>%
  mutate(mean_lfc = (abs(lfc_hypotonic) + abs(lfc_aa) + abs(lfc_hypoglycemic)) / 3) %>%
  arrange(desc(mean_lfc)) %>%
  slice(1:50)

lfc_mat <- as.matrix(top50[, c("lfc_hypotonic", "lfc_aa", "lfc_hypoglycemic")])
rownames(lfc_mat) <- top50$gene
colnames(lfc_mat) <- c("Hipotonicno", "AA stradanje", "Hipoglikemicno")

# Omeji vrednosti za boljšo vizualizacijo
lfc_mat_clipped <- pmin(pmax(lfc_mat, -15), 15)

png(file.path(VIZ_DIR, "heatmap_lfc_top50.png"), width = 800, height = 1200, res = 130)
pheatmap(lfc_mat_clipped,
         cluster_cols    = FALSE,
         clustering_method = "ward.D2",
         color           = colorRampPalette(c("#1E88E5", "white", "#E53935"))(100),
         breaks          = seq(-15, 15, length.out = 101),
         border_color    = NA,
         fontsize_row    = 7,
         main            = "log2FC top 50 skupnih DEG-ov\n(vrednosti omejene na ±15)")
dev.off()
cat("LFC heatmap shranjen.\n")

# --- 8. Povzetek -------------------------------------------------------------
cat("\n========== VIZUALIZACIJE SHRANJENE ==========\n")
cat("Mapa:", VIZ_DIR, "\n")
cat("  venn_diagram.png        — overlap DEG-ov med pogoji\n")
cat("  bubble_skupni_DEGs.png  — statistična moč skupnih DEG-ov\n")
cat("  stripplot_top20_DEGs.png — LFC po pogojih za top 20 genov\n")
cat("  heatmap_lfc_top50.png   — LFC heatmap top 50 skupnih DEG-ov\n")
