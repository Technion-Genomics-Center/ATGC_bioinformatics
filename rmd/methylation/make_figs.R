# Methylation (RRBS) example figures for the ATGC services page.
# Synthetic / representative data only. Mirrors the figure types in the
# Irit Ben Aharon follow-up report: per-sample DMR methylation heatmap,
# top-DMR violins, GO-BP pathway dot plot, and a pathway-gene chord.
# Output: rmd/methylation/figs/*.svg  (then copied to assets/).

setwd("C:/Users/HagayLadany/Desktop/webpage/rmd/methylation")
dir.create("figs", showWarnings = FALSE)
suppressPackageStartupMessages({
  library(ggplot2); library(svglite); library(pheatmap); library(circlize)
  library(dplyr); library(tidyr); library(scales); library(RColorBrewer)
})
set.seed(42)

grp_pal <- c(Control = "#289de0", Treatment = "#e82281")
dir_pal <- c(Hyper = "#e73334", Hypo = "#289de0")

## ---- synthetic per-sample methylation at significant DMRs ----
n_ctrl <- 4; n_dxr <- 4
samples <- c(paste0("Ctrl_", 1:n_ctrl), paste0("Treat_", 1:n_dxr))
group   <- factor(c(rep("Control", n_ctrl), rep("Treatment", n_dxr)), levels = c("Control", "Treatment"))
nD  <- 30
dir <- rep(c("Hyper", "Hypo"), length.out = nD)

M <- matrix(NA_real_, nrow = nD, ncol = length(samples),
            dimnames = list(sprintf("DMR%02d", 1:nD), samples))
delta <- numeric(nD)
for (i in 1:nD) {
  if (dir[i] == "Hyper") { base_ctrl <- runif(1, 15, 45); d <- runif(1, 18, 45) }   # gain in DXR
  else                   { base_ctrl <- runif(1, 55, 85); d <- -runif(1, 18, 45) }  # loss in DXR
  ctrl_vals <- pmin(100, pmax(0, rnorm(n_ctrl, base_ctrl, 6)))
  dxr_vals  <- pmin(100, pmax(0, rnorm(n_dxr, base_ctrl + d, 7)))
  M[i, ]    <- c(ctrl_vals, dxr_vals)
  delta[i]  <- mean(dxr_vals) - mean(ctrl_vals)
}
# a few CpG-coverage gaps (shown grey)
M[sample(length(M), 6)] <- NA

## ---- 1. heatmap (pheatmap) ----
ann_col <- data.frame(Group = group);                              rownames(ann_col) <- samples
ann_row <- data.frame(Direction = factor(dir, levels = c("Hyper","Hypo"))); rownames(ann_row) <- rownames(M)
ann_colors <- list(Group = grp_pal, Direction = dir_pal)
# cluster rows on a mean-imputed copy (NA breaks dist()), but display the NAs
M_imp <- M
for (i in 1:nD) M_imp[i, is.na(M_imp[i, ])] <- mean(M_imp[i, ], na.rm = TRUE)
hc_row <- hclust(dist(M_imp), method = "complete")

svglite("figs/meth-heatmap.svg", width = 8.0, height = 5.1)
pheatmap(M,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = seq(0, 100, length.out = 101),
         cluster_rows = hc_row, cluster_cols = FALSE,
         annotation_col = ann_col, annotation_row = ann_row, annotation_colors = ann_colors,
         show_rownames = FALSE, gaps_col = n_ctrl, na_col = "#e3e8ec",
         border_color = "#ffffff", fontsize = 9,
         main = "Per-sample % methylation at significant DMRs")
dev.off()

## ---- 2. top-DMR violins (ggplot2) ----
ord <- order(abs(delta), decreasing = TRUE)[1:2]
viol <- do.call(rbind, lapply(ord, function(i)
  data.frame(DMR = rownames(M)[i], group = group, value = M[i, ], dir = dir[i])))
viol <- viol[!is.na(viol$value), ]
viol$DMR <- factor(viol$DMR, levels = rownames(M)[ord])

p_violin <- ggplot(viol, aes(group, value, fill = group, colour = group)) +
  geom_violin(alpha = 0.20, colour = NA, width = 0.9, trim = FALSE) +
  geom_jitter(width = 0.09, size = 1.5, alpha = 0.85, stroke = 0) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "#33414c", linewidth = 0.4) +
  facet_wrap(~ DMR, nrow = 1) +
  scale_fill_manual(values = grp_pal) + scale_colour_manual(values = grp_pal) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = NULL, y = "% methylation",
       subtitle = "Top 2 DMRs by |delta methylation|") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "#eef3f7", colour = NA),
        strip.text = element_text(face = "bold", size = 9.5, colour = "#33414c"),
        plot.subtitle = element_text(size = 9.5, colour = "#6b7a86"),
        axis.text = element_text(colour = "#6b7a86"),
        axis.title = element_text(colour = "#42525e"))
svglite("figs/meth-violin.svg", width = 4.8, height = 3.5); print(p_violin); dev.off()

## ---- 3. GO-BP pathway dot plot (ggplot2) ----
paths <- data.frame(
  pathway = c("Cardiac muscle tissue development","Regulation of cardiac muscle contraction",
              "Cardiac conduction","Ventricular cardiac muscle development","Apoptotic signaling",
              "DNA methylation in gamete generation","Vascular endothelial cell proliferation",
              "Extracellular matrix organization","Response to oxidative stress",
              "Regulation of angiogenesis","Canonical NF-kB signaling","TNF cytokine production"),
  dir  = c("Hyper","Hyper","Hyper","Hyper","Hyper","Hyper","Hypo","Hypo","Hypo","Hypo","Hypo","Hypo"),
  padj = c(0.002,0.010,0.030,0.050,0.060,0.120,0.004,0.020,0.070,0.090,0.190,0.190),
  fold = c(4.8,3.6,3.1,2.7,2.8,2.3,4.2,3.4,2.9,2.4,2.1,2.0),
  hits = c(14,9,7,6,7,5,12,11,8,5,4,4))
paths$logp <- -log10(paths$padj)
paths$dir  <- factor(paths$dir, levels = c("Hyper","Hypo"))

p_dot <- ggplot(paths, aes(logp, reorder(pathway, logp))) +
  geom_vline(xintercept = -log10(0.15), linetype = "dashed", colour = "#b3c0cc") +
  geom_point(aes(size = hits, colour = fold)) +
  scale_colour_gradient(low = "#f6c6a8", high = "#d73027", name = "fold\nenrich.") +
  scale_size(range = c(2.5, 8), name = "gene hits") +
  facet_grid(dir ~ ., scales = "free_y", space = "free_y") +
  labs(x = "-log10(adj. P)", y = NULL,
       subtitle = "Curated GO-BP pathways (dashed line = adj. P 0.15)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#eef3f7", colour = NA),
        strip.text = element_text(face = "bold", colour = "#33414c"),
        plot.subtitle = element_text(size = 9.5, colour = "#6b7a86"),
        axis.text.y = element_text(size = 9, colour = "#42525e"),
        axis.text.x = element_text(colour = "#6b7a86"))
svglite("figs/meth-dotplot.svg", width = 7.4, height = 4.8); print(p_dot); dev.off()

## ---- 4. pathway-gene chord (circlize) ----
paths_c <- c("Cardiac development","NF-kB signaling","Angiogenesis",
             "Oxidative stress","ECM organization","Apoptosis")
genes   <- c("Gata4","Nppa","Myh6","Tnnt2","Rela","Nfkb1","Vegfa",
             "Kdr","Sod2","Nqo1","Col1a1","Mmp9","Casp3","Bax")
gene_dir <- c(Gata4="Hyper",Nppa="Hyper",Myh6="Hyper",Tnnt2="Hyper",Rela="Hypo",
              Nfkb1="Hypo",Vegfa="Hypo",Kdr="Hypo",Sod2="Hypo",Nqo1="Hypo",
              Col1a1="Hypo",Mmp9="Hypo",Casp3="Hyper",Bax="Hyper")
links <- list(
  "Cardiac development" = c("Gata4","Nppa","Myh6","Tnnt2"),
  "NF-kB signaling"     = c("Rela","Nfkb1","Mmp9"),
  "Angiogenesis"        = c("Vegfa","Kdr","Mmp9"),
  "Oxidative stress"    = c("Sod2","Nqo1","Bax"),
  "ECM organization"    = c("Col1a1","Mmp9","Tnnt2"),
  "Apoptosis"           = c("Casp3","Bax","Nfkb1"))
mat <- matrix(0, nrow = length(paths_c), ncol = length(genes),
              dimnames = list(paths_c, genes))
for (p in names(links)) for (g in links[[p]]) mat[p, g] <- 1 + runif(1, 0, 1.5)

path_col <- setNames(brewer.pal(length(paths_c), "Set2"), paths_c)
gene_col <- setNames(dir_pal[gene_dir[genes]], genes)
grid_col <- c(path_col, gene_col)

svglite("figs/meth-chord.svg", width = 6.4, height = 6.4)
circos.clear()
circos.par(gap.after = c(rep(2, length(paths_c) - 1), 10, rep(2, length(genes) - 1), 10),
           start.degree = 90)
chordDiagram(mat, grid.col = grid_col, transparency = 0.45,
             annotationTrack = "grid", preAllocateTracks = list(track.height = 0.10),
             link.sort = TRUE, link.decreasing = TRUE)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.4, CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.62)
}, bg.border = NA)
title("Pathway-gene chord (red = hyper, blue = hypo)", cex.main = 0.95)
circos.clear()
dev.off()

cat("DONE\n"); print(list.files("figs"))
