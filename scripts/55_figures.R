#!/usr/bin/env Rscript
# =============================================================================
# 55_figures.R  -  all data figures of the MAP manuscript (mSystems)
#
#   Fig2   The MAP pangenome is closed
#          A pangenome accumulation curves, four gene-presence methods
#          B pangenome composition (core / soft core / shell / cloud), donuts
#          C spurious genes in the within-study negative control
#          D genes per genome against assembly N50
#   Fig3   Lineage-specific gene content
#          A blocks of lineage-specific genes on K-10 and on Telford
#          B COG categories, genes and independent blocks
#   Fig4   Population structure (the only tree figure)
#          A circular cladogram of 1,792 genomes, rings: lineage, region, host
#          B regional composition of each lineage group
#          C 12-SNP typing scheme: marker-allele frequency per group
#   Fig5   Homoplasy and convergence
#          A changes per SNP allele, assembly calls vs read-verified calls
#          B explanations of the 21 SNPs first flagged as homoplastic
#          C gene-level convergence along the K-10 chromosome
#          D independent amino-acid changes in icl (lollipop plot)
#   Fig6   Sampling structure
#          A root-to-tip regression, B date permutations (ridges),
#          C within-study slopes (forest plot), D origin prediction (dumbbell)
#   FigS2  Scoary associations are batch artefacts
#   FigS3  Root-to-tip regression within single studies
#
# Fig 1 (workflow) and Fig S1 (cohort flow) are TikZ drawings in the
# manuscript folder and use the same colours.
#
# Usage   Rscript 55_figures.R [MAP_root] [out_dir]
#         defaults: /cfs/earth/scratch/tabasran/MAP  and  <root>/10_figures
# Needs   R >= 4.1; CRAN: ggplot2 (>= 3.4), patchwork, data.table, ape,
#         ggrepel, ggridges, scales
#   conda activate r_env
#   Rscript -e 'install.packages(c("ggplot2","patchwork","data.table","ape",
#                "ggrepel","ggridges","scales"), repos="https://cloud.r-project.org")'
# Input   files written by scripts 35/36, 46, 49, 50, 52 (paths below)
# Output  <out>/Fig2.pdf ... FigS3.pdf (vector, fonts embedded by cairo) and
#         a 600 dpi PNG of each for checking
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table); library(ape)
  library(ggrepel); library(ggridges); library(scales)
})

# unicode labels (\u2265, \u03b3, \u03c1) need a UTF-8 locale
if (!isTRUE(l10n_info()$`UTF-8`))
  for (lc in c("C.UTF-8", "en_US.UTF-8", "C.utf8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", lc)))) break
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "/cfs/earth/scratch/tabasran/MAP"
OUT  <- if (length(args) >= 2) args[2] else file.path(ROOT, "10_figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
f <- function(...) file.path(ROOT, ...)
set.seed(20261029)

# =============================================================================
# 1. COLOURS
# =============================================================================
# global_colors as given by the user
global_colors <- list(
  "Core (>=99%)"       = "#3D5A80",
  "Soft-Core (95-99%)" = "#81B29A",
  "Shell (15-95%)"     = "#F4A261",
  "Cloud (<15%)"       = "#E07A5F",
  LSDV = "#3D5A80", SPPV = "#E07A5F", GTPV = "#81B29A",
  Background = "white")

# colours measured in the reference figures (Figure 1-5, Fig S1-S5)
REF <- c(navy   = "#002F57",  # heatmap presence, map LSDV
         slate  = "#495C7C",  # ring 'viral species'
         steel  = "#90A0BC",  # core / LSDV / synonymous / 'no recombination'
         wine   = "#871746",  # shell / SPPV / non-synonymous / threshold line
         rose   = "#CB9A9F",  # cloud / SPPV ring
         blush  = "#E0D4DE",  # confidence band (Figure 5C)
         olive  = "#585D26",  # soft core / GTPV
         moss   = "#3F482C",  # Africa
         tan    = "#B69065",  # Europe / vaccine
         sand   = "#ECD7BC",  # Asia / wild
         apricot = "#F8BA79", # GTPV on the map
         sea    = "#D6EAF8")  # map background
NEUTRAL <- c(ink = "#222222", mid = "#6E6E6E", light = "#D9D9D9", faint = "#F2F2F2")

# Pangenome categories. Your list (global_colors) and Figure 1 of the
# reference set use different colours for the same categories; set
# PANGENOME_COLOURS to "reference" to use the Figure 1 colours instead.
PANGENOME_COLOURS <- "global_colors"
COL_PAN <- if (PANGENOME_COLOURS == "global_colors") {
  c(Core = global_colors[["Core (>=99%)"]], `Soft core` = global_colors[["Soft-Core (95-99%)"]],
    Shell = global_colors[["Shell (15-95%)"]], Cloud = global_colors[["Cloud (<15%)"]])
} else {
  c(Core = REF[["steel"]], `Soft core` = REF[["olive"]], Shell = REF[["wine"]], Cloud = REF[["rose"]])
}
PAN_LABEL <- c(Core = "Core (\u226599%)", `Soft core` = "Soft core (95\u201399%)",
               Shell = "Shell (15\u201395%)", Cloud = "Cloud (<15%)")

# gene-presence methods, least to most reliable (warm -> cool)
METHODS <- c("Annotation", "DNA alignment", "Reads, two states", "Reads, three states")
COL_METHOD <- setNames(c(global_colors[["Cloud (<15%)"]], global_colors[["Shell (15-95%)"]],
                         global_colors[["Soft-Core (95-99%)"]], global_colors[["Core (>=99%)"]]), METHODS)

GROUPS <- c("S/Type I", "S/Type III", "Bison", "C basal grade", "C clade D")
COL_GROUP <- setNames(c(REF[["wine"]], REF[["rose"]], REF[["olive"]], REF[["steel"]], REF[["navy"]]), GROUPS)
REGIONS <- c("Americas", "Europe", "Oceania", "Other/unknown")
COL_REGION <- setNames(c(REF[["moss"]], REF[["tan"]], REF[["apricot"]], NEUTRAL[["light"]]), REGIONS)
HOSTS <- c("Cattle", "Sheep", "Other ruminants", "Other/unknown")
COL_HOST <- setNames(c(REF[["slate"]], global_colors[["Cloud (<15%)"]], global_colors[["Soft-Core (95-99%)"]],
                       NEUTRAL[["light"]]), HOSTS)
COG_NAME <- c(J = "Translation, ribosomal structure", A = "RNA processing and modification",
  K = "Transcription", L = "Replication, recombination, repair", B = "Chromatin structure",
  D = "Cell cycle control, cell division", Y = "Nuclear structure", V = "Defense mechanisms",
  T = "Signal transduction", M = "Cell wall/membrane/envelope biogenesis", N = "Cell motility",
  Z = "Cytoskeleton", W = "Extracellular structures", U = "Intracellular trafficking, secretion",
  O = "Posttranslational modification, chaperones", X = "Mobilome: prophages, transposons",
  C = "Energy production and conversion", G = "Carbohydrate transport and metabolism",
  E = "Amino acid transport and metabolism", F = "Nucleotide transport and metabolism",
  H = "Coenzyme transport and metabolism", I = "Lipid transport and metabolism",
  P = "Inorganic ion transport and metabolism",
  Q = "Secondary metabolites biosynthesis, transport", R = "General function prediction only",
  S = "Function unknown", `-` = "No COG category or no hit")

# two-state contrasts (as Figure 5A of the reference set)
COL_PAIR <- c(REF[["steel"]], REF[["wine"]])

# =============================================================================
# 2. THEME
# =============================================================================
FONT <- "serif"          # the reference figures use a Times-like serif; "sans" for Helvetica
BASE <- 7.5              # pt at final size (17.4 cm wide)
theme_map <- function(base = BASE, family = FONT, grid = "xy") {
  t <- theme_bw(base_size = base, base_family = family) %+replace%
    theme(
      panel.border       = element_rect(colour = NEUTRAL[["ink"]], fill = NA, linewidth = 0.45),
      panel.background   = element_rect(fill = "white", colour = NA),
      panel.grid.major   = element_line(colour = "#E6E6E6", linewidth = 0.25),
      panel.grid.minor   = element_blank(),
      axis.title         = element_text(face = "bold", size = base, colour = NEUTRAL[["ink"]]),
      axis.title.x       = element_text(face = "bold", margin = margin(t = 3)),
      axis.title.y       = element_text(face = "bold", angle = 90, margin = margin(r = 3)),
      axis.text          = element_text(size = base - 1, colour = "#333333"),
      axis.ticks         = element_line(colour = NEUTRAL[["ink"]], linewidth = 0.3),
      axis.ticks.length  = unit(1.2, "mm"),
      legend.position    = "bottom",
      legend.direction   = "horizontal",
      legend.title       = element_text(face = "bold", size = base - 0.5),
      legend.text        = element_text(size = base - 1),
      legend.key         = element_blank(),
      legend.background  = element_blank(),
      legend.key.size    = unit(3, "mm"),
      legend.spacing.x   = unit(1.2, "mm"),
      legend.margin      = margin(0, 0, 0, 0),
      legend.box.margin  = margin(-2, 0, 0, 0),
      plot.title         = element_text(face = "bold", size = base + 0.5, hjust = 0, margin = margin(b = 4)),
      plot.subtitle      = element_text(size = base - 0.5, colour = NEUTRAL[["mid"]], hjust = 0, margin = margin(b = 3)),
      plot.tag           = element_text(face = "bold", size = base + 5, family = family),
      plot.tag.position  = c(0, 1),
      strip.background   = element_rect(fill = "#2D2D2D", colour = "#2D2D2D"),
      strip.text         = element_text(colour = "white", face = "bold", size = base - 0.5, margin = margin(2, 2, 2, 2)),
      plot.margin        = margin(6, 6, 4, 6))
  if (!grepl("x", grid)) t <- t + theme(panel.grid.major.x = element_blank())
  if (!grepl("y", grid)) t <- t + theme(panel.grid.major.y = element_blank())
  t
}
theme_void_map <- function(base = BASE, family = FONT)
  theme_void(base_size = base, base_family = family) +
  theme(plot.tag = element_text(face = "bold", size = base + 5, family = family),
        plot.tag.position = c(0, 1), plot.margin = margin(4, 4, 4, 4),
        legend.position = "none")
theme_set(theme_map())
update_geom_defaults("text", list(family = FONT, size = (BASE - 1) / .pt))
update_geom_defaults("label", list(family = FONT, size = (BASE - 1) / .pt))
TXT <- (BASE - 1) / .pt                 # geom_text size equal to axis text
tags <- function(p) p + plot_annotation(tag_levels = "A",
                                        theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

save_fig <- function(p, name, h_cm, w_cm = 17.4) {
  dev <- if (capabilities("cairo")) cairo_pdf else pdf
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w_cm, height = h_cm, units = "cm", device = dev, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w_cm, height = h_cm, units = "cm", dpi = 600, bg = "white")
  message("  written ", name, " (", w_cm, " x ", h_cm, " cm)")
}
exists_or <- function(p) { if (!file.exists(p)) message("  [!] missing ", p); file.exists(p) }
lab_n <- function(x) format(x, big.mark = ",", trim = TRUE)

# Manual legend block: several titled colour keys in one panel. Used where
# one panel carries several independent colour scales (tree rings).
legend_block <- function(blocks, ncol = 1, key = 0.9, gap = 0.5) {
  rows <- list(); y <- 0; col <- 0; x0 <- 0; per <- ceiling(length(blocks) / ncol); k <- 0
  for (nm in names(blocks)) {
    if (k == per) { k <- 0; col <- col + 1; y <- 0 }
    k <- k + 1
    v <- blocks[[nm]]
    rows[[length(rows) + 1]] <- data.table(x = col, y = y, lab = nm, fill = NA_character_, title = TRUE)
    y <- y - 1.15
    for (i in seq_along(v)) {
      rows[[length(rows) + 1]] <- data.table(x = col, y = y, lab = names(v)[i], fill = v[[i]], title = FALSE)
      y <- y - 1
    }
    y <- y - gap
  }
  d <- rbindlist(rows)
  ggplot(d) +
    geom_tile(data = d[title == FALSE], aes(x = x * 6, y = y, fill = I(fill)), width = key, height = key * 0.85,
              colour = "grey35", linewidth = 0.2) +
    geom_text(data = d[title == FALSE], aes(x = x * 6 + 0.8, y = y, label = lab), hjust = 0, size = TXT) +
    geom_text(data = d[title == TRUE], aes(x = x * 6 - 0.5, y = y, label = lab), hjust = 0, size = TXT * 1.05,
              fontface = "bold") +
    scale_x_continuous(limits = c(-0.6, max(d$x) * 6 + 5.5), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(add = 0.6)) +
    theme_void_map()
}

# =============================================================================
# 3. METADATA
# =============================================================================
read_meta <- function() {
  lin <- fread(f("05_phylogeny/typing/lineage_assignments.tsv"), colClasses = "character")
  sp  <- fread(f("05_phylogeny/typing/split_assignments.tsv"), colClasses = "character")
  tr  <- fread(f("06_pangwas/trait.csv"), colClasses = "character")
  m <- lin[, .(genome, lineage, sublineage)]
  m[, group := fifelse(lineage == "S", fifelse(sublineage == "III", "S/Type III", "S/Type I"),
                       fifelse(sublineage == "B", "Bison", "C basal grade"))]
  dside <- sp$genome[startsWith(sp$tree_side, "D ")]
  m[genome %in% dside, group := "C clade D"]
  m[, region := "Other/unknown"]
  for (r in c("Americas", "Europe", "Oceania")) {
    cc <- paste0("Region_", r)
    if (cc %in% names(tr)) m[genome %in% tr$Name[tr[[cc]] == "1"], region := r]
  }
  m[, host := "Other/unknown"]
  hc <- setNames(names(tr), tolower(names(tr)))
  for (k in list(c("hostother_ruminants", "Other ruminants"), c("hostsheep", "Sheep"), c("hostcattle", "Cattle"))) {
    if (k[1] %in% names(hc)) m[genome %in% tr$Name[tr[[hc[[k[1]]]]] == "1"], host := k[2]]
  }
  m[, group := factor(group, GROUPS)]
  m[, region := factor(region, REGIONS)]
  m[, host := factor(host, HOSTS)]
  message(sprintf("  metadata: %d genomes | %s", nrow(m),
                  paste(names(table(m$group)), table(m$group), sep = " ", collapse = ", ")))
  m[]
}

# =============================================================================
# 4. FIGURE 2  -  closed pangenome
# =============================================================================
read_presence <- function(path, ternary = FALSE) {
  x <- fread(path)
  g <- x[[1]]; x[, 1 := NULL]
  M <- as.matrix(x)
  if (!is.numeric(M)) storage.mode(M) <- "numeric"
  rownames(M) <- g
  M
}

accumulation <- function(M, n_perm = 50, step = 8) {   # 50 orders, as in Table 1
  P <- M > 0.5; P[is.na(P)] <- FALSE
  n <- ncol(P); ks <- unique(c(1:10, seq(10, n, by = step), n))
  pans <- sapply(seq_len(n_perm), function(i) {
    o <- sample.int(n)
    first <- apply(P[, o, drop = FALSE], 1, function(v) { w <- which(v); if (length(w)) w[1] else NA_integer_ })
    cs <- cumsum(tabulate(first[!is.na(first)], nbins = n))
    cs[ks]
  })
  d <- data.table(k = ks, mean = rowMeans(pans), lo = apply(pans, 1, min), hi = apply(pans, 1, max))
  fit <- lm(log(mean) ~ log(k), data = d[k >= 10])
  list(curve = d, gamma = unname(coef(fit)[2]))
}

categories <- function(M) {
  fr <- rowMeans(M > 0.5, na.rm = TRUE)
  fr <- fr[!is.na(fr) & fr > 0]
  cat <- cut(fr, c(-Inf, 0.15, 0.95, 0.99, Inf), right = FALSE, labels = c("Cloud", "Shell", "Soft core", "Core"))
  as.data.table(table(category = cat))
}

fig2 <- function(meta) {
  src <- list(
    "Annotation"          = f("04_pangenome/panaroo/gene_presence_absence.Rtab"),
    "DNA alignment"       = f("04_pangenome/dna_presence/presence_final.Rtab"),
    "Reads, two states"   = f("08_reads/summary/presence_reads.Rtab"),
    "Reads, three states" = f("08_reads/summary/presence_reads_ternary.tsv"))
  src <- src[vapply(src, exists_or, logical(1))]
  mats <- lapply(src, read_presence)
  # same genomes in every matrix
  ids <- Reduce(intersect, lapply(mats, colnames))
  mats <- lapply(mats, function(M) M[, ids, drop = FALSE])
  message("  presence matrices: ", paste(names(mats), sapply(mats, nrow), collapse = "; "), " genes; ",
          length(ids), " genomes")

  # ---- A accumulation
  acc <- lapply(mats, accumulation)
  A <- rbindlist(lapply(names(acc), function(m) acc[[m]]$curve[, method := m]))
  A[, method := factor(method, METHODS)]
  gam <- sapply(acc, `[[`, "gamma")
  message("  Heaps gamma: ", paste(names(gam), sprintf("%.4f", gam), collapse = "; "),
          "  (compare Table 1: 0.029 / 0.004 / 0.002)")
  ends <- A[, .SD[which.max(k)], by = method]
  ends[, lab := sprintf("%s\n\u03b3 = %.3f", method, gam[as.character(method)])]
  pA <- ggplot(A, aes(k, mean, colour = method, fill = method)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_text_repel(data = ends, aes(x = k, y = mean, label = lab), hjust = 0, lineheight = 0.9, family = FONT,
                    size = TXT, show.legend = FALSE, direction = "y", nudge_x = max(A$k) * 0.05,
                    xlim = c(max(A$k) * 1.04, NA), segment.size = 0.25, segment.colour = "grey55",
                    min.segment.length = 0, box.padding = 0.15, seed = 1) +
    scale_colour_manual(values = COL_METHOD, guide = "none") +
    scale_fill_manual(values = COL_METHOD, guide = "none") +
    scale_x_continuous(labels = label_comma(), expand = expansion(mult = c(0.01, 0.30))) +
    scale_y_continuous(labels = label_comma()) +
    labs(x = "Number of genomes", y = "Pangenome size (gene clusters)") +
    coord_cartesian(clip = "off")

  # ---- B composition donuts
  C <- rbindlist(lapply(names(mats), function(m) categories(mats[[m]])[, method := m]))
  C[, category := factor(category, rev(c("Core", "Soft core", "Shell", "Cloud")))]
  C[, method := factor(method, METHODS)]
  setorder(C, method, category)
  C[, `:=`(frac = N / sum(N), total = sum(N)), by = method]
  C[, `:=`(ymax = cumsum(frac)), by = method]
  C[, ymin := ymax - frac]
  message("  categories (core/soft/shell/cloud):\n",
          paste(capture.output(print(dcast(C, method ~ category, value.var = "N"))), collapse = "\n"))
  acc_lab <- C[category != "Core", .(acc = sum(N), total = total[1]), by = method]
  pB <- ggplot(C) +
    geom_rect(aes(xmin = 2.4, xmax = 4, ymin = ymin, ymax = ymax, fill = category), colour = "white", linewidth = 0.3) +
    geom_text(data = acc_lab, aes(x = 0, y = 0, label = sprintf("%s\naccessory", lab_n(acc))),
              size = TXT * 1.05, fontface = "bold", lineheight = 0.9) +
    coord_polar(theta = "y") +
    facet_wrap(~method, nrow = 2, labeller = label_wrap_gen(14)) +
    scale_x_continuous(limits = c(0, 4)) +
    scale_fill_manual(values = COL_PAN, breaks = names(PAN_LABEL), labels = PAN_LABEL, name = NULL) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_void_map() +
    theme(legend.position = "bottom", legend.text = element_text(size = BASE - 1),
          legend.key.size = unit(2.8, "mm"),
          strip.text = element_text(face = "bold", size = BASE - 0.5, margin = margin(b = 1)),
          panel.spacing = unit(1, "mm"))

  # ---- C negative control (values from 26_call_validation, 35, 36 logs)
  FP <- data.table(method = factor(METHODS, METHODS),
                   `Split at study median N50` = c(278, 115, 212, 23),
                   `Split at study median read depth` = c(NA, 63, 260, 34))
  FPl <- melt(FP, id.vars = "method", variable.name = "control", value.name = "genes")
  FPl[, side := if (.N == 2) ifelse(genes == min(genes, na.rm = TRUE), 1, 0) else 0, by = method]
  FPl[is.na(side), side := 0]
  pC <- ggplot(FPl[!is.na(genes)], aes(genes, method)) +
    geom_line(aes(group = method), colour = "grey70", linewidth = 0.8) +
    geom_point(aes(fill = control), shape = 21, size = 2.6, colour = "white", stroke = 0.4) +
    geom_text(aes(label = genes, colour = control, hjust = side, x = genes + (0.5 - side) * 14),
              size = TXT, show.legend = FALSE) +
    scale_fill_manual(values = setNames(COL_PAIR, levels(FPl$control)), name = NULL) +
    scale_colour_manual(values = setNames(COL_PAIR, levels(FPl$control))) +
    scale_y_discrete(limits = rev(METHODS)) +
    scale_x_continuous(limits = c(0, 300), expand = expansion(mult = c(0.02, 0.04))) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = "Spurious significant genes", y = NULL) +
    theme_map(grid = "x")

  # ---- D genes per genome vs N50
  q <- fread(f("02_qc/genome_qc_master.tsv"))
  gcol <- if ("genome" %in% names(q)) "genome" else names(q)[1]
  ncol_ <- grep("^n50$", names(q), ignore.case = TRUE, value = TRUE)[1]
  pD <- NULL
  if (!is.na(ncol_)) {
    n50 <- setNames(as.numeric(q[[ncol_]]), q[[gcol]])
    D <- rbindlist(lapply(c("Annotation", "Reads, three states"), function(m) {
      if (!m %in% names(mats)) return(NULL)
      data.table(genome = ids, genes = colSums(mats[[m]] > 0.5, na.rm = TRUE), method = m)
    }))
    D[, N50 := n50[genome]]
    D <- D[!is.na(N50) & N50 > 0]
    D[, method := factor(method, METHODS)]
    rho <- D[, .(rho = cor(genes, N50, method = "spearman")), by = method]
    rho[, lab := sprintf("%s: \u03c1 = %.2f", method, rho)]
    message("  Spearman genes vs N50: ", paste(rho$lab, collapse = "; "))
    pD <- ggplot(D, aes(N50 / 1000, genes, colour = method)) +
      geom_point(size = 0.45, alpha = 0.45, stroke = 0) +
      geom_smooth(method = "loess", se = FALSE, linewidth = 0.7, formula = y ~ x, span = 0.8) +
      scale_x_log10(breaks = c(5, 10, 20, 50, 100, 200, 500), labels = label_comma()) +
      scale_y_continuous(labels = label_comma()) +
      scale_colour_manual(values = COL_METHOD, labels = setNames(rho$lab, rho$method), name = NULL) +
      guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.6, alpha = 1, linewidth = 0))) +
      labs(x = "Assembly N50 (kb, log scale)", y = "Genes present per genome")
  } else message("  [!] no N50 column in genome_qc_master.tsv; panel D skipped")

  top <- pA + pB + plot_layout(widths = c(1.35, 1))
  bot <- if (is.null(pD)) pC else (pC | pD) + plot_layout(widths = c(1, 1.15))
  p <- tags((top / bot) + plot_layout(heights = c(1.15, 1)))
  save_fig(p, "Fig2", h_cm = 15.5)
}

# =============================================================================
# 5b. FIGURE 3  -  lineage-specific gene content (written by 56_lineage_gene_function.py)
# =============================================================================
CONTRAST_LAB <- c(S_only = "Present only in sheep-type (S)", S_absent = "Present only in cattle-type (C)")
COL_CONTRAST <- setNames(c(REF[["wine"]], REF[["navy"]]), CONTRAST_LAB)

fasta_len <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  x <- readLines(path, warn = FALSE)
  sum(nchar(x[!startsWith(x, ">")]))
}

fig3_lineage <- function() {
  dir <- f("07_functional/lineage_genes")
  if (!exists_or(file.path(dir, "blocks.tsv")) || !exists_or(file.path(dir, "lineage_specific_genes.tsv"))) return()
  B <- fread(file.path(dir, "blocks.tsv"))
  G <- fread(file.path(dir, "lineage_specific_genes.tsv"))
  B <- B[contrast %in% names(CONTRAST_LAB) & !is.na(start)]
  B[, lab := CONTRAST_LAB[contrast]]
  B[, lab := factor(lab, CONTRAST_LAB)]
  B[, refname := fifelse(ref == "Telford", "Telford (S-type reference)", "K-10 (C-type reference)")]
  B[, refname := factor(refname, c("Telford (S-type reference)", "K-10 (C-type reference)"))]
  B[, mid := (start + end) / 2 / 1e6]

  # chromosome backbone: real length if the reference is available, else the data range
  len <- c("Telford (S-type reference)" = fasta_len(f("05_phylogeny/typing/external/Telford.clean.fasta")),
           "K-10 (C-type reference)" = fasta_len(f("03_genomes/fasta/GCA_000007865.1.fasta")))
  BB <- data.table(refname = factor(names(len), levels(B$refname)), len = unname(len) / 1e6)
  BB[is.na(len), len := B[, max(end) / 1e6 * 1.02]]
  message("  reference lengths (Mb): ", paste(BB$refname, round(BB$len, 2), collapse = "; "))

  B[, blab := as.character(n_genes)]
  pA <- ggplot(B) +
    geom_segment(data = BB, aes(x = 0, xend = len, y = 0, yend = 0), inherit.aes = FALSE,
                 linewidth = 1.6, colour = NEUTRAL[["light"]], lineend = "round") +
    geom_segment(aes(x = mid, xend = mid, y = 0, yend = n_genes, colour = lab), linewidth = 0.4,
                 show.legend = FALSE) +
    geom_point(aes(mid, n_genes, fill = lab), shape = 21, size = 2.4, colour = "white", stroke = 0.35,
               show.legend = FALSE) +
    geom_text_repel(aes(mid, n_genes, label = blab, colour = lab), size = TXT * 0.85, family = FONT,
                    lineheight = 0.85, direction = "both", nudge_y = 3, segment.size = 0.2,
                    min.segment.length = 0.4, box.padding = 0.18, max.overlaps = 40, seed = 1,
                    show.legend = FALSE) +
    facet_wrap(~refname, ncol = 1) +
    scale_colour_manual(values = COL_CONTRAST, name = NULL, drop = FALSE) +
    scale_fill_manual(values = COL_CONTRAST, name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = breaks_pretty(6), expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
    guides(colour = "none", fill = "none") +
    labs(x = "Position on the reference chromosome (Mb)", y = "Genes in the block") +
    theme_map(grid = "y")

  # ---- B COG categories
  G <- G[contrast %in% names(CONTRAST_LAB)]
  C <- G[, {k <- strsplit(fifelse(is.na(eggnog_COG_category) | eggnog_COG_category == "", "-",
                                  eggnog_COG_category), "")[[1]]
            .(COG = k)}, by = .(contrast, gene, block)]
  CS <- C[, .(genes = uniqueN(gene), blocks = uniqueN(block)), by = .(contrast, COG)]
  CS[, lab := factor(CONTRAST_LAB[contrast], CONTRAST_LAB)]
  CS[, name := fifelse(COG %in% names(COG_NAME), COG_NAME[COG], COG)]
  CS[, name := sprintf("%s  %s", COG, name)]
  ord <- CS[, .(n = sum(genes)), by = name][order(n)]
  CS[, name := factor(name, ord$name)]
  message("  COG categories drawn: ", nrow(CS))
  pB <- ggplot(CS, aes(genes, name, fill = lab)) +
    geom_col(position = position_dodge2(preserve = "single", padding = 0.1), width = 0.75,
             colour = NEUTRAL[["ink"]], linewidth = 0.15) +
    geom_text(aes(label = sprintf("%d (%d)", genes, blocks)),
              position = position_dodge2(width = 0.75, preserve = "single", padding = 0.1),
              hjust = -0.15, size = TXT * 0.85) +
    scale_fill_manual(values = COL_CONTRAST, name = NULL, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    guides(fill = guide_legend(nrow = 1, keywidth = unit(3, "mm"), keyheight = unit(3, "mm"))) +
    labs(x = "Genes (blocks)", y = NULL) +
    theme_map(grid = "x")

  p <- tags((pA / pB) + plot_layout(heights = c(1, 1.4)) &
              theme(legend.position = "bottom", legend.justification = "center",
                    legend.box.margin = margin(t = 1, b = 0)))
  save_fig(p, "Fig3", h_cm = 17)
  message("  blocks: ", paste(B[, .N, by = lab][, paste0(lab, " ", N)], collapse = "; "))
}

# =============================================================================
# 5. FIGURE 3  -  population structure (circular cladogram)
# =============================================================================
circ_layout <- function(phy, gap = 0.06) {
  phy <- ladderize(phy, right = FALSE)
  n <- Ntip(phy); N <- n + phy$Nnode
  ord <- phy$edge[phy$edge[, 2] <= n, 2]          # plotting order of tips
  x <- numeric(N); x[ord] <- seq_len(n)
  h <- numeric(N); xmin <- rep(Inf, N); xmax <- rep(-Inf, N)
  e <- reorder(phy, "postorder")$edge
  for (i in seq_len(nrow(e))) {
    p <- e[i, 1]; c <- e[i, 2]
    if (c > n) x[c] <- (xmin[c] + xmax[c]) / 2
    h[p] <- max(h[p], h[c] + 1)
    xmin[p] <- min(xmin[p], x[c]); xmax[p] <- max(xmax[p], x[c])
  }
  rt <- n + 1; x[rt] <- (xmin[rt] + xmax[rt]) / 2
  r <- (1 - h / max(h))^0.7                        # cladogram, tips on the unit circle
  # angle: clockwise from 12 o'clock, with a gap of `gap` (fraction of the circle) at the top
  th <- function(x) pi / 2 - 2 * pi * (gap / 2 + (1 - gap) * (x - 0.5) / n)
  edges <- data.table(x = r[e[, 1]] * cos(th(x[e[, 2]])), y = r[e[, 1]] * sin(th(x[e[, 2]])),
                      xend = r[e[, 2]] * cos(th(x[e[, 2]])), yend = r[e[, 2]] * sin(th(x[e[, 2]])))
  inner <- (n + 1):N
  arcs <- rbindlist(lapply(inner, function(v) {
    k <- max(2L, ceiling((xmax[v] - xmin[v]) / n * 240))
    a <- th(seq(xmin[v], xmax[v], length.out = k))
    data.table(id = v, x = r[v] * cos(a), y = r[v] * sin(a))
  }))
  list(phy = phy, n = n, x = x, r = r, xmin = xmin, xmax = xmax, th = th, edges = edges, arcs = arcs,
       tips = data.table(label = phy$tip.label, x = x[seq_len(n)]))
}
# one wedge per tip between radii r0 and r1
ring_poly <- function(L, x, fill, r0, r1, pad = 0.15) {
  a0 <- L$th(x - 0.5 - pad); a1 <- L$th(x + 0.5 + pad); id <- seq_along(x)
  data.table(id = rep(id, 4), fill = rep(fill, 4),
             x = c(r0 * cos(a0), r1 * cos(a0), r1 * cos(a1), r0 * cos(a1)),
             y = c(r0 * sin(a0), r1 * sin(a0), r1 * sin(a1), r0 * sin(a1)))
}
band_poly <- function(L, x0, x1, r0, r1) {
  a <- L$th(seq(x0, x1, length.out = 200))
  data.table(x = c(r0 * cos(a), rev(r1 * cos(a))), y = c(r0 * sin(a), rev(r1 * sin(a))))
}

fig4_tree <- function(meta) {
  phy <- read.tree(f("05_phylogeny/iqtree/tree_final.nwk"))
  s_tips <- intersect(phy$tip.label, meta[lineage == "S", genome])
  phy <- tryCatch(root(phy, outgroup = s_tips, resolve.root = TRUE),
                  error = function(e) { message("  [i] S not monophyletic; midpoint root"); phytools_mid(phy) })
  L <- circ_layout(phy)
  n <- L$n
  T <- merge(L$tips, meta, by.x = "label", by.y = "genome", all.x = TRUE)
  T[is.na(group), group := "C basal grade"]
  T[is.na(region), region := "Other/unknown"]; T[is.na(host), host := "Other/unknown"]
  RW <- 0.065; r1 <- 1.035; r2 <- r1 + RW + 0.012; r3 <- r2 + RW + 0.012
  rings <- rbind(ring_poly(L, T$x, COL_GROUP[as.character(T$group)], r1, r1 + RW)[, ring := 1],
                 ring_poly(L, T$x, COL_REGION[as.character(T$region)], r2, r2 + RW)[, ring := 2],
                 ring_poly(L, T$x, COL_HOST[as.character(T$host)], r3, r3 + RW)[, ring := 3])
  rings[, gid := paste(ring, id)]
  # clade D: shaded sector from its MRCA to the tips
  mD <- getMRCA(L$phy, T[group == "C clade D", label])
  shade <- band_poly(L, L$xmin[mD] - 0.5, L$xmax[mD] + 0.5, L$r[mD], 1.0)
  # group labels outside the rings
  T[, lgroup := fifelse(startsWith(as.character(group), "S/"), "S lineage", as.character(group))]
  LB <- T[, .(x = mean(range(x)), n = .N, col = COL_GROUP[as.character(group[1])]), by = lgroup]
  LB[, lab := sprintf("%s\n(n = %s)", sub("C clade D", "Clade D", lgroup), lab_n(n))]
  setorder(LB, x)
  sep <- n * 0.04                                  # keep neighbouring labels apart
  for (i in seq_len(nrow(LB))[-1]) LB$x[i] <- max(LB$x[i], LB$x[i - 1] + sep)
  rl <- r3 + RW + 0.06
  LB[, `:=`(px = rl * cos(L$th(x)), py = rl * sin(L$th(x)))]
  LB[, hj := fifelse(px > 0.15, 0, fifelse(px < -0.15, 1, 0.5))]
  LB[, vj := fifelse(py > 0.4, 0, fifelse(py < -0.4, 1, 0.5))]
  ring_lab <- data.table(y = c(r1, r2, r3) + RW / 2, lab = c("Lineage", "Region", "Host"))
  ink <- "#3A3A3A"
  pA <- ggplot() +
    geom_polygon(data = shade, aes(x, y), fill = REF[["navy"]], alpha = 0.07) +
    geom_segment(data = L$edges, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.12, colour = ink,
                 lineend = "round") +
    geom_path(data = L$arcs, aes(x, y, group = id), linewidth = 0.12, colour = ink) +
    geom_polygon(data = rings, aes(x, y, group = gid, fill = I(fill)), colour = NA) +
    geom_text(data = ring_lab, aes(x = 0, y = y, label = lab), size = TXT * 0.85, fontface = "bold",
              colour = NEUTRAL[["ink"]]) +
    geom_text(data = LB, aes(px, py, label = lab, hjust = hj, vjust = vj, colour = I(col)), size = TXT, lineheight = 0.85,
              fontface = "bold") +
    coord_equal(xlim = c(-1.42, 1.42), ylim = c(-1.32, 1.32), expand = FALSE, clip = "off") +
    theme_void_map()

  leg <- legend_block(list(Lineage = COL_GROUP, Region = COL_REGION, Host = COL_HOST), ncol = 1)

  # ---- B composition
  B <- meta[, .N, by = .(group, region)]
  B[, frac := N / sum(N), by = group]
  tot <- meta[, .N, by = group]
  pB <- ggplot(B, aes(frac, group, fill = region)) +
    geom_col(width = 0.68, colour = "white", linewidth = 0.25, position = position_fill(reverse = TRUE)) +
    geom_text(data = tot, aes(x = 1.02, y = group, label = lab_n(N)), inherit.aes = FALSE, hjust = 0, size = TXT) +
    scale_fill_manual(values = COL_REGION, name = "Region", drop = FALSE) +
    scale_x_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0.2)), breaks = c(0, .5, 1)) +
    scale_y_discrete(limits = rev(GROUPS)) +
    labs(x = "Genomes", y = NULL) +
    guides(fill = guide_legend(nrow = 2)) +
    theme_map(grid = "") + theme(legend.position = "none") +
    coord_cartesian(clip = "off")

  # ---- C typing scheme
  K <- data.table(level = rep(c("Level 1: S lineage", "Level 2: S type", "Level 3: Bison", "Level 4: clade D"),
                              c(4, 4, 2, 2)),
                  pos = c(9124, 5365, 5751, 6610, 8170, 7037, 9288, 6883, 7470, 8955, 405941, 444953),
                  gene = c("gyrA", "gyrB", "gyrB A165V", "gyrB", "gyrA", "gyrB E594K", "gyrA", "gyrB",
                           "gyrA", "gyrA", "rsfB M28L", "S195G"),
                  allele = c("A", "C", "T", "C", "A", "A", "T", "T", "G", "A", "A", "T"),
                  marks = c(rep("S", 4), "Type I", "Type I", "Type III", "Type III", "Bison", "Bison", "D", "D"))
  pat <- paste0("$2==", K$pos, collapse = " || ")
  ct <- fread(cmd = sprintf("awk -F'\\t' 'NR==1 || %s' %s", pat, shQuote(f("05_phylogeny/core_masked/core.tab"))),
              colClasses = "character")
  gcols <- intersect(names(ct), meta$genome)
  H <- rbindlist(lapply(seq_len(nrow(K)), function(i) {
    row <- ct[as.integer(POS) == K$pos[i]]
    if (!nrow(row)) return(data.table(marker = i, group = GROUPS, frac = NA_real_, n = 0L))
    al <- toupper(unlist(row[1, ..gcols]))
    d <- data.table(genome = gcols, al = al)[meta, on = "genome", nomatch = 0]
    d[!al %in% c("N", "-", "?", ""), .(frac = mean(al == K$allele[i]), n = .N), by = group][, marker := i]
  }), fill = TRUE)
  H <- merge(CJ(marker = seq_len(nrow(K)), group = factor(GROUPS, GROUPS)), H, by = c("marker", "group"), all.x = TRUE)
  H <- cbind(H, K[H$marker, .(level, pos, gene, allele)])
  H[, mlab := sprintf("%s\n%s", gene, lab_n(pos))]
  H[, mlab := factor(mlab, unique(sprintf("%s\n%s", K$gene, lab_n(K$pos))))]
  H[, level := factor(level, unique(K$level))]
  pC <- ggplot(H, aes(mlab, group, fill = frac)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(is.na(frac), "\u2013", ifelse(frac > 0 & frac < 0.995 & frac > 0.005,
                                                               sprintf("%.0f%%", 100 * frac),
                                                               ifelse(frac >= 0.995, "100%", ""))),
                  colour = !is.na(frac) & frac > 0.5), size = TXT) +
    facet_grid(. ~ level, scales = "free_x", space = "free_x") +
    scale_fill_gradientn(colours = c("#F4F6F9", REF[["steel"]], REF[["navy"]]), limits = c(0, 1),
                         labels = label_percent(), name = "Genomes carrying the marker allele",
                         na.value = NEUTRAL[["faint"]]) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = NEUTRAL[["ink"]]), guide = "none") +
    scale_y_discrete(limits = rev(GROUPS)) +
    guides(fill = guide_colourbar(barwidth = unit(35, "mm"), barheight = unit(2.2, "mm"),
                                  title.position = "left", title.vjust = 1)) +
    labs(x = NULL, y = NULL) +
    theme_map(grid = "") +
    theme(panel.border = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(size = BASE - 1.5, lineheight = 0.9, face = "italic"),
          panel.spacing.x = unit(1.5, "mm"),
          axis.text.y = element_text(face = "bold"))

  right <- (leg / (pB + labs(tag = "B"))) + plot_layout(heights = c(1.3, 1))
  top <- ((pA + labs(tag = "A")) | right) + plot_layout(widths = c(2.05, 1))
  p <- (top / (pC + labs(tag = "C"))) + plot_layout(heights = c(2.5, 1)) +
    plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))
  save_fig(p, "Fig4", h_cm = 18)
  message("  group x region:\n", paste(capture.output(print(dcast(B, group ~ region, value.var = "N", fill = 0))),
                                        collapse = "\n"))
}
phytools_mid <- function(phy) { if (!requireNamespace("phangorn", quietly = TRUE)) return(phy); phangorn::midpoint(phy) }

# =============================================================================
# 6. FIGURE 4  -  homoplasy and convergence
# =============================================================================
fig5_homoplasy <- function() {
  A <- fread(f("06_pangwas/convergence_v2/site_events.tsv"))
  binf <- function(v) factor(ifelse(v >= 6, "\u22656", as.character(v)), c(1:5, "\u22656"))
  H <- rbind(A[carriers_raw > 0 & fitch_raw > 0, .(n = .N), by = .(b = binf(fitch_raw))][, calls := "Assembly calls"],
             A[carriers > 0 & fitch > 0, .(n = .N), by = .(b = binf(fitch))][, calls := "Read-verified calls"])
  H[, calls := factor(calls, c("Assembly calls", "Read-verified calls"))]
  ci <- A[carriers > 0, .N / sum(pmax(fitch, 1))]
  one <- A[carriers >= 2, mean(fitch == 1)]
  pA <- ggplot(H, aes(b, n, fill = calls)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.74, colour = NEUTRAL[["ink"]], linewidth = 0.2) +
    geom_text(aes(label = lab_n(n)), position = position_dodge(width = 0.78), vjust = -0.4, size = TXT * 0.9) +
    annotate("text", x = 6.45, y = max(H$n) * 2.2, hjust = 1, vjust = 1, size = TXT, lineheight = 0.95,
             label = sprintf("Consistency index %.3f\n%.1f%% of alleles: one change", ci, 100 * one)) +
    scale_y_log10(labels = label_comma(), expand = expansion(mult = c(0, 0.25))) +
    scale_fill_manual(values = setNames(COL_PAIR, levels(H$calls)), name = NULL) +
    labs(x = "Changes on the tree per SNP allele", y = "SNP alleles (log scale)") +
    theme_map(grid = "y")

  # ---- B 21 flagged SNPs (45_discordant_tips)
  EX <- data.table(why = c("Mixed-strain samples", "Single terminal change", "Assembly without reads",
                           "Candidate parallel mutation"), n = c(16, 3, 1, 1))
  EX[, why := factor(why, why)]
  W <- EX[rep(seq_len(.N), n)][, i := .I - 1][, `:=`(cx = i %/% 7, cy = -(i %% 7))]
  COL_WHY <- setNames(c(REF[["rose"]], REF[["steel"]], REF[["tan"]], REF[["wine"]]), levels(EX$why))
  pB <- ggplot(W, aes(cx, cy, fill = why)) +
    geom_tile(width = 0.86, height = 0.86, colour = NA) +
    coord_equal() +
    scale_fill_manual(values = COL_WHY, name = NULL,
                      labels = sprintf("%s (%d)", levels(EX$why), EX$n)) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    labs(subtitle = "21 SNPs first flagged as homoplastic") +
    theme_void_map() +
    theme(legend.position = "right", legend.text = element_text(size = BASE - 1),
          legend.key.size = unit(2.8, "mm"),
          plot.subtitle = element_text(size = BASE - 0.5, face = "bold", hjust = 0.5, margin = margin(b = 3)))

  # ---- C gene-level convergence along the chromosome
  G <- fread(f("06_pangwas/convergence_v2/gene_branch_convergence.tsv"))
  G[, name := fifelse(!is.na(gene) & gene != "", gene, locus_tag)]
  G[, sig := q < 0.05]
  thr <- if (any(G$sig)) max(G[sig == TRUE, p]) else NA
  lab <- G[sig == TRUE][order(p)][seq_len(min(.N, 7))]
  pC <- ggplot(G, aes(start / 1e6, -log10(p))) +
    geom_point(data = G[sig == FALSE], aes(size = branch_events), shape = 16, colour = REF[["steel"]], alpha = 0.45) +
    { if (!is.na(thr)) geom_hline(yintercept = -log10(thr), linetype = "22", colour = REF[["wine"]], linewidth = 0.4) } +
    geom_point(data = G[sig == TRUE], aes(size = branch_events), shape = 21, fill = REF[["wine"]], colour = "white",
               stroke = 0.3) +
    geom_text_repel(data = lab, aes(label = name), size = TXT, fontface = "italic", family = FONT,
                    min.segment.length = 0,
                    segment.size = 0.2, segment.colour = "grey40", box.padding = 0.3, max.overlaps = 30, seed = 1) +
    { if (!is.na(thr)) annotate("text", x = 0, y = -log10(thr), label = "FDR 5%", vjust = -0.5,
                                hjust = 0, size = TXT * 0.9, colour = REF[["wine"]], fontface = "italic") } +
    scale_size_area(max_size = 3.2, breaks = c(1, 5, 10, 15), name = "Branch events") +
    scale_x_continuous(breaks = 0:5, expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = "Position on the K-10 chromosome (Mb)", y = expression(bold(-log[10]~italic(P)))) +
    theme_map(grid = "y")

  # ---- D icl lollipop
  pD <- NULL
  ef <- f("06_pangwas/convergence_v2/gene_events_detail.tsv")
  if (exists_or(ef)) {
    E <- fread(ef)
    E <- E[gene == "icl" | locus_tag == "BDOFFIBO_04182"]
    if (nrow(E)) {
      len_aa <- G[name == "icl" | locus_tag == "BDOFFIBO_04182", length][1] / 3 - 1
      E[, aa_pos := as.integer(gsub("^[A-Z*]+([0-9]+)[A-Z*]+$", "\\1", aa))]
      E[, lin := fifelse(startsWith(lineages, "S/"), "S lineage",
                         fifelse(startsWith(lineages, "C/B"), "Bison", "C lineage (non-Bison)"))]
      E[, lin := factor(lin, c("S lineage", "Bison", "C lineage (non-Bison)"))]
      setorder(E, aa_pos, lin)
      E[, k := seq_len(.N), by = aa]
      top <- E[, .(k = max(k), aa_pos = aa_pos[1], effect = effect[1]), by = aa]
      COL_LIN <- c(`S lineage` = REF[["wine"]], Bison = REF[["olive"]], `C lineage (non-Bison)` = REF[["navy"]])
      pD <- ggplot() +
        annotate("rect", xmin = 1, xmax = len_aa, ymin = -0.55, ymax = 0, fill = REF[["sand"]],
                 colour = REF[["tan"]], linewidth = 0.3) +
        annotate("text", x = len_aa / 2, y = -0.28, label = sprintf("Icl (isocitrate lyase), %d aa", round(len_aa)),
                 size = TXT * 0.9, fontface = "italic") +
        geom_segment(data = top, aes(x = aa_pos, xend = aa_pos, y = 0, yend = k), colour = "grey45", linewidth = 0.35) +
        geom_point(data = E, aes(aa_pos, k, fill = lin), shape = 21, size = 2.2, colour = "white", stroke = 0.35) +
        geom_text_repel(data = top, aes(aa_pos, k, label = aa), size = TXT * 0.95, fontface = "bold", family = FONT,
                        nudge_y = 1.1, direction = "x", segment.size = 0.2, segment.colour = "grey55",
                        min.segment.length = 0.3, box.padding = 0.15, point.padding = 0.2, seed = 1) +
        scale_fill_manual(values = COL_LIN, name = NULL, drop = FALSE) +
        scale_x_continuous(limits = c(-25, len_aa + 15), expand = c(0, 0)) +
        scale_y_continuous(breaks = seq(0, 14, 2), limits = c(-0.6, max(top$k) + 2.2), expand = c(0, 0)) +
        labs(x = "Amino-acid position", y = "Independent changes") +
        theme_map(grid = "y")
    }
  }
  top <- (pA | pB) + plot_layout(widths = c(1.3, 1))
  bot <- if (is.null(pD)) pC else (pC | pD) + plot_layout(widths = c(1.3, 1))
  p <- tags((top / bot) + plot_layout(heights = c(1, 1)))
  save_fig(p, "Fig5", h_cm = 14.5)
  message(sprintf("  verified CI %.3f; alleles >= 2 carriers with one change %.1f%%; genes q<0.05: %d",
                  ci, 100 * one, sum(G$sig)))
}

# =============================================================================
# 7. FIGURE 5  -  sampling structure
# =============================================================================
fig6_confounding <- function() {
  R <- fread(f("09_dating/fig/rtt_C_nonB.tsv"))
  R[, region := factor(fifelse(region %in% REGIONS[1:3], region, "Other/unknown"), REGIONS)]
  fit <- lm(rtt ~ date, R); r2 <- summary(fit)$r.squared
  pA <- ggplot(R, aes(date, rtt)) +
    geom_point(aes(fill = region), shape = 21, colour = "white", stroke = 0.15, size = 1.3, alpha = 0.9) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", linetype = "22", linewidth = 0.6,
                fill = "grey60", alpha = 0.25) +
    annotate("text", x = min(R$date), y = max(R$rtt), hjust = 0, vjust = 1, size = TXT,
             label = sprintf("slope = %.2f per year\nR\u00b2 = %.3f; n = %s", coef(fit)[2], r2, lab_n(nrow(R)))) +
    scale_fill_manual(values = COL_REGION, name = "Region", drop = TRUE) +
    guides(fill = guide_legend(override.aes = list(size = 2.2))) +
    labs(x = "Sampling date", y = "Root-to-tip distance")

  P <- fread(f("09_dating/fig/perm_C_nonB.tsv"))
  obs <- P[scheme == "observed", rate][1]
  SCH <- c(free = "All dates permuted", within_region = "Within regions", within_study = "Within studies")
  Pp <- P[scheme %in% names(SCH)][, scheme := factor(SCH[scheme], rev(SCH))]
  pv <- Pp[, .(p = (sum(rate >= obs) + 1) / (.N + 1)), by = scheme]
  pv[, lab := ifelse(p < 0.001, "P < 0.001", sprintf("P = %.2f", p))]
  pB <- ggplot(Pp, aes(rate, scheme, fill = scheme)) +
    geom_density_ridges(scale = 1.6, colour = "white", linewidth = 0.3, alpha = 0.95, rel_min_height = 0.005) +
    geom_vline(xintercept = obs, colour = REF[["wine"]], linetype = "22", linewidth = 0.5) +
    annotate("text", x = obs, y = 3.95, label = "observed", colour = REF[["wine"]], hjust = -0.08, size = TXT,
             fontface = "italic") +
    geom_text(data = pv, aes(x = Inf, y = scheme, label = lab), hjust = 1.1, vjust = -0.8, size = TXT,
              fontface = "italic", inherit.aes = FALSE) +
    scale_fill_manual(values = setNames(c(REF[["steel"]], REF[["tan"]], REF[["slate"]]), SCH), guide = "none") +
    scale_y_discrete(expand = expansion(add = c(0.25, 1.1))) +
    labs(x = "Root-to-tip slope (per year)", y = NULL) +
    theme_map(grid = "x")

  B <- fread(f("09_dating/fig/rtt_by_study.tsv"))
  S <- B[, {m <- lm(rtt ~ date); ci <- confint(m)[2, ]
            .(slope = coef(m)[2], lo = ci[1], hi = ci[2], r2 = summary(m)$r.squared, n = .N,
              span = sprintf("%d\u2013%d", floor(min(date)), floor(max(date))))}, by = study]
  S[, study := factor(study, S[order(slope), study])]
  logx <- all(S$lo > 0)
  pC <- ggplot(S, aes(slope, study)) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, colour = "grey40", linewidth = 0.4) +
    geom_point(aes(size = n), shape = 21, fill = REF[["navy"]], colour = "white", stroke = 0.4) +
    geom_text(aes(x = hi, label = sprintf("n = %s, %s", lab_n(n), span)), hjust = -0.12, size = TXT * 0.9) +
    { if (logx) scale_x_log10(expand = expansion(mult = c(0.08, 0.55))) else
      scale_x_continuous(expand = expansion(mult = c(0.08, 0.6))) } +
    scale_size_area(max_size = 3.4, guide = "none") +
    labs(x = if (logx) "Within-study slope (per year, log scale)" else "Within-study slope (per year)", y = NULL) +
    theme_map(grid = "x")

  M <- fread(f("07_ml/C_nonB/metrics.tsv"))
  message("  ML labels: ", paste(unique(M$label), collapse = ", "), " | schemes: ", paste(unique(M$scheme), collapse = ", "))
  MOD <- c(tree_1NN = "Nearest relative (tree)", snp_1NN = "Nearest relative (SNPs)",
           logreg_L1 = "L1 logistic regression", rf = "Random forest")
  SC <- c(random_cv = "Random cross-validation", leave_study_out = "Leave one study out")
  D <- M[label == "region" & scheme %in% names(SC) & model %in% names(MOD)]
  D[, `:=`(model = factor(MOD[model], rev(MOD)), scheme = factor(SC[scheme], SC))]
  pD <- ggplot(D, aes(balanced_accuracy, model)) +
    geom_vline(xintercept = 1 / 3, linetype = "22", colour = "grey45", linewidth = 0.4) +
    annotate("text", x = 1 / 3, y = 4.45, label = "chance", hjust = 1.1, size = TXT * 0.9, colour = "grey35",
             fontface = "italic") +
    geom_line(aes(group = model), colour = "grey70", linewidth = 0.9) +
    geom_point(aes(fill = scheme), shape = 21, size = 2.6, colour = "white", stroke = 0.4) +
    scale_fill_manual(values = setNames(COL_PAIR, SC), name = NULL) +
    scale_x_continuous(limits = c(0.25, 1), breaks = seq(0.3, 1, 0.1)) +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = "Balanced accuracy (region)", y = NULL) +
    theme_map(grid = "x")

  p <- tags(((pA | pB) / (pC | pD)) + plot_layout(heights = c(1.1, 1)))
  save_fig(p, "Fig6", h_cm = 14.5)
  message(sprintf("  rtt: n %d, slope %.3f, R2 %.4f; observed perm slope %.3f", nrow(R), coef(fit)[2], r2, obs))
  print(pv)
}

# =============================================================================
# 8. SUPPLEMENTAL FIGURES
# =============================================================================
figS2 <- function() {
  # 21_batch_check: median |frequency difference| of hit genes; 22_contig_edge
  BATCH <- data.table(trait = c("Americas", "Europe", "Oceania", "Environment", "Excretory", "Invasive", "Cattle"),
                      `All C-lineage genomes` = c(0.07, 0.09, 0.14, 0.38, 0.11, 0.13, 0.23),
                      `Within one platform` = c(0.01, 0.00, 0.09, 0.04, 0.05, 0.07, 0.04))
  Bl <- melt(BATCH, id.vars = "trait", variable.name = "set", value.name = "d")
  Bl[, trait := factor(trait, BATCH[order(`All C-lineage genomes`), trait])]
  pA <- ggplot(Bl, aes(d, trait)) +
    geom_line(aes(group = trait), colour = "grey70", linewidth = 0.9) +
    geom_point(aes(fill = set), shape = 21, size = 2.6, colour = "white", stroke = 0.4) +
    scale_fill_manual(values = setNames(COL_PAIR, levels(Bl$set)), name = NULL) +
    scale_x_continuous(limits = c(0, 0.4)) +
    labs(x = "Median absolute difference in gene frequency", y = NULL) +
    theme_map(grid = "x")
  EDGE <- data.table(set = c("Core genes", "Accessory genes", "Hits, non-robust (enriched)",
                             "Hits, non-robust (depleted)", "Hits, robust"),
                     pct = c(2.6, 8.2, 15.8, 16.0, 56.9),
                     kind = c("Reference", "Reference", "Associated", "Associated", "Associated"))
  EDGE[, set := factor(set, rev(set))]
  pB <- ggplot(EDGE, aes(pct, set, fill = kind)) +
    geom_col(width = 0.66, colour = NEUTRAL[["ink"]], linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.15, size = TXT) +
    scale_fill_manual(values = c(Reference = REF[["steel"]], Associated = REF[["wine"]]), name = NULL,
                      breaks = c("Reference", "Associated"),
                      labels = c("All genes (reference)", "Genes associated with a trait")) +
    scale_x_continuous(limits = c(0, 66), expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(x = "Gene copies within 100 bp of a contig end", y = NULL) +
    theme_map(grid = "x")
  p <- tags(pA | pB)
  save_fig(p, "FigS2", h_cm = 7.5)
}

figS3 <- function() {
  B <- fread(f("09_dating/fig/rtt_by_study.tsv"))
  S <- B[, {m <- lm(rtt ~ date); .(slope = coef(m)[2], r2 = summary(m)$r.squared, n = .N,
                                    x = min(date), y = max(rtt))}, by = study]
  S[, lab := sprintf("slope %.2f per year\nR\u00b2 = %.3f; n = %d", slope, r2, n)]
  p <- ggplot(B, aes(date, rtt)) +
    geom_point(shape = 21, fill = REF[["steel"]], colour = "white", stroke = 0.15, size = 1.3) +
    geom_smooth(method = "lm", formula = y ~ x, colour = REF[["wine"]], fill = REF[["blush"]], linewidth = 0.6,
                alpha = 0.6) +
    geom_text(data = S, aes(x = x, y = y, label = lab), hjust = 0, vjust = 1, size = TXT, lineheight = 0.95) +
    facet_wrap(~study, scales = "free") +
    scale_x_continuous(breaks = function(l) unique(round(pretty(l)))) +
    labs(x = "Sampling date", y = "Root-to-tip distance")
  save_fig(p, "FigS3", h_cm = 6.5 * ceiling(uniqueN(B$study) / 3))
}

# =============================================================================
# 9. RUN
# =============================================================================
ONLY <- strsplit(Sys.getenv("MAP_FIGS"), ",")[[1]]     # e.g. MAP_FIGS=Fig3,Fig4 to redraw a subset
run <- function(nm, fn) {
  if (length(ONLY) && !nm %in% ONLY) return(invisible())
  message("== ", nm)
  tryCatch(fn(), error = function(e) message("  [!] ", nm, " failed: ", conditionMessage(e)))
}
meta <- read_meta()
run("Fig2", function() fig2(meta))
run("Fig3", fig3_lineage)
run("Fig4", function() fig4_tree(meta))
run("Fig5", fig5_homoplasy)
run("Fig6", fig6_confounding)
run("FigS2", figS2)
run("FigS3", figS3)
message("written: ", OUT)
