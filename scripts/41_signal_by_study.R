#!/usr/bin/env Rscript
# =============================================================================
# 41_signal_by_study.R
#
# 40 showed that the C_nonB temporal signal survives date permutation within
# regions but NOT within studies (p = 0.98): the regression is carried by
# differences BETWEEN studies (which sampled different clades in different
# years). Whether MAP has a measurable clock must therefore be tested where
# one study sampled the same population over many years. For every study
# (and the subtree of its genomes) with >= --min-n dated C/non-B genomes
# spanning >= --min-span years:
#   root chosen by initRoot, root-to-tip rate, root date, R^2, and a free
#   date permutation test (1000; meaningful here because all tips come from
#   one study). Rates are substitutions per genome per year (40: total tree
#   length 35,267 vs 27,679 Gubbins sites -> branch lengths are per genome).
#
#   Rscript 41_signal_by_study.R [min_n=40] [min_span=8]
# =============================================================================
suppressPackageStartupMessages({ library(ape); library(BactDating) })
a <- commandArgs(trailingOnly = TRUE)
min_n <- ifelse(length(a) >= 1, as.numeric(a[1]), 40)
min_span <- ifelse(length(a) >= 2, as.numeric(a[2]), 8)
root <- "/cfs/earth/scratch/tabasran/MAP"
set.seed(1)
tree0 <- loadGubbins(file.path(root, "05_phylogeny/gubbins/gubbins"))
tree0$unrec <- NULL
d   <- read.delim(file.path(root, "09_dating/dates.tsv"), stringsAsFactors = FALSE)
lin <- read.delim(file.path(root, "05_phylogeny/typing/lineage_assignments.tsv"), stringsAsFactors = FALSE)
ena <- read.delim(file.path(root, "01_metadata/ena_raw.tsv"), stringsAsFactors = FALSE)
use <- d[d$used %in% c(TRUE, "True", "TRUE"), ]
use <- use[use$genome %in% lin$genome[lin$lineage == "C" & lin$sublineage == "non-B"], ]
use$study <- ena$study_accession[match(use$genome, ena$run_accession)]
use$study[is.na(use$study)] <- "assembly"

cat("per-study date spans (dated C/non-B genomes):\n")
sp <- do.call(rbind, lapply(split(use, use$study), function(x)
  data.frame(study = x$study[1], n = nrow(x), first = min(x$mid), last = max(x$mid),
             span = max(x$mid) - min(x$mid), day_precision = mean(x$precision == "day"))))
sp <- sp[order(-sp$n), ]
print(format(sp, digits = 4), row.names = FALSE)

res <- list()
for (s in sp$study[sp$n >= min_n & sp$span >= min_span & sp$study != "assembly"]) {
  g <- use$genome[use$study == s]
  t <- drop.tip(tree0, setdiff(tree0$tip.label, g))
  x <- use$mid[match(t$tip.label, use$genome)]
  names(x) <- t$tip.label
  rt <- initRoot(t, x)
  rtt <- node.depth.edgelength(rt)[1:Ntip(rt)]
  xx <- x[rt$tip.label]
  f <- lm(rtt ~ xx)
  slope <- coef(f)[2]
  perm <- replicate(1000, coef(lm(rtt ~ sample(xx)))[2])
  res[[s]] <- data.frame(study = s, n = length(g), span = round(max(xx) - min(xx), 1),
                         rate_per_genome_year = round(slope, 3),
                         root_date = round(-coef(f)[1] / slope, 1),
                         R2 = round(summary(f)$r.squared, 3),
                         p_perm = mean(perm >= slope))
  pdf(file.path(root, "09_dating", paste0("rtt_", s, ".pdf")), width = 6, height = 5)
  plot(xx, rtt, pch = 16, cex = .6, xlab = "sampling date", ylab = "root-to-tip",
       main = sprintf("%s  n=%d  R2=%.3f", s, length(g), summary(f)$r.squared))
  abline(f)
  dev.off()
}
cat("\nroot-to-tip within single studies (subtree of that study's genomes):\n")
if (length(res)) print(do.call(rbind, res), row.names = FALSE) else
  cat("  no study meets n >=", min_n, "and span >=", min_span, "years\n")
cat("\nA study supports a measurable clock if rate > 0 with p_perm < 0.05 and a",
    "non-trivial R^2; consistent rates across studies support a common clock.\n")
