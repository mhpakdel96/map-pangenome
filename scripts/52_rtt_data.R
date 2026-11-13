#!/usr/bin/env Rscript
# =============================================================================
# 52_rtt_data.R
#
# Data behind the temporal-signal panels (Fig. 5B and Fig. S4), written as
# tables so that 53_figures.py can draw them. Same steps as 40 and 41:
#   C/non-B genomes with a usable date; Gubbins tree pruned to them; root
#   chosen by BactDating::initRoot; root-to-tip distance per tip; date
#   permutations free, within study and within region (1,000 each).
#   Within-study regressions for the studies that passed 41 (>= 40 dated
#   C/non-B genomes spanning >= 8 years).
# Outputs (09_dating/fig/):
#   rtt_C_nonB.tsv        genome, date, rtt, region, study
#   perm_C_nonB.tsv       scheme, rate  (permuted slopes) + observed slope
#   rtt_by_study.tsv      study, genome, date, rtt
#   Rscript 52_rtt_data.R          (bactdating_env)
# =============================================================================
suppressPackageStartupMessages({ library(ape); library(BactDating) })
root <- "/cfs/earth/scratch/tabasran/MAP"
out  <- file.path(root, "09_dating/fig")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(1)
tree0 <- loadGubbins(file.path(root, "05_phylogeny/gubbins/gubbins"))
tree0$unrec <- NULL
d   <- read.delim(file.path(root, "09_dating/dates.tsv"), stringsAsFactors = FALSE)
lin <- read.delim(file.path(root, "05_phylogeny/typing/lineage_assignments.tsv"), stringsAsFactors = FALSE)
ena <- read.delim(file.path(root, "01_metadata/ena_raw.tsv"), stringsAsFactors = FALSE)
tr  <- read.csv(file.path(root, "06_pangwas/trait.csv"), stringsAsFactors = FALSE)

use <- d[d$used %in% c(TRUE, "True", "TRUE"), ]
use <- use[use$genome %in% lin$genome[lin$lineage == "C" & lin$sublineage == "non-B"], ]
use$study <- ena$study_accession[match(use$genome, ena$run_accession)]
use$study[is.na(use$study)] <- "assembly"
use$region <- "unknown"
for (c in c("Region_Americas", "Region_Europe", "Region_Oceania")) {
  use$region[tr[[c]][match(use$genome, tr$Name)] %in% 1] <- sub("Region_", "", c)
}

rtt_of <- function(g) {
  t <- drop.tip(tree0, setdiff(tree0$tip.label, g))
  x <- use$mid[match(t$tip.label, use$genome)]
  names(x) <- t$tip.label
  rt <- initRoot(t, x)
  r <- node.depth.edgelength(rt)[1:Ntip(rt)]
  data.frame(genome = rt$tip.label, date = x[rt$tip.label], rtt = r, stringsAsFactors = FALSE)
}

# ---- all C/non-B
R <- rtt_of(use$genome)
R$region <- use$region[match(R$genome, use$genome)]
R$study  <- use$study[match(R$genome, use$genome)]
write.table(R, file.path(out, "rtt_C_nonB.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
obs <- coef(lm(rtt ~ date, R))[2]
perm <- function(groups, n = 1000) {
  sapply(1:n, function(k) {
    xp <- R$date
    for (g in unique(groups)) { w <- which(groups == g); if (length(w) > 1) xp[w] <- sample(R$date[w]) }
    coef(lm(R$rtt ~ xp))[2]
  })
}
P <- rbind(data.frame(scheme = "free", rate = perm(rep("all", nrow(R)))),
           data.frame(scheme = "within_study", rate = perm(R$study)),
           data.frame(scheme = "within_region", rate = perm(R$region)),
           data.frame(scheme = "observed", rate = obs))
write.table(P, file.path(out, "perm_C_nonB.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("C/non-B: %d tips, observed rate %.3f, R2 %.4f\n", nrow(R), obs, summary(lm(rtt ~ date, R))$r.squared))

# ---- within studies (same selection as 41)
S <- do.call(rbind, lapply(split(use, use$study), function(x)
  data.frame(study = x$study[1], n = nrow(x), span = max(x$mid) - min(x$mid))))
keep <- S$study[S$n >= 40 & S$span >= 8 & S$study != "assembly"]
B <- do.call(rbind, lapply(keep, function(s) cbind(study = s, rtt_of(use$genome[use$study == s]))))
write.table(B, file.path(out, "rtt_by_study.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("studies:", paste(keep, collapse = ", "), "\nwritten:", out, "\n")
