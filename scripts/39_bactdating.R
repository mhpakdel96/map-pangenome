#!/usr/bin/env Rscript
# =============================================================================
# 39_bactdating.R
#
# Dated phylogeny of MAP with BactDating (Didelot et al. 2018, NAR), on the
# Gubbins recombination-corrected tree. Two stages, run separately:
#
#   Rscript 39_bactdating.R signal  [subset]
#       Temporal signal: root-to-tip regression with a date-permutation test
#       (BactDating::roottotip, 10000 permutations). If there is no signal
#       (p >= 0.05 or a negative rate), dating is NOT done for that subset.
#   Rscript 39_bactdating.R date    [subset] [nbIts]
#       BactDating MCMC, strict ('poisson') and relaxed ('relaxedgamma')
#       clock, with Gubbins' unrecombined substitutions (useRec = TRUE),
#       compared by DIC (modelcompare); convergence by effective sample size.
#
# subset: all | C_nonB | C | S   (lineages from 16_lineage_typing)
# Tips without a usable date (09_dating/dates.tsv, used == FALSE: no date,
# culture-collection strains, possible mixed infections) are dropped.
# Dates: midpoint of the reported interval (year-only dates -> mid-year);
# if this BactDating version accepts a 2-column date matrix, the interval is
# used instead (reported in the log).
# Outputs: 09_dating/<subset>/...
# =============================================================================
suppressPackageStartupMessages({
  library(ape)
  library(BactDating)
  library(coda)
})
args <- commandArgs(trailingOnly = TRUE)
mode   <- ifelse(length(args) >= 1, args[1], "signal")
subset <- ifelse(length(args) >= 2, args[2], "all")
nbIts  <- ifelse(length(args) >= 3, as.numeric(args[3]), 1e5)
root   <- "/cfs/earth/scratch/tabasran/MAP"
prefix <- file.path(root, "05_phylogeny/gubbins/gubbins")
out    <- file.path(root, "09_dating", subset)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(1)
cat("BactDating", as.character(packageVersion("BactDating")), "| ape",
    as.character(packageVersion("ape")), "| mode", mode, "| subset", subset, "\n")

# ---------------------------------------------------------------- tree + dates
tree <- loadGubbins(prefix)
cat("Gubbins tree:", Ntip(tree), "tips; unrecombined substitutions available:",
    !is.null(tree$unrec), "\n")
d   <- read.delim(file.path(root, "09_dating/dates.tsv"), stringsAsFactors = FALSE)
lin <- read.delim(file.path(root, "05_phylogeny/typing/lineage_assignments.tsv"),
                  stringsAsFactors = FALSE)
lin$lab <- paste(lin$lineage, lin$sublineage, sep = "/")
keep <- d$genome[d$used %in% c(TRUE, "True", "TRUE")]
if (subset == "C_nonB") keep <- intersect(keep, lin$genome[lin$lab == "C/non-B"])
if (subset == "C")      keep <- intersect(keep, lin$genome[lin$lineage == "C"])
if (subset == "S")      keep <- intersect(keep, lin$genome[lin$lineage == "S"])
keep <- intersect(keep, tree$tip.label)
drop <- setdiff(tree$tip.label, keep)
cat("tips kept:", length(keep), "| dropped:", length(drop), "\n")
if (length(drop)) {
  # ape::drop.tip re-indexes edges, so Gubbins' per-branch unrecombined
  # fractions ($unrec) cannot be carried over; the Gubbins tree is already
  # built from non-recombinant sites, so dating proceeds with useRec = FALSE
  tree <- ape::drop.tip(tree, drop)
  tree$unrec <- NULL
  cat("[i] tips dropped -> $unrec removed; useRec = FALSE",
      "(branch lengths are Gubbins' recombination-corrected lengths)\n")
}
dd  <- d[match(tree$tip.label, d$genome), ]
mid <- dd$mid
names(mid) <- dd$genome
rng <- cbind(dd$lower, dd$upper)
cat("dates: ", sprintf("%.1f", min(mid)), "-", sprintf("%.1f", max(mid)),
    "| exact-day tips:", sum(dd$precision == "day"), "| year/month/range tips:",
    sum(dd$precision != "day"), "\n")

# ---------------------------------------------------------------- signal
if (mode == "signal") {
  pdf(file.path(out, "roottotip.pdf"), width = 8, height = 6)
  r <- roottotip(tree, mid, permTest = 10000, showFig = TRUE)
  dev.off()
  sink(file.path(out, "roottotip.txt"))
  print(r)
  sink()
  print(r)
  cat("\nInterpretation: rate > 0 and permutation p < 0.05 -> temporal signal;",
      "otherwise do not date this subset.\n")
  quit(status = 0)
}

# ---------------------------------------------------------------- dating
use_rec <- !is.null(tree$unrec)
dates_in <- tryCatch({
  # interval dates, if supported by this version (2-column matrix)
  bactdate(tree, rng, nbIts = 100, model = "poisson", useRec = use_rec, showProgress = FALSE)
  cat("[i] date intervals accepted: using lower/upper bounds\n")
  rng
}, error = function(e) {
  cat("[i] date intervals not accepted (", conditionMessage(e), ") -> using midpoints\n")
  mid
})
fits <- list()
for (m in c("poisson", "relaxedgamma")) {
  t0 <- Sys.time()
  fits[[m]] <- bactdate(tree, dates_in, nbIts = nbIts, model = m, useRec = use_rec,
                        showProgress = FALSE)
  cat("\n=== model", m, "(", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
      "min ) ===\n")
  print(fits[[m]])
  mc <- as.mcmc.resBactDating(fits[[m]], burnin = 0.5)
  cat("effective sample sizes:\n")
  print(round(effectiveSize(mc)))
  saveRDS(fits[[m]], file.path(out, paste0("bactdate_", m, ".rds")))
  pdf(file.path(out, paste0("trace_", m, ".pdf")), width = 9, height = 7)
  try(plot(fits[[m]], "trace"))
  dev.off()
  pdf(file.path(out, paste0("dated_tree_", m, ".pdf")), width = 10, height = 30)
  try(plot(fits[[m]], "treeCI"))
  dev.off()
}
cat("\n=== model comparison (DIC; lower is better) ===\n")
print(modelcompare(fits[["poisson"]], fits[["relaxedgamma"]]))
cat("\nAccept a result only if ESS >= 200 for mu, sigma (relaxed), alpha and the",
    "root date; otherwise rerun with more iterations.\n")
