#!/usr/bin/env Rscript
# =============================================================================
# 40_temporal_signal.R
#
# Stricter temporal-signal test for the subset that passed 39 (C_nonB).
# The plain date permutation of 39 can be significant only because sampling
# time differs between clades/studies/regions (here: Oceania median 2005,
# Americas 2012, Europe 2016). Clustered permutations shuffle dates only
# WITHIN a group, so a clade-level date difference cannot create the signal
# (Duchene et al. 2015, Mol Biol Evol; Murray et al. 2016, Methods Ecol Evol).
#
#   1. root chosen by BactDating::initRoot (best root-to-tip fit to the dates)
#   2. root-to-tip regression: rate, x-intercept (root date), R^2
#   3. permutations (1000 each): free; within study; within region
#   4. units check: total tree length vs the number of Gubbins SNP sites, to
#      state the rate as substitutions per GENOME per year correctly
#
#   Rscript 40_temporal_signal.R [subset]      (default C_nonB)
# =============================================================================
suppressPackageStartupMessages({ library(ape); library(BactDating) })
args   <- commandArgs(trailingOnly = TRUE)
subset <- ifelse(length(args) >= 1, args[1], "C_nonB")
root   <- "/cfs/earth/scratch/tabasran/MAP"
out    <- file.path(root, "09_dating", subset)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

tree <- loadGubbins(file.path(root, "05_phylogeny/gubbins/gubbins"))
d    <- read.delim(file.path(root, "09_dating/dates.tsv"), stringsAsFactors = FALSE)
lin  <- read.delim(file.path(root, "05_phylogeny/typing/lineage_assignments.tsv"), stringsAsFactors = FALSE)
lin$lab <- paste(lin$lineage, lin$sublineage, sep = "/")
ena  <- read.delim(file.path(root, "01_metadata/ena_raw.tsv"), stringsAsFactors = FALSE)
tr   <- read.csv(file.path(root, "06_pangwas/trait.csv"), stringsAsFactors = FALSE)

keep <- d$genome[d$used %in% c(TRUE, "True", "TRUE")]
if (subset == "C_nonB") keep <- intersect(keep, lin$genome[lin$lab == "C/non-B"])
if (subset == "C")      keep <- intersect(keep, lin$genome[lin$lineage == "C"])
tree <- drop.tip(tree, setdiff(tree$tip.label, keep))
tree$unrec <- NULL
dates <- d$mid[match(tree$tip.label, d$genome)]
names(dates) <- tree$tip.label

# units: Gubbins polymorphic sites vs total branch length
nsites <- NA
fp <- file.path(root, "05_phylogeny/gubbins/gubbins.filtered_polymorphic_sites.fasta")
if (file.exists(fp)) {
  l <- readLines(fp, n = 2000)
  i <- which(startsWith(l, ">"))
  nsites <- sum(nchar(l[(i[1] + 1):(i[2] - 1)]))
}
cat("subset", subset, "| tips", Ntip(tree), "| total tree length", round(sum(tree$edge.length), 2),
    "| Gubbins polymorphic sites", nsites, "\n")
cat("-> if total length is in the thousands it is in substitutions per genome;",
    "if < 1 it is per site and must be multiplied by genome length\n")

rt <- initRoot(tree, dates)
rtt <- node.depth.edgelength(rt)[1:Ntip(rt)]
names(rtt) <- rt$tip.label
x <- dates[rt$tip.label]
fit <- lm(rtt ~ x)
slope <- coef(fit)[2]
r2 <- summary(fit)$r.squared
cat(sprintf("root-to-tip after initRoot: rate %.4f per year, root date %.1f, R^2 %.4f\n",
            slope, -coef(fit)[1] / slope, r2))

study <- ena$study_accession[match(rt$tip.label, ena$run_accession)]
study[is.na(study)] <- "assembly"
region <- rep("unknown", length(rt$tip.label))
for (c in c("Region_Americas", "Region_Europe", "Region_Oceania")) {
  region[tr[[c]][match(rt$tip.label, tr$Name)] %in% 1] <- sub("Region_", "", c)
}
perm <- function(groups, n = 1000) {
  s <- numeric(n)
  for (k in 1:n) {
    xp <- x
    for (g in unique(groups)) {
      w <- which(groups == g)
      if (length(w) > 1) xp[w] <- sample(x[w])
    }
    s[k] <- coef(lm(rtt ~ xp))[2]
  }
  s
}
res <- list(free = perm(rep("all", length(x))), within_study = perm(study), within_region = perm(region))
cat("\npermutation tests (1000 each): one-sided p = P(permuted rate >= observed)\n")
for (n in names(res)) {
  cat(sprintf("  %-14s p = %.4f   permuted rate: median %.4f, 95%% upper %.4f\n", n,
              mean(res[[n]] >= slope), median(res[[n]]), quantile(res[[n]], .95)))
}
pdf(file.path(out, "temporal_signal.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
cols <- c(Americas = "#1f77b4", Europe = "#d62728", Oceania = "#2ca02c", unknown = "grey")
plot(x, rtt, col = cols[region], pch = 16, cex = .5, xlab = "sampling date",
     ylab = "root-to-tip distance", main = sprintf("%s: R2 %.3f", subset, r2))
abline(fit)
legend("topleft", names(cols), col = cols, pch = 16, cex = .7)
boxplot(res, ylab = "permuted rate", main = "date permutations")
abline(h = slope, col = "red")
dev.off()
saveRDS(rt, file.path(out, "tree_initroot.rds"))
cat("\nwritten:", file.path(out, "temporal_signal.pdf"), "\n")
