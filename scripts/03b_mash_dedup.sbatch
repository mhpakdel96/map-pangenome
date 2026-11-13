#!/bin/bash
#
#SBATCH --job-name=mash_dedup_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=64GB
#SBATCH --time=06:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/02_qc/mash/slurm_%j.out
# =============================================================================
# 03b_mash_dedup.sbatch
#
# Detect near-identical genomes WITHIN the cohort.
#
# The GCA audit showed this project has a systematic redundancy problem: the
# same assembly appears under GenBank and RefSeq accessions, and under
# versioned and unversioned names. GCA_002208705 was caught because both
# versions carried the same accession stem. A duplicate that entered under two
# unrelated run accessions (one SRR, one ERR) would NOT have been caught, and
# duplicate genomes inflate core-gene support, distort branch lengths, and
# break the independence assumption Scoary relies on.
#
# Mash sketches every genome and computes pairwise distances. For MAP, a
# monomorphic organism, unrelated isolates still sit well above zero; distances
# near zero mean the same biological sample sequenced or assembled twice.
#
# SETUP:
#     conda create -y -n mash_env -c conda-forge -c bioconda --override-channels mash
#
# Usage:  sbatch 03b_mash_dedup.sbatch
# =============================================================================

set -euo pipefail

module purge
module load slurm

source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
conda activate mash_env

ROOT="/cfs/earth/scratch/tabasran/MAP"
FASTA="${ROOT}/03_genomes/fasta"
OUT="${ROOT}/02_qc/mash"

mkdir -p "${OUT}"
cd "${OUT}"

echo "=== Mash version ==="
mash --version
echo "===================="

N=$(ls -1 "${FASTA}"/*.fasta | wc -l)
echo "Sketching ${N} genomes..."

# -s 10000 : larger sketch than the default 1000, so distances between very
#            similar genomes stay resolvable instead of collapsing to 0.
# -k 21    : standard for bacterial genomes.
ls -1 "${FASTA}"/*.fasta > genome_list.txt
mash sketch -p 30 -s 10000 -k 21 -o cohort -l genome_list.txt

echo "Computing pairwise distances..."
mash dist -p 30 cohort.msh cohort.msh > mash_dist_all.tsv

# Columns: ref, query, distance, p-value, shared-hashes
# Keep only distinct pairs below a permissive threshold, then inspect.
awk -F'\t' '$1 != $2 && $3 < 0.001' mash_dist_all.tsv \
    | sort -t$'\t' -k3,3g > mash_candidate_pairs.tsv

echo
echo "==============================================================="
echo " Genomes sketched      : ${N}"
echo " Pairs with dist<0.001 : $(wc -l < mash_candidate_pairs.tsv)"
echo " (each pair listed twice, A-B and B-A)"
echo "==============================================================="
echo
echo "Closest 40 pairs:"
head -40 mash_candidate_pairs.tsv \
    | awk -F'\t' '{gsub(/.*\//,"",$1); gsub(/.*\//,"",$2);
                   printf "  %-22s %-22s %s\n", $1, $2, $3}'
echo
echo "Full output: ${OUT}/mash_candidate_pairs.tsv"
echo "Note: mash_dist_all.tsv is large; compress or delete it once reviewed."
