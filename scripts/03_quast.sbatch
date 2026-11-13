#!/bin/bash
#
#SBATCH --job-name=quast_MAP_all
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/02_qc/quast/slurm_%j.out
# =============================================================================
# 03_quast.sbatch
#
# Assembly statistics for EVERY genome in the cohort (de novo + GCA), from the
# single clean set built by 01_collect_genomes.sh.
#
# Reference-free mode on purpose: no -r / -g. Aligning ~1,800 genomes to
# K-10 would add days of runtime and the manuscript only reports contiguity,
# length and GC, none of which need a reference.
#
# --min-contig 500 matches the previously reported methodology.
# =============================================================================

set -euo pipefail

module purge
module load slurm

source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
conda activate quast_env          # <-- must contain QUAST; see note below

ROOT="/cfs/earth/scratch/tabasran/MAP"
FASTA="${ROOT}/03_genomes/fasta"
OUTDIR="${ROOT}/02_qc/quast/all_genomes"

mkdir -p "${OUTDIR}"
cd "${FASTA}"

# Record the exact version — the report previously claimed v5.2.0 while the
# log said 5.0.2. Whatever it prints here is what goes into the manuscript.
echo "=== QUAST version ==="
quast.py --version
echo "====================="

N=$(ls -1 "${FASTA}"/*.fasta | wc -l)
echo "Running QUAST on ${N} genomes..."

# Pass files explicitly rather than with a glob so labels stay clean and
# QUAST does not append "1"/"2" suffixes to duplicate basenames.
quast.py \
    --threads 30 \
    --min-contig 500 \
    --no-plots \
    --no-html \
    --no-icarus \
    --silent \
    -o "${OUTDIR}" \
    "${FASTA}"/*.fasta

echo "QUAST finished. Key outputs:"
echo "  ${OUTDIR}/report.tsv"
echo "  ${OUTDIR}/transposed_report.tsv"
