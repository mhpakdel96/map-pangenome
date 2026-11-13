#!/bin/bash
#
#SBATCH --job-name=checkm2_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=128GB
#SBATCH --time=1-00:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/02_qc/checkm2/slurm_%j.out
# =============================================================================
# 04_checkm2.sbatch
#
# Genome completeness and contamination for every genome in the cohort.
# This is the metric the project currently has no substitute for: DIAMOND-vs-nr
# tells you what a genome IS, not how complete or how contaminated it is, and
# reviewers of a ~1,800-genome pan-genome paper expect it.
#
# CheckM2 is used rather than CheckM1 — it is much faster on large sets and
# does not need a lineage-specific marker placement step.
#
# ONE-TIME SETUP (login node, needs internet):
#     conda create -n checkm2_env -c bioconda -c conda-forge checkm2
#     conda activate checkm2_env
#     checkm2 database --download --path /cfs/earth/scratch/tabasran/db/checkm2
#
# If the compute nodes have no internet, the database MUST be downloaded on
# the login node first; this job only reads it.
# =============================================================================

set -euo pipefail

module purge
module load slurm

source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
conda activate checkm2_env

ROOT="/cfs/earth/scratch/tabasran/MAP"
FASTA="${ROOT}/03_genomes/fasta"
OUTDIR="${ROOT}/02_qc/checkm2/all_genomes"
DB="/cfs/earth/scratch/tabasran/db/checkm2/CheckM2_database/uniref100.KO.1.dmnd"

export CHECKM2DB="${DB}"

mkdir -p "${OUTDIR}"

echo "=== CheckM2 version ==="
checkm2 --version
echo "======================="

if [[ ! -f "${DB}" ]]; then
    echo "[!] CheckM2 database not found at ${DB}"
    echo "    Run the one-time setup described in the header first."
    exit 1
fi

N=$(ls -1 "${FASTA}"/*.fasta | wc -l)
echo "Running CheckM2 on ${N} genomes..."

# --lowmem halves the DIAMOND memory footprint; drop it if the node has room.
checkm2 predict \
    --threads 30 \
    --input "${FASTA}" \
    --extension .fasta \
    --output-directory "${OUTDIR}" \
    --database_path "${DB}" \
    --force

echo "CheckM2 finished. Key output:"
echo "  ${OUTDIR}/quality_report.tsv"
