#!/bin/bash
#
#SBATCH --job-name=roary_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=128GB
#SBATCH --time=3-00:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/04_pangenome/slurm_roary_%j.out
# =============================================================================
# 11_roary.sbatch
#
# Roary on the SAME locked cohort and with the SAME core threshold as Panaroo
# (-cd 99, identity 95), so that the Roary-vs-Panaroo comparison in the paper
# reflects the algorithms and nothing else. The previous comparison mixed
# -cd 95 (Roary) with --core_threshold 0.99 (Panaroo).
#
# Roary is used only for that comparison; the core alignment, tree and Scoary
# all come from Panaroo. Hence no -e/--mafft here: the alignment step would add
# days of runtime for an output nothing downstream uses.
# =============================================================================

set -euo pipefail

ROOT="/cfs/earth/scratch/tabasran/MAP"
LIST="${ROOT}/04_pangenome/panaroo_input_gffs.txt"   # written by 10_panaroo
COHORT="${ROOT}/01_metadata/cohort_final.txt"
OUT="${ROOT}/04_pangenome/roary"

module purge
module load slurm
unset PERL5LIB PERLLIB
source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
conda activate roary_env

command -v roary >/dev/null || { echo "[!] roary not on PATH"; exit 1; }
echo "=== Roary version ==="; roary --version 2>&1 | tail -1; echo

# Build the list here too, in case Roary runs before or without Panaroo.
if [[ ! -s "${LIST}" ]]; then
    sed "s|^|${ROOT}/03_genomes/gff/|; s|$|.gff|" "${COHORT}" | grep -v '/\.gff$' > "${LIST}"
fi
N_COHORT=$(grep -c . "${COHORT}")
N_LIST=$(wc -l < "${LIST}")
[[ "${N_LIST}" -eq "${N_COHORT}" ]] || { echo "[!] list ${N_LIST} != cohort ${N_COHORT}"; exit 1; }
[[ -e "${OUT}" ]] && { echo "[!] ${OUT} exists — move it away first"; exit 1; }

# Roary refuses to write into an existing directory and needs its own cwd
mkdir -p "$(dirname "${OUT}")"
cd "$(dirname "${OUT}")"

# shellcheck disable=SC2046
roary -p 30 -i 95 -cd 99 -f "${OUT}" $(cat "${LIST}")

GPA="${OUT}/gene_presence_absence.csv"
[[ -s "${GPA}" ]] || { echo "[!] ${GPA} not produced"; exit 1; }
# Roary: 14 annotation columns, then one per genome
N_GENOMES=$(head -1 "${GPA}" | awk -F'","' '{print NF-14}')
echo
echo "==============================================================="
echo " genomes in gene_presence_absence.csv : ${N_GENOMES}  (expected ${N_COHORT})"
echo "==============================================================="
cat "${OUT}/summary_statistics.txt"
[[ "${N_GENOMES}" -eq "${N_COHORT}" ]] && echo "OK" || { echo "NOT OK — genome count mismatch"; exit 1; }
