#!/bin/bash
#
#SBATCH --job-name=panaroo_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=128GB
#SBATCH --time=3-00:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/04_pangenome/slurm_panaroo_%j.out
# =============================================================================
# 10_panaroo.sbatch
#
# Pan-genome of the locked cohort (09_apply_filter.py --apply) with Panaroo.
#
# Parameters are kept identical to the previous run so results stay comparable,
# with one difference: the input is now an explicit file list built from
# 01_metadata/cohort_final.txt, not a `*.gff` glob in a directory that also
# held outputs and logs.
#
#   --clean-mode moderate    default graph cleaning
#   -c 0.95                  sequence identity threshold (matches Roary -i 95)
#   --core_threshold 0.99    core = present in >= 99% of genomes
#                            (Roary below uses -cd 99 so the two are comparable)
#   --remove-invalid-genes   drop CDS with internal stops / bad lengths
#   -a core --aligner mafft  core-gene alignment, input for Gubbins
#
# Set PANAROO_ENV to the conda environment that provides panaroo. If panaroo
# comes from a module instead, set PANAROO_ENV="" and add the module load line.
# =============================================================================

set -euo pipefail

PANAROO_ENV="panaroo"          # <-- EDIT after checking `conda env list`

ROOT="/cfs/earth/scratch/tabasran/MAP"
GFF="${ROOT}/03_genomes/gff"
COHORT="${ROOT}/01_metadata/cohort_final.txt"
OUT="${ROOT}/04_pangenome/panaroo"
LIST="${ROOT}/04_pangenome/panaroo_input_gffs.txt"

module purge
module load slurm
source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
[[ -n "${PANAROO_ENV}" ]] && conda activate "${PANAROO_ENV}"

# ------------------------------------------------------------- pre-flight
command -v panaroo >/dev/null || { echo "[!] panaroo not on PATH"; exit 1; }
echo "=== Panaroo version ==="; panaroo --version 2>&1 | head -1; echo

[[ -s "${COHORT}" ]] || { echo "[!] ${COHORT} missing — run 09_apply_filter.py --apply"; exit 1; }
[[ -e "${OUT}" ]] && { echo "[!] ${OUT} exists — move it away first"; exit 1; }

: > "${LIST}"
MISSING=0
while read -r ID; do
    [[ -z "${ID}" ]] && continue
    if [[ -s "${GFF}/${ID}.gff" ]]; then
        echo "${GFF}/${ID}.gff" >> "${LIST}"
    else
        echo "[!] missing GFF: ${ID}"; MISSING=$((MISSING + 1))
    fi
done < "${COHORT}"

N_COHORT=$(grep -c . "${COHORT}")
N_LIST=$(wc -l < "${LIST}")
N_DIR=$(ls -1 "${GFF}"/*.gff | wc -l)
echo "cohort_final.txt : ${N_COHORT}"
echo "GFFs in list     : ${N_LIST}"
echo "GFFs in gff/     : ${N_DIR}"
if [[ "${MISSING}" -ne 0 || "${N_LIST}" -ne "${N_COHORT}" || "${N_DIR}" -ne "${N_COHORT}" ]]; then
    echo "[!] cohort and GFF directory disagree — stopping."; exit 1
fi

# ------------------------------------------------------------------- run
panaroo -i "${LIST}" -o "${OUT}" \
        --clean-mode moderate \
        -c 0.95 \
        --core_threshold 0.99 \
        --remove-invalid-genes \
        -a core --aligner mafft \
        -t 30

# ------------------------------------------------------------ post-check
GPA="${OUT}/gene_presence_absence.csv"
[[ -s "${GPA}" ]] || { echo "[!] ${GPA} not produced"; exit 1; }
# columns: Gene, Non-unique Gene name, Annotation, then one per genome
N_GENOMES=$(head -1 "${GPA}" | awk -F',' '{print NF-3}')
echo
echo "==============================================================="
echo " genomes in gene_presence_absence.csv : ${N_GENOMES}  (expected ${N_COHORT})"
echo "==============================================================="
cat "${OUT}/summary_statistics.txt"
ls -la "${OUT}/core_gene_alignment.aln" 2>/dev/null || echo "[!] core alignment missing"

[[ "${N_GENOMES}" -eq "${N_COHORT}" ]] && echo "OK" || { echo "NOT OK — genome count mismatch"; exit 1; }
