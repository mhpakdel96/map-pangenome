#!/bin/bash
#
#SBATCH --job-name=scoary_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --time=08:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/06_pangwas/scoary_runs/logs/scoary_%A_%a.out
# =============================================================================
# 19_scoary_run.sbatch
#
# One array task = one line of 06_pangwas/scoary_runs/manifest.tsv (written by
# 18_scoary_prep.py): one isolate subset x one trait, with its own gene
# matrix, trait file and pruned tree.
#
# Scoary is run with the individual (naive) p-value filter only (-c I -p 0.05)
# so that EVERY gene with naive p < 0.05 is written with all of Scoary's
# columns (Bonferroni_p, Benjamini_H_p, Best/Worst_pairwise_comp_p). The
# reported filter (Bonferroni < 0.05 AND worst pairwise < 0.05) is applied
# afterwards in 20_scoary_collect.py, so the threshold is visible in code and
# can be changed without re-running Scoary.
#
# The tree passed with -n is the IQ-TREE ML tree pruned to the run's isolates;
# Scoary otherwise builds its own UPGMA tree from gene presence/absence.
#
# Restartable: a run directory that already holds a *.results.csv is skipped.
#
# SUBMIT (N = number of runs = lines in manifest minus header):
#   mkdir -p 06_pangwas/scoary_runs/logs
#   N=$(( $(wc -l < 06_pangwas/scoary_runs/manifest.tsv) - 1 ))
#   sbatch --array=1-${N}%20 scripts/19_scoary_run.sbatch
# =============================================================================

set -uo pipefail
ROOT="/cfs/earth/scratch/tabasran/MAP"
RUNS="${ROOT}/06_pangwas/scoary_runs"
MAN="${RUNS}/manifest.tsv"
TASK="${SLURM_ARRAY_TASK_ID:?run as an array job}"

module purge
module load slurm
source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
conda activate scoary_env
command -v scoary >/dev/null || { echo "[!] scoary not on PATH"; exit 1; }

LINE=$(awk -F'\t' -v t="${TASK}" 'NR == t + 1' "${MAN}")
[[ -n "${LINE}" ]] || { echo "[!] no manifest line for task ${TASK}"; exit 1; }
RUN=$(cut -f1 <<< "${LINE}"); TRAIT=$(cut -f3 <<< "${LINE}")
D="${RUNS}/${RUN}"
echo "task ${TASK}: ${LINE}"
echo "scoary: $(scoary --version 2>&1 | tail -1)"

if ls "${D}"/*.results.csv >/dev/null 2>&1; then
    echo "results already present — skipping"; exit 0
fi
for f in genes.csv trait.csv tree.nwk; do
    [[ -s "${D}/${f}" ]] || { echo "[!] missing ${D}/${f}"; exit 1; }
done

# use --no-time only if this Scoary version has it
EXTRA=()
scoary --help 2>&1 | grep -q -- '--no-time' && EXTRA+=(--no-time)

cd "${D}"
start=$(date +%s)
scoary -g genes.csv -t trait.csv -n tree.nwk -s 15 \
       -c I -p 0.05 --threads "${SLURM_CPUS_PER_TASK:-2}" \
       -o "${D}" "${EXTRA[@]}" > scoary.log 2>&1
rc=$?
echo "exit ${rc} after $(( $(date +%s) - start )) s"
tail -n 5 scoary.log
ls "${D}"/*.results.csv >/dev/null 2>&1 \
    && echo "OK: $(ls "${D}"/*.results.csv)" \
    || { echo "[!] no results.csv for ${RUN} (${TRAIT})"; exit 1; }
