#!/bin/bash
#
#SBATCH --job-name=pyseer_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24GB
#SBATCH --time=06:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/06_pangwas/pyseer/logs/pyseer_%A_%a.out
# =============================================================================
# 28_pyseer_run.sbatch
#
# One array task = one line of 06_pangwas/pyseer/manifest.tsv (27_pyseer_prep.py):
# pyseer linear mixed model, phylogenetic similarity as random effect,
# covariates log10 N50 (2q) and platform (3) unless the manifest says 'none'
# (negative control). Variants with 1-99% frequency in the run's samples.
# --output-patterns is written so the number of unique patterns (and hence
# the Bonferroni threshold) can be computed in 29_pyseer_collect.py.
# Restartable: a run whose result file has a header and the log says done is
# skipped.
#
# SUBMIT:
#   mkdir -p 06_pangwas/pyseer/logs 06_pangwas/pyseer/runs
#   N=$(( $(wc -l < 06_pangwas/pyseer/manifest.tsv) - 1 ))
#   sbatch --array=1-${N}%10 scripts/28_pyseer_run.sbatch
# =============================================================================

set -uo pipefail
D="/cfs/earth/scratch/tabasran/MAP/06_pangwas/pyseer"
TASK="${SLURM_ARRAY_TASK_ID:?run as an array job}"
CPU="${SLURM_CPUS_PER_TASK:-4}"

module purge
module load slurm
source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
conda activate pyseer_env
command -v pyseer >/dev/null || { echo "[!] pyseer not on PATH"; exit 1; }

LINE=$(awk -F'\t' -v t="${TASK}" 'NR == t + 1' "${D}/manifest.tsv")
[[ -n "${LINE}" ]] || { echo "[!] no manifest line ${TASK}"; exit 1; }
IFS=$'\t' read -r RUN SUB TRAIT VAR POS NEG COVS <<< "${LINE}"
echo "task ${TASK}: ${RUN} (pos ${POS}, neg ${NEG}, covariates: ${COVS})"
echo "pyseer: $(pyseer --version 2>&1 | tail -1)"

OUT="${D}/runs/${RUN}.tsv"
if [[ -s "${OUT}" ]] && grep -q '__RUN_DONE__' "${D}/runs/${RUN}.log" 2>/dev/null; then
    echo "already done — skipping"; exit 0
fi

P="${D}/pheno/${SUB}__${TRAIT}.tsv"
V="${D}/variants_${VAR}.Rtab"
K="${D}/kinship_${SUB}.tsv"
for f in "$P" "$V" "$K"; do [[ -s "$f" ]] || { echo "[!] missing $f"; exit 1; }; done

# the option name changed between pyseer versions (--similarity / --similarities)
SIM=$(pyseer --help 2>&1 | grep -oE -- '--similarit(y|ies)' | head -1)
[[ -n "${SIM}" ]] || { echo "[!] cannot find the similarity option in pyseer --help"; exit 1; }

COVARGS=()
if [[ "${COVS}" != "none" ]]; then
    read -r -a UC <<< "${COVS}"
    COVARGS=(--covariates "${D}/covariates.tsv" --use-covariates "${UC[@]}")
fi

start=$(date +%s)
pyseer --lmm --phenotypes "$P" --pres "$V" "${SIM}" "$K" "${COVARGS[@]}" \
       --min-af 0.01 --max-af 0.99 \
       --output-patterns "${D}/runs/${RUN}.patterns" --cpu "${CPU}" \
       > "${OUT}" 2> "${D}/runs/${RUN}.log"
rc=$?
echo "exit ${rc} after $(( $(date +%s) - start )) s"
tail -n 8 "${D}/runs/${RUN}.log"
[[ ${rc} -eq 0 && -s "${OUT}" ]] && { echo "__RUN_DONE__" >> "${D}/runs/${RUN}.log"; echo "OK: $(($(wc -l < "${OUT}") - 1)) variants tested"; } \
    || { echo "[!] pyseer failed for ${RUN}"; exit 1; }
