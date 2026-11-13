#!/bin/bash
#
#SBATCH --job-name=snippy_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=12:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --array=1-60%30
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/05_phylogeny/logs/snippy_%A_%a.out
# =============================================================================
# 12_snippy.sbatch
#
# Map every genome of the locked cohort onto the K-10 reference with Snippy,
# so that the core alignment used for Gubbins, IQ-TREE, dating and SNP
# distances is in REFERENCE COORDINATES. The Panaroo core alignment is a
# concatenation of genes in arbitrary order, which breaks Gubbins' detection of
# recombination from spatially clustered SNPs.
#
# Assemblies are used for all genomes (--ctgs): this treats de novo and
# reference genomes identically and every input already passed CheckM2.
#
# Array layout: the cohort is split into NCHUNK contiguous chunks, one per
# array task. NCHUNK must equal the upper bound of --array above.
# The script is RESTARTABLE: a genome whose snps.aligned.fa already exists is
# skipped, so failed tasks can simply be resubmitted.
#
# BEFORE SUBMITTING (the log directory must exist, or SLURM kills the job
# silently):
#     mkdir -p /cfs/earth/scratch/tabasran/MAP/05_phylogeny/logs
#
# ONE-TIME SETUP:
#     conda create -y -n snippy_env -c conda-forge -c bioconda \
#           --override-channels snippy=4.6.0
# =============================================================================

set -uo pipefail     # not -e: one failing genome must not stop the others

NCHUNK=60            # <-- must match --array=1-NCHUNK

ROOT="/cfs/earth/scratch/tabasran/MAP"
COHORT="${ROOT}/01_metadata/cohort_final.txt"
FASTA="${ROOT}/03_genomes/fasta"
REF="${FASTA}/GCA_000007865.1.fasta"          # K-10
OUT="${ROOT}/05_phylogeny/snippy"
CPUS="${SLURM_CPUS_PER_TASK:-4}"
TASK="${SLURM_ARRAY_TASK_ID:?run this as an array job}"

module purge
module load slurm
# The cluster's PERL5LIB points conda's Perl at incompatible system modules
# ("Cwd.c: loadable library and perl binaries are mismatched"). Snippy is Perl.
unset PERL5LIB PERLLIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
conda activate snippy_env

command -v snippy >/dev/null || { echo "[!] snippy not on PATH"; exit 1; }
# Fail fast if the environment is broken, instead of failing every genome.
if ! snippy --check > /dev/null 2>&1; then
    echo "[!] 'snippy --check' failed — environment is broken:"
    snippy --check 2>&1 | tail -15
    exit 1
fi
[[ -s "${COHORT}" ]] || { echo "[!] ${COHORT} missing"; exit 1; }
[[ -s "${REF}" ]]    || { echo "[!] reference ${REF} missing"; exit 1; }
mkdir -p "${OUT}"

if [[ "${TASK}" -eq 1 ]]; then
    echo "=== $(snippy --version 2>&1 | head -1) ==="
    echo "reference : ${REF}"
    echo "records   : $(grep -c '>' "${REF}")"
    grep '>' "${REF}" | head -3
fi

N=$(grep -c . "${COHORT}")
PER=$(( (N + NCHUNK - 1) / NCHUNK ))
START=$(( (TASK - 1) * PER + 1 ))
END=$(( START + PER - 1 ))
(( END > N )) && END=${N}
if (( START > N )); then echo "task ${TASK}: nothing to do"; exit 0; fi
echo "task ${TASK}/${NCHUNK}: genomes ${START}-${END} of ${N}"

DONE=0; SKIP=0; FAIL=0
while read -r ID; do
    [[ -z "${ID}" ]] && continue
    D="${OUT}/${ID}"
    if [[ -s "${D}/snps.aligned.fa" ]]; then SKIP=$((SKIP + 1)); continue; fi

    rm -rf "${D}"
    if [[ ! -s "${FASTA}/${ID}.fasta" ]]; then
        echo "[!] ${ID}: input FASTA missing"; FAIL=$((FAIL + 1)); continue
    fi
    snippy --cpus "${CPUS}" --ram 14 --quiet --cleanup \
           --ref "${REF}" --ctgs "${FASTA}/${ID}.fasta" --outdir "${D}" \
           > "${OUT}/${ID}.snippy.log" 2>&1

    if [[ -s "${D}/snps.aligned.fa" ]]; then
        DONE=$((DONE + 1)); rm -f "${OUT}/${ID}.snippy.log"
    else
        echo "[!] ${ID}: snippy failed — see ${OUT}/${ID}.snippy.log"
        FAIL=$((FAIL + 1))
    fi
done < <(sed -n "${START},${END}p" "${COHORT}")

echo "task ${TASK}: done ${DONE}, skipped ${SKIP}, failed ${FAIL}"
[[ "${FAIL}" -eq 0 ]] || exit 1
