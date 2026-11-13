#!/bin/bash
#
#SBATCH --job-name=lineage_typing_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/05_phylogeny/typing/typing_%j.out
# =============================================================================
# 16_lineage_typing.sbatch
#
# 1. Snippy for the external Type I control Telford (CP033688.1) against K-10,
#    with exactly the settings of 12_snippy.sbatch. Telford is NOT added to the
#    cohort; it is only a positive control for the gyrA/gyrB test.
# 2. 16_lineage_typing.py: lineages from the ML tree + gyrA/gyrB alleles.
#
# BEFORE SUBMITTING (login node; compute nodes may have no internet):
#   mkdir -p 05_phylogeny/typing/external
#   curl -s "https://www.ebi.ac.uk/ena/browser/api/fasta/CP033688.1" \
#        > 05_phylogeny/typing/external/Telford_CP033688.1.fasta
# =============================================================================

set -uo pipefail
ROOT="/cfs/earth/scratch/tabasran/MAP"
EXT="${ROOT}/05_phylogeny/typing/external"
TEL_FA="${EXT}/Telford_CP033688.1.fasta"
REF="${ROOT}/03_genomes/fasta/GCA_000007865.1.fasta"

module purge
module load slurm
unset PERL5LIB PERLLIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh

[[ -s "${TEL_FA}" ]] || { echo "[!] ${TEL_FA} missing — download it first (see header)"; exit 1; }
# ENA writes ">ENA|CP033688|CP033688.1 ..."; check the first header, then
# write a copy with a plain ">CP033688.1" header (no '|' in contig names).
head -1 "${TEL_FA}" | grep -q 'CP033688' || { echo "[!] ${TEL_FA} is not CP033688:"; head -1 "${TEL_FA}"; exit 1; }
[[ "$(grep -c '>' "${TEL_FA}")" -eq 1 ]] || { echo "[!] expected one record in ${TEL_FA}"; exit 1; }
TEL_CLEAN="${EXT}/Telford.clean.fasta"
awk 'NR==1{print ">CP033688.1"; next}{print}' "${TEL_FA}" > "${TEL_CLEAN}"
echo "Telford: 1 record, $(grep -v '>' "${TEL_CLEAN}" | tr -d '\n' | wc -c) bp"
TEL_FA="${TEL_CLEAN}"

conda activate snippy_env
if [[ ! -s "${EXT}/Telford/snps.aligned.fa" ]]; then
    snippy --ctgs "${TEL_FA}" --ref "${REF}" --outdir "${EXT}/Telford" \
           --cpus "${SLURM_CPUS_PER_TASK:-4}" --ram 14 --cleanup --force \
           > "${EXT}/snippy_Telford.log" 2>&1 \
        || { echo "[!] snippy failed"; tail -20 "${EXT}/snippy_Telford.log"; exit 1; }
fi
echo "Telford vs K-10: $(grep -vc '^#' "${EXT}/Telford/snps.vcf" 2>/dev/null) variant records"
conda deactivate

conda activate checkm2_env
python "${ROOT}/scripts/16_lineage_typing.py" \
       --external "Telford_TypeI=${EXT}/Telford/snps.aligned.fa"
