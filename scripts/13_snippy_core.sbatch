#!/bin/bash
#
#SBATCH --job-name=snippy_core_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/05_phylogeny/logs/snippy_core_%j.out
# =============================================================================
# 13_snippy_core.sbatch
#
# Combine the per-genome Snippy runs (12_snippy.sbatch) into one whole-genome
# alignment in K-10 coordinates, ready for Gubbins.
#
# Steps
#   1. completeness check: every cohort genome must have snps.aligned.fa
#   2. mask BED from the K-10 Prokka GFF: PE/PPE families, transposases /
#      insertion sequences, integrases and phage genes. In mycobacteria these
#      repetitive regions attract mis-mapped fragments and produce false SNPs
#      (IS900 alone has 16-22 copies in MAP). Masking them is standard for
#      M. tuberculosis and is reported in full below for review.
#   3. snippy-core with the mask
#   4. snippy-clean_full_aln (non-ACGT -> N), as recommended for Gubbins
#   5. drop the 'Reference' row snippy-core adds: K-10 is already in the cohort
#      as GCA_000007865.1, and tree tips must match the Scoary trait table
#   6. per-genome alignment QC from core.txt
#
# Set MASK=0 to build the unmasked alignment instead (e.g. for a sensitivity
# check); outputs then go to a separate directory.
# =============================================================================

set -uo pipefail

MASK=1

ROOT="/cfs/earth/scratch/tabasran/MAP"
COHORT="${ROOT}/01_metadata/cohort_final.txt"
REF="${ROOT}/03_genomes/fasta/GCA_000007865.1.fasta"
REF_GFF="${ROOT}/03_genomes/gff/GCA_000007865.1.gff"
SNIPPY="${ROOT}/05_phylogeny/snippy"
if [[ "${MASK}" -eq 1 ]]; then CORE="${ROOT}/05_phylogeny/core_masked"; else CORE="${ROOT}/05_phylogeny/core_unmasked"; fi

module purge
module load slurm
# The cluster's PERL5LIB breaks conda's Perl (snippy-core is Perl).
unset PERL5LIB PERLLIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
conda activate snippy_env

for t in snippy-core snippy-clean_full_aln; do
    command -v "$t" >/dev/null || { echo "[!] $t not on PATH"; exit 1; }
done
snippy-core --version > /dev/null 2>&1 \
    || { echo "[!] snippy-core does not run:"; snippy-core --version 2>&1 | tail -5; exit 1; }
[[ -e "${CORE}" ]] && { echo "[!] ${CORE} exists — move it away first"; exit 1; }
mkdir -p "${CORE}"
cd "${CORE}"

# ------------------------------------------------------------ 1. inputs
N_COHORT=$(grep -c . "${COHORT}")
: > dirs.txt; : > missing.txt
while read -r ID; do
    [[ -z "${ID}" ]] && continue
    if [[ -s "${SNIPPY}/${ID}/snps.aligned.fa" ]]; then
        echo "${SNIPPY}/${ID}" >> dirs.txt
    else
        echo "${ID}" >> missing.txt
    fi
done < "${COHORT}"
N_OK=$(wc -l < dirs.txt); N_MISS=$(wc -l < missing.txt)
echo "cohort ${N_COHORT} | snippy outputs found ${N_OK} | missing ${N_MISS}"
if [[ "${N_MISS}" -gt 0 ]]; then
    echo "[!] missing genomes (resubmit 12_snippy.sbatch; it skips finished ones):"
    head -30 missing.txt | sed 's/^/   /'
    exit 1
fi

# ------------------------------------------------------------- 2. mask
REF_IDS=$(grep '>' "${REF}" | sed 's/^>//; s/[[:space:]].*//' | sort -u)
MASK_ARG=()
if [[ "${MASK}" -eq 1 ]]; then
    [[ -s "${REF_GFF}" ]] || { echo "[!] ${REF_GFF} missing"; exit 1; }
    awk -F'\t' '
        /^##FASTA/ {exit}
        /^#/ || NF < 9 {next}
        $3 == "CDS" {
            p = $9; sub(/.*product=/, "", p); sub(/;.*/, "", p)
            l = tolower(p)
            if (l ~ /(^|[^a-z])ppe([^a-z]|$)/ || l ~ /(^|[^a-z])pe([^a-z]|$)/ ||
                l ~ /pe[_-]pgrs/ || l ~ /pe family/ || l ~ /ppe family/ ||
                l ~ /transposase/ || l ~ /insertion (sequence|element)/ ||
                l ~ /integrase/ || l ~ /phage/)
                print $1 "\t" ($4 - 1) "\t" $5 "\t" p
        }' "${REF_GFF}" | sort -k1,1 -k2,2n > mask_genes.tsv

    GFF_IDS=$(cut -f1 mask_genes.tsv | sort -u)
    for s in ${GFF_IDS}; do
        grep -qx "${s}" <<< "${REF_IDS}" || {
            echo "[!] GFF seqid '${s}' is not a sequence in ${REF}"
            echo "    reference ids: ${REF_IDS}"
            exit 1; }
    done

    # merge overlaps -> BED for snippy-core
    awk -F'\t' 'BEGIN{OFS="\t"}
        { if ($1==c && $2<=e) { if ($3>e) e=$3 }
          else { if (c!="") print c,s,e; c=$1; s=$2; e=$3 } }
        END { if (c!="") print c,s,e }' mask_genes.tsv > mask.bed

    REF_LEN=$(grep -v '>' "${REF}" | tr -d '\n\r' | wc -c)
    MASK_BP=$(awk '{t+=$3-$2} END{print t+0}' mask.bed)
    echo
    echo "--- mask (from ${REF_GFF##*/}) ---"
    echo " genes masked : $(wc -l < mask_genes.tsv)"
    echo " regions      : $(wc -l < mask.bed)"
    echo " bp masked    : ${MASK_BP} of ${REF_LEN} ($(awk -v a="${MASK_BP}" -v b="${REF_LEN}" 'BEGIN{printf "%.2f", 100*a/b}')%)"
    echo " products:"
    cut -f4 mask_genes.tsv | sort | uniq -c | sort -rn | head -25 | sed 's/^/   /'
    [[ "$(wc -l < mask.bed)" -gt 0 ]] || { echo "[!] empty mask — check the GFF"; exit 1; }
    MASK_ARG=(--mask mask.bed)
fi

# ------------------------------------------------------- 3. snippy-core
echo
echo "=== $(snippy-core --version 2>&1 | head -1) ==="
# shellcheck disable=SC2046
snippy-core --ref "${REF}" --prefix core "${MASK_ARG[@]}" $(cat dirs.txt) \
    > snippy_core.log 2>&1 || { echo "[!] snippy-core failed"; tail -20 snippy_core.log; exit 1; }
for f in core.full.aln core.aln core.txt; do
    [[ -s "$f" ]] || { echo "[!] $f not produced"; tail -20 snippy_core.log; exit 1; }
done

# ------------------------------------------------- 4-5. clean, drop ref
snippy-clean_full_aln core.full.aln > clean.full.with_ref.aln
awk '/^>/ {keep = ($0 != ">Reference")} keep' clean.full.with_ref.aln > clean.full.aln

N_SEQ=$(grep -c '>' clean.full.aln)
LENS=$(awk '/^>/ {if (n) print n; n=0; next} {n += length($0)} END {print n}' clean.full.aln | sort -u)
N_LENS=$(wc -l <<< "${LENS}")

# --------------------------------------------------------- 6. QC report
# core.txt: ID LENGTH ALIGNED UNALIGNED VARIANT HET MASKED LOWCOV
awk -F'\t' 'NR==1 || $1=="Reference" {next}
            {printf "%s\t%.2f\t%s\t%s\n", $1, 100*$3/$2, $5, $6}' core.txt \
    | sort -t$'\t' -k2,2g > per_genome_aligned.tsv
echo
echo "--- per-genome alignment (lowest aligned % first) ---"
printf '   %-18s %8s %8s %6s\n' genome aligned% variants het
head -15 per_genome_aligned.tsv | awk -F'\t' '{printf "   %-18s %8s %8s %6s\n",$1,$2,$3,$4}'
echo
awk -F'\t' '{a=$2+0; if(a<90)x++; if(a<95)y++; if(a<75)z++} END{
    printf " aligned < 95%%: %d   < 90%%: %d   < 75%% (Gubbins default drop): %d\n", y, x, z}' per_genome_aligned.tsv
echo " variant (SNP-only) sites in core.aln: $(awk '/^>/{if(s)exit; s=1; next} {n+=length($0)} END{print n}' core.aln)"

echo
echo "==============================================================="
echo " sequences in clean.full.aln : ${N_SEQ}  (expected ${N_COHORT})"
echo " distinct sequence lengths   : ${N_LENS}  (${LENS//$'\n'/ })"
echo " output                      : ${CORE}/clean.full.aln"
echo "==============================================================="
if [[ "${N_SEQ}" -eq "${N_COHORT}" && "${N_LENS}" -eq 1 ]]; then
    echo "OK"
else
    echo "NOT OK"; exit 1
fi
