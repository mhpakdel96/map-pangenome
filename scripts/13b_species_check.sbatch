#!/bin/bash
#
#SBATCH --job-name=species_check_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=06:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/02_qc/species_check/slurm_%j.out
# =============================================================================
# 13b_species_check.sbatch
#
# Confirm that every genome in the cohort is actually M. avium subsp.
# paratuberculosis. CheckM2 measures completeness at genus level and cannot
# tell MAP from other M. avium subspecies; SRR11839093 (deposited as MAP,
# strain field "Avium") passed CheckM2 at 100% but lacks IS900 and sits at
# Mash distance 0.0128 from K-10.
#
# Two independent tests per genome:
#   1. IS900 by SEQUENCE: a K-10 IS900 copy is aligned (bwa mem) to the
#      genome's contigs. IS900 is the defining MAP insertion element (16-22
#      copies) and is absent from other M. avium subspecies. Annotation names
#      are NOT used: Prokka labelled IS900 in only part of the cohort.
#   2. Mash distance to K-10. Within MAP, the most divergent lineage (S-type)
#      sits around 0.002; SRR11839093 sits at 0.0128.
#
# Flags (for review, nothing is removed here):
#   no_IS900  primary alignment < 1000 of 1200 bp, or > 10 mismatches
#   distant   Mash distance > 0.005
#
# BEFORE SUBMITTING:
#     mkdir -p /cfs/earth/scratch/tabasran/MAP/02_qc/species_check
# =============================================================================

set -uo pipefail

ROOT="/cfs/earth/scratch/tabasran/MAP"
COHORT="${ROOT}/01_metadata/cohort_final.txt"
FASTA="${ROOT}/03_genomes/fasta"
REF="${FASTA}/GCA_000007865.1.fasta"
MASKG="${ROOT}/05_phylogeny/core_masked/mask_genes.tsv"
OUT="${ROOT}/02_qc/species_check"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
MASH_MAX=0.005
IS_MIN_ALN=1000
IS_MAX_NM=10

module purge
module load slurm
unset PERL5LIB PERLLIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh

mkdir -p "${OUT}"
cd "${OUT}"
TMPBASE="${OUT}/tmp"; rm -rf "${TMPBASE}"; mkdir -p "${TMPBASE}"

# ---------------------------------------------------- 1. IS900 query
conda activate snippy_env                      # provides bwa and samtools
for t in bwa samtools; do command -v $t >/dev/null || { echo "[!] $t missing"; exit 1; }; done

cp "${REF}" K10.fa && samtools faidx K10.fa
R=$(awk -F'\t' '$4 ~ /IS900/ {print $1":"$2+1"-"$3; exit}' "${MASKG}")
[[ -n "${R}" ]] || { echo "[!] no IS900 coordinates in ${MASKG}"; exit 1; }
samtools faidx K10.fa "${R}" > IS900.fa
IS_LEN=$(grep -v '>' IS900.fa | tr -d '\n' | wc -c)
echo "IS900 query: ${R} (${IS_LEN} bp)"

check_one () {
    local g="$1" d
    d=$(mktemp -d "${TMPBASE}/${g}.XXXX")
    ln -s "$(readlink -f "${FASTA}/${g}.fasta")" "${d}/g.fa"
    bwa index -p "${d}/idx" "${d}/g.fa" >/dev/null 2>&1
    # primary alignments only (-F 2308 drops unmapped, secondary, supplementary)
    bwa mem -t 1 "${d}/idx" IS900.fa 2>/dev/null | samtools view -F 2308 - \
      | awk -v g="${g}" 'BEGIN{best=0; nm="NA"}
          { c=$6; m=0
            while (match(c, /[0-9]+M/)) { m += substr(c, RSTART, RLENGTH-1); c = substr(c, RSTART+RLENGTH) }
            n="NA"; for (i=12; i<=NF; i++) if ($i ~ /^NM:i:/) n=substr($i,6)
            if (m > best) { best=m; nm=n } }
          END { printf "%s\t%d\t%s\n", g, best, nm }'
    rm -rf "${d}"
}
export -f check_one
export FASTA TMPBASE

echo "aligning IS900 to $(grep -c . "${COHORT}") genomes with ${THREADS} parallel jobs..."
grep . "${COHORT}" | xargs -P "${THREADS}" -I{} bash -c 'check_one "$@"' _ {} \
    | sort > is900.tsv
conda deactivate

# ---------------------------------------------------- 2. Mash
conda activate mash_env
command -v mash >/dev/null || { echo "[!] mash missing"; exit 1; }
grep . "${COHORT}" | sed "s|^|${FASTA}/|; s|$|.fasta|" > genome_paths.txt
mash sketch -s 10000 -o k10 K10.fa > /dev/null 2>&1
mash sketch -p "${THREADS}" -s 10000 -o cohort -l genome_paths.txt > /dev/null 2>&1
mash dist -p "${THREADS}" k10.msh cohort.msh \
    | awk -F'\t' '{q=$2; sub(/.*\//,"",q); sub(/\.fasta$/,"",q); print q "\t" $3}' \
    | sort > mash.tsv
conda deactivate

# ---------------------------------------------------- 3. combine
join -t$'\t' is900.tsv mash.tsv \
  | awk -F'\t' -v amin="${IS_MIN_ALN}" -v nmax="${IS_MAX_NM}" -v dmax="${MASH_MAX}" '
      BEGIN { OFS="\t"; print "genome","is900_aligned_bp","is900_mismatches","mash_dist_K10","flag" }
      { f=""
        if ($2 < amin || $3 == "NA" || $3 > nmax) f = "no_IS900"
        if ($4 > dmax) f = (f == "" ? "distant" : f ",distant")
        print $1, $2, $3, $4, (f == "" ? "ok" : f) }' > species_check.tsv
rm -rf "${TMPBASE}"

N_COHORT=$(grep -c . "${COHORT}")
N_ROWS=$(( $(wc -l < species_check.tsv) - 1 ))

echo
echo "--- Mash distance to K-10 (percentiles) ---"
tail -n +2 species_check.tsv | cut -f4 | sort -g | awk '{v[NR]=$1} END{
    split("1 5 25 50 75 95 99 100", p, " ")
    for (i=1; i<=8; i++) { k=int(NR*p[i]/100); if (k<1) k=1; printf "   p%-3s %s\n", p[i], v[k] } }'
echo
echo "--- IS900 alignment ---"
tail -n +2 species_check.tsv | awk -F'\t' '{b=($2>=1150)?"full (>=1150 bp)":($2>=1000)?"partial (1000-1149)":($2>0)?"fragment (<1000)":"none"; c[b]++}
    END{for (k in c) printf "   %-22s %d\n", k, c[k]}' | sort
echo
echo "--- flagged genomes ---"
awk -F'\t' 'NR==1 || $5 != "ok" {printf "   %-18s %10s %10s %14s  %s\n", $1, $2, $3, $4, $5}' species_check.tsv
echo
echo "--- sanity checks ---"
awk -F'\t' '$1=="GCA_000007865.1"{print "   K-10 itself     : aligned "$2" bp, mismatches "$3", mash "$4" -> "$5}' species_check.tsv
awk -F'\t' '$1=="SRR11839093"{print "   SRR11839093     : aligned "$2" bp, mismatches "$3", mash "$4" -> "$5}' species_check.tsv
echo
echo "==============================================================="
echo " genomes checked : ${N_ROWS}  (expected ${N_COHORT})"
echo " flagged         : $(awk -F'\t' 'NR>1 && $5!="ok"' species_check.tsv | wc -l)"
echo " table           : ${OUT}/species_check.tsv"
echo "==============================================================="
[[ "${N_ROWS}" -eq "${N_COHORT}" ]] || { echo "NOT OK — row count mismatch"; exit 1; }
# Positive control: K-10 must carry a full-length, perfect IS900 hit and zero
# Mash distance. If not, the test itself is broken and no flag can be trusted.
K10_FLAG=$(awk -F'\t' '$1=="GCA_000007865.1"{print $5}' species_check.tsv)
[[ "${K10_FLAG}" == "ok" ]] || { echo "NOT OK — positive control K-10 flagged '${K10_FLAG}'; test is broken"; exit 1; }
echo "OK"
