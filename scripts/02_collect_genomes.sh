#!/bin/bash
# =============================================================================
# 02_collect_genomes.sh   (v3 — explicit paths, no recursive find)
#
# Build ONE clean directory of genome FASTAs — exactly one per genome ID in the
# GFF inventory — so QUAST, Mash and CheckM2 all run on an identical,
# unambiguous cohort.
#
# WHY EXPLICIT PATHS MATTER HERE
# ------------------------------
# The SPAdes output tree contains per-k-mer intermediate assemblies:
#     Spades/<ID>_Spades_Out/K33/final_contigs.fasta
#     Spades/<ID>_Spades_Out/K77/scaffolds.fasta
#     Spades/<ID>_Spades_Out/K127/before_rr.fasta
# Any of these would be silently accepted by a recursive find, and QUAST would
# then report statistics for a partial intermediate graph rather than the final
# assembly. Only the named final outputs below are ever accepted.
#
# Three historical layouts exist because the rename step was applied partially:
#     A  Final/<ID>/Spades_fa/<ID>_contigs.fasta
#     B  Final/<ID>/Spades/<ID>_Spades_Out/<ID>_contigs.fasta
#     C  GCA/<ID>/<ID>.fasta  or  GCA/<ID>/<ID>*_genomic.fna
# The summary reports how many genomes came from each, so any mixing of
# sources is visible rather than hidden.
#
# Usage:  bash 02_collect_genomes.sh
# =============================================================================

set -uo pipefail        # NOT -e: an unresolved genome must be reported, never
                        # abort the run silently.

ROOT="/cfs/earth/scratch/tabasran/MAP"
SPADES_ROOT="${ROOT}/Final"
GCA_ROOT="${ROOT}/GCA"
OUT="${ROOT}/03_genomes/fasta"
META="${ROOT}/01_metadata"
INV="${META}/gff_inventory.txt"

mkdir -p "${OUT}" "${META}"

if [[ ! -s "${INV}" ]]; then
    echo "[!] ${INV} missing. Run 00_archive_and_restructure.sh first."
    exit 1
fi

find "${OUT}" -maxdepth 1 -name '*.fasta' -delete 2>/dev/null
: > "${META}/missing_fasta.txt"
printf 'genome\tsource\tlayout\tresolved_path\n' > "${META}/genome_paths.tsv"

declare -A LAYOUT_COUNT=( [A_spades_fa]=0 [B_spades_out]=0 \
                          [C_gca_fasta]=0 [C_gca_genomic]=0 )

resolve_denovo () {
    # Echoes "<layout>\t<path>" or nothing. Only final assemblies are accepted.
    local id="$1" p
    p="${SPADES_ROOT}/${id}/Spades_fa/${id}_contigs.fasta"
    if [[ -s "${p}" ]]; then printf 'A_spades_fa\t%s' "${p}"; return 0; fi

    p="${SPADES_ROOT}/${id}/Spades/${id}_Spades_Out/${id}_contigs.fasta"
    if [[ -s "${p}" ]]; then printf 'B_spades_out\t%s' "${p}"; return 0; fi

    p="${SPADES_ROOT}/${id}/Spades/${id}_Spades_Out/contigs.fasta"
    if [[ -s "${p}" ]]; then printf 'B_spades_out\t%s' "${p}"; return 0; fi

    return 0
}

resolve_reference () {
    local id="$1" base="${1%.*}" d p
    for d in "${GCA_ROOT}/${id}" "${GCA_ROOT}/${base}"; do
        [[ -d "${d}" ]] || continue
        p="${d}/${id}.fasta"
        if [[ -s "${p}" ]]; then printf 'C_gca_fasta\t%s' "${p}"; return 0; fi
        p="${d}/${base}.fasta"
        if [[ -s "${p}" ]]; then printf 'C_gca_fasta\t%s' "${p}"; return 0; fi
        # NCBI naming: <ACC>_<assembly name>_genomic.fna, directly inside the
        # genome directory only (depth 1).
        p=$(find "${d}" -maxdepth 1 -type f -name '*_genomic.fna' \
              2>/dev/null | sort | head -n1)
        if [[ -n "${p}" && -s "${p}" ]]; then printf 'C_gca_genomic\t%s' "${p}"; return 0; fi
    done
    return 0
}

echo "Resolving one FASTA per genome ID..."
LINE=0
while read -r ID; do
    [[ -z "${ID}" ]] && continue
    LINE=$((LINE + 1))
    (( LINE % 400 == 0 )) && echo "  ... ${LINE} processed"

    if [[ "${ID}" == GC[AF]_* ]]; then
        RES=$(resolve_reference "${ID}"); SRC="reference"
    else
        RES=$(resolve_denovo "${ID}");    SRC="de_novo"
    fi

    if [[ -n "${RES}" ]]; then
        LAYOUT="${RES%%$'\t'*}"
        FA="${RES#*$'\t'}"
        ln -sf "${FA}" "${OUT}/${ID}.fasta"
        printf '%s\t%s\t%s\t%s\n' "${ID}" "${SRC}" "${LAYOUT}" "${FA}" \
            >> "${META}/genome_paths.tsv"
        LAYOUT_COUNT[${LAYOUT}]=$(( ${LAYOUT_COUNT[${LAYOUT}]} + 1 ))
    else
        echo "${ID}" >> "${META}/missing_fasta.txt"
    fi
done < "${INV}"
echo "  ... ${LINE} processed (done)"

# Guard: no accepted path may sit inside a per-k-mer intermediate directory.
KMER_HITS=$(awk -F'\t' 'NR>1 && $4 ~ /\/K[0-9]+\//' "${META}/genome_paths.tsv" | wc -l)

BROKEN=0
shopt -s nullglob
for f in "${OUT}"/*.fasta; do
    if [[ ! -s "${f}" ]]; then
        echo "[!] empty/broken: $(basename "${f}")"
        BROKEN=$((BROKEN + 1))
    fi
done
shopt -u nullglob

DUPSRC=$(awk -F'\t' 'NR>1{print $4}' "${META}/genome_paths.tsv" | sort | uniq -d | wc -l)

N_WANT=$(wc -l < "${INV}")
N_GOT=$(find "${OUT}" -maxdepth 1 -name '*.fasta' | wc -l)
N_MISS=$(wc -l < "${META}/missing_fasta.txt")
N_REF=$(awk -F'\t' '$2=="reference"' "${META}/genome_paths.tsv" | wc -l)
N_DEN=$(awk -F'\t' '$2=="de_novo"'   "${META}/genome_paths.tsv" | wc -l)

echo
echo "==============================================================="
echo " Genome IDs expected      : ${N_WANT}"
echo " Resolved                 : ${N_GOT}  (de novo ${N_DEN}, reference ${N_REF})"
echo " Missing                  : ${N_MISS}"
echo " Broken/empty links       : ${BROKEN}"
echo " Same file via two IDs    : ${DUPSRC}"
echo " Paths inside a K* dir    : ${KMER_HITS}   (must be 0)"
echo "---------------------------------------------------------------"
echo " Layout breakdown:"
for k in A_spades_fa B_spades_out C_gca_fasta C_gca_genomic; do
    printf '   %-16s %s\n' "${k}" "${LAYOUT_COUNT[${k}]}"
done
echo "==============================================================="

if [[ "${N_MISS}" -gt 0 ]]; then
    echo
    echo " Unresolved IDs (${META}/missing_fasta.txt):"
    cat "${META}/missing_fasta.txt" | sed 's/^/   /'
fi
if [[ "${DUPSRC}" -gt 0 ]]; then
    echo
    echo " [!] Same source file reached via two IDs:"
    awk -F'\t' 'NR>1{print $4}' "${META}/genome_paths.tsv" | sort | uniq -d | sed 's/^/   /'
fi
if [[ "${KMER_HITS}" -gt 0 ]]; then
    echo
    echo " [!] Intermediate k-mer assemblies were selected — DO NOT PROCEED:"
    awk -F'\t' 'NR>1 && $4 ~ /\/K[0-9]+\//{print "   "$1"  "$4}' "${META}/genome_paths.tsv"
fi

if [[ "${N_GOT}" -eq "${N_WANT}" && "${N_MISS}" -eq 0 && "${BROKEN}" -eq 0 \
      && "${DUPSRC}" -eq 0 && "${KMER_HITS}" -eq 0 ]]; then
    echo
    echo " OK — cohort complete and unambiguous. Safe to submit QUAST/Mash/CheckM2."
    exit 0
else
    echo
    echo " NOT OK — resolve the issues above before submitting any jobs."
    exit 1
fi
