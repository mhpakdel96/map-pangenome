#!/bin/bash
# =============================================================================
# 33_read_inventory.sh
#
# Find the trimmed reads (Trimmomatic output, Final/<ID>/Clean/) for every
# cohort genome and prepare the shared references for 34_read_mapping.sbatch.
#
# Writes (08_reads/)
#   read_paths.tsv   genome, R1_paired, R2_paired, R1_unpaired, R2_unpaired
#                    (only genomes whose two paired files exist and are non-empty)
#   no_reads.txt     cohort genomes without usable reads (the NCBI assemblies
#                    and anything missing) — reported, never silently dropped
#   pan_ref.fa       Panaroo pan_genome_reference.fa (all clusters, so reads go
#                    to their best-matching gene, core or accessory)
#   k10_ref.fa       copy of the K-10 reference (indexed by 34 on first use)
#   core_sites.tsv   CHR<TAB>POS of every core SNP site (snippy-core core.tab),
#                    for allele counts from reads at these positions
# =============================================================================
set -uo pipefail
ROOT="/cfs/earth/scratch/tabasran/MAP"
OUT="${ROOT}/08_reads"
COHORT="${ROOT}/01_metadata/cohort_final.txt"
mkdir -p "${OUT}"

printf "genome\tR1\tR2\tU1\tU2\n" > "${OUT}/read_paths.tsv"
: > "${OUT}/no_reads.txt"
while read -r g; do
    [[ -z "$g" ]] && continue
    c="${ROOT}/Final/${g}/Clean"
    # Trimmomatic output is *.fastq.gz for most runs, *.fq.gz for a few
    ext=fastq.gz
    [[ -s "${c}/${g}_1_paired.fastq.gz" ]] || ext=fq.gz
    r1="${c}/${g}_1_paired.${ext}"; r2="${c}/${g}_2_paired.${ext}"
    u1="${c}/${g}_1_unpaired.${ext}"; u2="${c}/${g}_2_unpaired.${ext}"
    if [[ -s "$r1" && -s "$r2" ]]; then
        printf "%s\t%s\t%s\t%s\t%s\n" "$g" "$r1" "$r2" \
            "$([[ -s $u1 ]] && echo "$u1" || echo NA)" "$([[ -s $u2 ]] && echo "$u2" || echo NA)" \
            >> "${OUT}/read_paths.tsv"
    else
        echo "$g" >> "${OUT}/no_reads.txt"
    fi
done < "${COHORT}"

cp "${ROOT}/04_pangenome/panaroo/pan_genome_reference.fa" "${OUT}/pan_ref.fa"
# private, indexed copy of K-10 (the fasta in 03_genomes is a symlink)
# copied only once, so its .fai (made on the login node) never becomes older
# than the fasta (bcftools would otherwise rebuild it from every array task)
[[ -s "${OUT}/k10_ref.fa" ]] || cp -L "${ROOT}/03_genomes/fasta/GCA_000007865.1.fasta" "${OUT}/k10_ref.fa"
tail -n +2 "${ROOT}/05_phylogeny/core_masked/core.tab" | cut -f1,2 > "${OUT}/core_sites.tsv"

N=$(( $(wc -l < "${OUT}/read_paths.tsv") - 1 ))
echo "cohort            : $(grep -c . "${COHORT}")"
echo "with paired reads : ${N}"
echo "without reads     : $(wc -l < "${OUT}/no_reads.txt")  (of which GCA_*: $(grep -c '^GCA_' "${OUT}/no_reads.txt"))"
grep -v '^GCA_' "${OUT}/no_reads.txt" | head -20 | sed 's/^/   missing: /'
echo "pan_ref.fa        : $(grep -c '>' "${OUT}/pan_ref.fa") sequences"
echo "core_sites.tsv    : $(wc -l < "${OUT}/core_sites.tsv") sites"
echo "disk used by paired reads (sample of 20):"
tail -n +2 "${OUT}/read_paths.tsv" | head -20 | cut -f2,3 | tr '\t' '\n' | xargs du -ch 2>/dev/null | tail -1
