#!/bin/bash
#
#SBATCH --job-name=readmap_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12GB
#SBATCH --time=02:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/08_reads/logs/readmap_%A_%a.out
# =============================================================================
# 34_read_mapping.sbatch
#
# Gene presence and allele counts from READS, independent of assembly.
# One array task = one genome (line of 08_reads/read_paths.tsv).
#
#  A. reads (paired) -> pan_ref.fa (all Panaroo clusters) with minimap2 -x sr;
#     samtools coverage per cluster: reads, covered bases, breadth %, mean
#     depth, mean MAPQ  ->  08_reads/pan_cov/<ID>.tsv
#  B. reads (paired) -> K-10 genome; genome-wide mean depth (for normalising A)
#     and allele depths (AD) at every core SNP site (bcftools mpileup,
#     base and mapping quality >= 20)  ->  08_reads/k10/<ID>.depth,
#     08_reads/k10/<ID>.ad.tsv
#     The AD table shows whether a genome's call at a site is clean (one
#     allele) or MIXED (two alleles with substantial support), which tests
#     the 'homoplastic' SNP patterns and screens for mixed infections.
# BAMs are temporary and deleted. Restartable (finished genomes skipped).
#
# BEFORE SUBMITTING: bash scripts/33_read_inventory.sh ; mkdir -p 08_reads/logs
#                    (snippy_env) samtools faidx 08_reads/k10_ref.fa
# SUBMIT (MaxArraySize: see scontrol show config | grep MaxArraySize):
#         sbatch --array=1-999%60 scripts/34_read_mapping.sbatch
#         sbatch --export=ALL,OFFSET=999 --array=1-<N-999>%60 scripts/34_read_mapping.sbatch
# =============================================================================
set -uo pipefail
ROOT="/cfs/earth/scratch/tabasran/MAP"
D="${ROOT}/08_reads"
K10="${D}/k10_ref.fa"
# The cluster limits array indices (MaxArraySize), so large sets are submitted
# in parts; OFFSET (default 0) shifts the index: line = task + OFFSET
TASK=$(( ${SLURM_ARRAY_TASK_ID:?run as an array job} + ${OFFSET:-0} ))
CPU="${SLURM_CPUS_PER_TASK:-4}"

module purge
module load slurm
unset PERL5LIB PERLLIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
conda activate snippy_env
for t in minimap2 samtools bcftools; do command -v $t >/dev/null || { echo "[!] $t missing"; exit 1; }; done

IFS=$'\t' read -r G R1 R2 U1 U2 <<< "$(awk -F'\t' -v t="${TASK}" 'NR == t + 1' "${D}/read_paths.tsv")"
[[ -n "${G:-}" ]] || { echo "[!] no line ${TASK}"; exit 1; }
mkdir -p "${D}/pan_cov" "${D}/k10"
# the index is made once on the login node (see BEFORE SUBMITTING), never by
# concurrent array tasks
[[ -s "${K10}.fai" ]] || { echo "[!] ${K10}.fai missing — run: samtools faidx ${K10}"; exit 1; }
if [[ -s "${D}/pan_cov/${G}.tsv" && -s "${D}/k10/${G}.ad.tsv" && -s "${D}/k10/${G}.depth" ]]; then
    echo "${G}: done already"; exit 0
fi
TMP=$(mktemp -d "${TMPDIR:-/tmp}/rm_${G}.XXXX")
trap 'rm -rf "${TMP}"' EXIT
echo "${G}: $(date)  minimap2 $(minimap2 --version)  samtools $(samtools --version | head -1 | cut -d' ' -f2)  bcftools $(bcftools --version | head -1 | cut -d' ' -f2)"

# A. pan-genome clusters
minimap2 -ax sr -t "${CPU}" "${D}/pan_ref.fa" "${R1}" "${R2}" 2> "${TMP}/a.log" \
  | samtools sort -@ 2 -o "${TMP}/pan.bam" - \
  && samtools index "${TMP}/pan.bam" \
  && samtools coverage "${TMP}/pan.bam" > "${TMP}/pan_cov.tsv" \
  || { echo "[!] pan mapping failed"; tail "${TMP}/a.log"; exit 1; }

# B. K-10: depth + allele counts at core SNP sites
minimap2 -ax sr -t "${CPU}" "${K10}" "${R1}" "${R2}" 2> "${TMP}/b.log" \
  | samtools sort -@ 2 -o "${TMP}/k10.bam" - \
  && samtools index "${TMP}/k10.bam" \
  || { echo "[!] K-10 mapping failed"; tail "${TMP}/b.log"; exit 1; }
samtools coverage "${TMP}/k10.bam" | awk 'NR==2{print $7"\t"$6}' > "${TMP}/depth"   # meandepth, breadth%
bcftools mpileup -f "${K10}" -T "${D}/core_sites.tsv" -a AD,DP -q 20 -Q 20 -d 10000 \
         "${TMP}/k10.bam" 2> "${TMP}/mp.log" \
  | bcftools query -f '%POS\t%REF\t%ALT\t[%AD]\n' > "${TMP}/ad.tsv" \
  || { echo "[!] mpileup failed"; tail "${TMP}/mp.log"; exit 1; }

mv "${TMP}/pan_cov.tsv" "${D}/pan_cov/${G}.tsv"
mv "${TMP}/ad.tsv" "${D}/k10/${G}.ad.tsv"
mv "${TMP}/depth" "${D}/k10/${G}.depth"
echo "OK ${G}: K-10 depth/breadth $(cat "${D}/k10/${G}.depth") | clusters with reads $(awk 'NR>1 && $4>0' "${D}/pan_cov/${G}.tsv" | wc -l) | AD sites $(wc -l < "${D}/k10/${G}.ad.tsv")"
