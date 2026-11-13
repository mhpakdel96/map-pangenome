#!/bin/bash
#
#SBATCH --job-name=gubbins_iqtree_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=128GB
#SBATCH --time=3-00:00:00
#SBATCH --partition=earth-1
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/05_phylogeny/logs/gubbins_iqtree_%j.out
# =============================================================================
# 14_gubbins_iqtree.sbatch
#
# Recombination filtering and maximum-likelihood phylogeny of the locked cohort.
#
# Input: 05_phylogeny/core_masked/clean.full.aln (13_snippy_core.sbatch), a
# whole-genome alignment in K-10 coordinates with PE/PPE and IS regions masked.
# Because positions follow the chromosome, Gubbins can detect recombination
# from spatially clustered SNPs, which it could not do on Panaroo's
# concatenated gene alignment.
#
# Outputs used downstream
#   gubbins/gubbins.final_tree.tre   branch lengths over the whole alignment;
#                                    input for BactDating / dating
#   iqtree/tree_final.nwk            ML topology with UFBoot support; used for
#                                    Scoary, typing (C/S) and figures
#
# IQ-TREE runs on variable sites only, so +ASC (ascertainment-bias correction)
# is required. The previous run used GTR+ASC with no rate heterogeneity; +G4 is
# added here. Recombination masking can leave a Gubbins "polymorphic" site with
# a single nucleotide state; such sites are removed before IQ-TREE, otherwise
# +ASC refuses the alignment.
#
# The job is restartable: if Gubbins already finished, it is not re-run.
#
# Environment: iqtree_env (the environment the previous run used for both
# run_gubbins.py and iqtree). Versions are printed for the Methods.
# =============================================================================

set -uo pipefail

ROOT="/cfs/earth/scratch/tabasran/MAP"
ALN="${ROOT}/05_phylogeny/core_masked/clean.full.aln"
COHORT="${ROOT}/01_metadata/cohort_final.txt"
GDIR="${ROOT}/05_phylogeny/gubbins"
TDIR="${ROOT}/05_phylogeny/iqtree"
THREADS="${SLURM_CPUS_PER_TASK:-30}"
PFX="gubbins"

module purge
module load slurm
source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
conda activate iqtree_env

# ------------------------------------------------------------- pre-flight
command -v run_gubbins.py >/dev/null || { echo "[!] run_gubbins.py not on PATH"; exit 1; }
command -v python3 >/dev/null        || { echo "[!] python3 not on PATH"; exit 1; }
IQ=""
for b in iqtree3 iqtree2 iqtree; do command -v $b >/dev/null && { IQ=$b; break; }; done
[[ -n "${IQ}" ]] || { echo "[!] no iqtree binary on PATH"; exit 1; }

echo "=== versions ==="
run_gubbins.py --version 2>&1 | head -1
${IQ} --version 2>&1 | head -1
echo "================"

[[ -s "${ALN}" ]] || { echo "[!] ${ALN} missing — run 13_snippy_core.sbatch"; exit 1; }
N_COHORT=$(grep -c . "${COHORT}")
N_ALN=$(grep -c '>' "${ALN}")
echo "alignment: ${N_ALN} sequences (cohort ${N_COHORT})"
[[ "${N_ALN}" -eq "${N_COHORT}" ]] || { echo "[!] alignment and cohort disagree"; exit 1; }

# ---------------------------------------------------------------- Gubbins
mkdir -p "${GDIR}"
cd "${GDIR}"
if [[ -s "${PFX}.final_tree.tre" && -s "${PFX}.filtered_polymorphic_sites.fasta" ]]; then
    echo "Gubbins output already present — skipping Gubbins."
else
    echo "Gubbins started: $(date)"
    run_gubbins.py --threads "${THREADS}" --prefix "${PFX}" "${ALN}" \
        > gubbins.log 2>&1
    echo "Gubbins finished: $(date)"
fi
for f in final_tree.tre filtered_polymorphic_sites.fasta recombination_predictions.gff per_branch_statistics.csv; do
    [[ -s "${PFX}.${f}" ]] || { echo "[!] ${PFX}.${f} not produced"; tail -30 gubbins.log 2>/dev/null; exit 1; }
done

N_G=$(grep -c '>' "${PFX}.filtered_polymorphic_sites.fasta")
N_REC=$(grep -vc '^#' "${PFX}.recombination_predictions.gff")
echo "Gubbins: ${N_G} taxa retained, ${N_REC} recombination blocks predicted"
if [[ "${N_G}" -ne "${N_COHORT}" ]]; then
    echo "[!] Gubbins dropped $((N_COHORT - N_G)) taxa (>25% missing data):"
    grep '>' "${PFX}.filtered_polymorphic_sites.fasta" | sed 's/>//' | sort > retained.txt
    sort "${COHORT}" | comm -23 - retained.txt | sed 's/^/   /'
fi

# ------------------------------------------------ variable sites for +ASC
mkdir -p "${TDIR}"
cd "${TDIR}"
python3 - "${GDIR}/${PFX}.filtered_polymorphic_sites.fasta" snps_variable.fasta <<'PY'
import sys
names, seqs, cur = [], [], []
with open(sys.argv[1]) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if names:
                seqs.append("".join(cur))
            names.append(line[1:].split()[0])
            cur = []
        else:
            cur.append(line)
    if names:
        seqs.append("".join(cur))
L = len(seqs[0])
if any(len(s) != L for s in seqs):
    sys.exit("[!] sequences differ in length")
up = [s.upper() for s in seqs]
keep = []
for i in range(L):
    seen = set()
    for s in up:
        c = s[i]
        if c in "ACGT":
            seen.add(c)
            if len(seen) > 1:
                keep.append(i)
                break
with open(sys.argv[2], "w") as out:
    for n, s in zip(names, seqs):
        out.write(">" + n + "\n" + "".join(s[i] for i in keep) + "\n")
print("Gubbins polymorphic sites: {}  ->  variable (>=2 nucleotide states): {}".format(L, len(keep)))
PY
[[ -s snps_variable.fasta ]] || { echo "[!] variable-site alignment not written"; exit 1; }

# ---------------------------------------------------------------- IQ-TREE
echo "IQ-TREE started: $(date)"
${IQ} -s snps_variable.fasta -m GTR+G4+ASC -B 1000 -T "${THREADS}" \
      --prefix iqtree -redo > iqtree.stdout 2>&1
echo "IQ-TREE finished: $(date)"
[[ -s iqtree.treefile ]] || { echo "[!] IQ-TREE produced no tree"; tail -30 iqtree.log 2>/dev/null; exit 1; }
cp iqtree.treefile tree_final.nwk

N_TIPS=$(python3 -c "import re,sys; t=open(sys.argv[1]).read(); print(len(re.findall(r'[(,]([^(),:;]+):', t)))" tree_final.nwk)
echo
echo "--- IQ-TREE ---"
grep -E "^(Log-likelihood of the tree|Best-fit model|Total tree length|Number of parsimony informative)" iqtree.iqtree 2>/dev/null | sed 's/^/   /'
echo
echo "==============================================================="
echo " cohort                 : ${N_COHORT}"
echo " tips in tree_final.nwk : ${N_TIPS}"
echo " Gubbins tree           : ${GDIR}/${PFX}.final_tree.tre"
echo " ML tree                : ${TDIR}/tree_final.nwk"
echo "==============================================================="
[[ "${N_TIPS}" -eq "${N_G}" ]] && echo "OK" || { echo "NOT OK — tip count mismatch"; exit 1; }
