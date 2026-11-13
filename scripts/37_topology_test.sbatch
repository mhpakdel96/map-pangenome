#!/bin/bash
#
#SBATCH --job-name=topotest_MAP
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=xtah@zhaw.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=/cfs/earth/scratch/tabasran/MAP/05_phylogeny/topology_test/logs/topo_%A_%a.out
# =============================================================================
# 37_topology_test.sbatch
#
# Are the 'homoplastic' core SNPs real homoplasy, or is the ML tree wrong?
# Reads (35_read_presence.py, PART 4) showed the alleles are clean (median
# allele fraction 1.0 in carriers, ~no mixed calls), so mis-calls and mixed
# infections are excluded. Gubbins found no recombination at these sites.
# Remaining question: can the carriers of a pattern form ONE clade without a
# significantly worse fit to the data? If yes -> the ML tree is merely
# uncertain there and the SNP marks a clade; if the constrained tree is
# rejected -> the conflict is real (convergence or undetected recombination).
#
# MODE=constraint (array, one task per variant in VARIANTS below):
#   carriers (minor allele, all 1792 genomes) of the variant are forced to be
#   monophyletic: IQ-TREE -g on the same alignment and model as 14
#   (Gubbins-filtered variable sites, GTR+G4+ASC).
# MODE=au (single job, after all constraint tasks):
#   AU test (and KH/SH) of the unconstrained ML tree vs every constrained
#   tree, 10000 RELL replicates.
#
# SUBMIT:
#   mkdir -p 05_phylogeny/topology_test/logs
#   J=$(sbatch --parsable --export=ALL,MODE=constraint --array=1-4 scripts/37_topology_test.sbatch)
#   sbatch --export=ALL,MODE=au --dependency=afterok:$J --cpus-per-task=16 --time=06:00:00 \
#          scripts/37_topology_test.sbatch
# =============================================================================
set -uo pipefail
ROOT="/cfs/earth/scratch/tabasran/MAP"
D="${ROOT}/05_phylogeny/topology_test"
ALN="${ROOT}/05_phylogeny/iqtree/snps_variable.fasta"
MLTREE="${ROOT}/05_phylogeny/iqtree/tree_final.nwk"
SNPS="${ROOT}/06_pangwas/pyseer/variants_snps.Rtab"
MODEL="GTR+G4+ASC"
THREADS="${SLURM_CPUS_PER_TASK:-30}"
# the largest conflicting patterns (31_homoplasy_check.py): the 7-site European
# group, and the three large Americas-dominated patterns
VARIANTS=(AE016958.1_374739_T_C AE016958.1_3312429_G_C AE016958.1_2992754_C_G AE016958.1_529627_C_G)

module purge
module load slurm
source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
conda activate iqtree_env
IQ=""
for b in iqtree3 iqtree2 iqtree; do command -v $b >/dev/null && { IQ=$b; break; }; done
[[ -n "${IQ}" ]] || { echo "[!] no iqtree"; exit 1; }
mkdir -p "${D}"
cd "${D}"
echo "MODE=${MODE:?set MODE=constraint or MODE=au}  $(${IQ} --version 2>&1 | head -1)"

if [[ "${MODE}" == "constraint" ]]; then
    V="${VARIANTS[$(( ${SLURM_ARRAY_TASK_ID:?array} - 1 ))]}"
    P="c_${V}"
    if [[ -s "${P}.treefile" ]]; then echo "${P} done already"; exit 0; fi
    # constraint: (carriers) + every other taxon unresolved; Python 3.5-safe
    python3 - "${SNPS}" "${V}" "${ALN}" "${P}.constraint.nwk" <<'PY'
import sys
snps, var, aln, out = sys.argv[1:5]
taxa = [l[1:].split()[0] for l in open(aln) if l.startswith(">")]
with open(snps) as fh:
    head = fh.readline().rstrip("\n").split("\t")[1:]
    row = None
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[0] == var:
            row = f[1:]
            break
if row is None:
    sys.exit("[!] variant not found: " + var)
x = dict(zip(head, row))
ones = [t for t in taxa if x.get(t) == "1"]
zeros = [t for t in taxa if x.get(t) == "0"]
car, oth = (ones, zeros) if len(ones) <= len(zeros) else (zeros, ones)
miss = [t for t in taxa if t not in x]
if miss:
    sys.exit("[!] taxa without genotype: {}".format(miss[:5]))
with open(out, "w") as o:
    o.write("((" + ",".join(car) + ")," + ",".join(oth) + ");\n")
print("{}: {} carriers (minor allele), {} others".format(var, len(car), len(oth)))
PY
    [[ -s "${P}.constraint.nwk" ]] || exit 1
    ${IQ} -s "${ALN}" -m "${MODEL}" -g "${P}.constraint.nwk" -T "${THREADS}" \
          --prefix "${P}" -redo > "${P}.stdout" 2>&1
    [[ -s "${P}.treefile" ]] || { echo "[!] constrained search failed"; tail -20 "${P}.log"; exit 1; }
    grep -E "^(BEST SCORE FOUND|Log-likelihood of the tree)" "${P}.log" "${P}.iqtree" 2>/dev/null | tail -2
    echo "OK ${P}"
elif [[ "${MODE}" == "au" ]]; then
    cp "${MLTREE}" trees.nwk
    for V in "${VARIANTS[@]}"; do
        [[ -s "c_${V}.treefile" ]] || { echo "[!] missing c_${V}.treefile"; exit 1; }
        cat "c_${V}.treefile" >> trees.nwk
    done
    printf "tree\n1\tML (unconstrained)\n" > trees_index.tsv
    k=2; for V in "${VARIANTS[@]}"; do printf "%s\tconstrained: %s\n" $k "$V" >> trees_index.tsv; k=$((k+1)); done
    # parameters estimated on the ML tree (-te); no new search (-n 0)
    ${IQ} -s "${ALN}" -m "${MODEL}" -z trees.nwk -te "${MLTREE}" -n 0 -zb 10000 -au \
          -T "${THREADS}" --prefix au -redo > au.stdout 2>&1 \
      || ${IQ} -s "${ALN}" -m "${MODEL}" -z trees.nwk -n 0 -zb 10000 -au \
          -T "${THREADS}" --prefix au -redo > au.stdout 2>&1
    cat trees_index.tsv
    sed -n '/USER TREES/,/^$/p;/Tree      logL/,/^$/p' au.iqtree
    grep -A12 'Tree      logL' au.iqtree | head -14
else
    echo "[!] MODE must be constraint or au"; exit 1
fi
