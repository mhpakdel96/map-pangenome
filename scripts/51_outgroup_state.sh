#!/bin/bash
# =============================================================================
# 51_outgroup_state.sh
#
# Ancestral (outgroup) state of convergent SNPs: which base does
# M. avium subsp. hominissuis 104 (CP000479.1, an outgroup to MAP) carry at
# the K-10 position? Tells whether recurrent changes are new mutations or
# returns to the ancestral base.
# Run on the LOGIN node (needs internet for the ENA download).
#   bash scripts/51_outgroup_state.sh [K10_POS ...]     (default: icl site)
# =============================================================================
set -uo pipefail
ROOT=/cfs/earth/scratch/tabasran/MAP
OUT=${ROOT}/10_outgroup
K10=${ROOT}/03_genomes/fasta/GCA_000007865.1.fasta
POSITIONS=("${@:-4421235}")
mkdir -p "${OUT}"; cd "${OUT}"
source /cfs/earth/scratch/tabasran/miniconda3/etc/profile.d/conda.sh
conda activate snippy_env
unset PERL5LIB PERLLIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
if [[ ! -s MAH104.fa ]]; then
    curl -s "https://www.ebi.ac.uk/ena/browser/api/fasta/CP000479.1" > MAH104.fa
fi
grep -q '^>' MAH104.fa || { echo "[!] download failed"; exit 1; }
echo "outgroup: $(head -1 MAH104.fa | cut -c1-100) | $(grep -v '>' MAH104.fa | tr -d '\n' | wc -c) bp"
for P in "${POSITIONS[@]}"; do
    S=$((P - 300)); E=$((P + 300))
    samtools faidx "${K10}" "AE016958.1:${S}-${E}" > win_${P}.fa
    minimap2 -c --cs MAH104.fa win_${P}.fa 2>/dev/null > win_${P}.paf
    python3 - "${P}" "${S}" win_${P}.paf win_${P}.fa <<'PY'
import sys, re
pos, start, paf, fa = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
q = "".join(l.strip() for l in open(fa) if not l.startswith(">")).upper()
off = pos - start                         # 0-based offset of the site in the query window
comp = str.maketrans("ACGTN", "TGCAN")
lines = [l.rstrip("\n").split("\t") for l in open(paf) if l.strip()]
if not lines:
    print("{}: no alignment to the outgroup".format(pos)); sys.exit()
best = max(lines, key=lambda f: int(f[9]))
qlen, qs, qe, strand, ts = int(best[1]), int(best[2]), int(best[3]), best[4], int(best[7])
cs = [f[5:] for f in best[12:] if f.startswith("cs:Z:")][0]
# walk the cs string along the target; query coordinate in forward-query space
qpos = qs if strand == "+" else qe - 1
step = 1 if strand == "+" else -1
tpos, found = ts, None
for op, val in re.findall(r"([:*+\-~])([A-Za-z0-9]+)", cs):
    if op == ":":
        n = int(val)
        for _ in range(n):
            if qpos == off:
                found = ("match", None)
            qpos += step; tpos += 1
    elif op == "*":
        tb, qb = val[0].upper(), val[1].upper()
        if qpos == off:
            found = ("subst", tb if strand == "+" else tb.translate(comp))
        qpos += step; tpos += 1
    elif op == "+":
        for _ in val:
            qpos += step
    elif op == "-":
        tpos += len(val)
    if found:
        break
ident = int(best[9]) / int(best[10])
if found is None:
    print("{}: site not covered by the alignment (query {}-{} of {})".format(pos, qs, qe, qlen))
elif found[0] == "match":
    print("{}: outgroup MAH104 has the K-10 base {} (aligned identity {:.3f}, strand {})".format(pos, q[off], ident, strand))
else:
    print("{}: outgroup MAH104 has {} where K-10 has {} (aligned identity {:.3f}, strand {})".format(
        pos, found[1], q[off], ident, strand))
print("  K-10 context: {}[{}]{}".format(q[off-10:off], q[off], q[off+1:off+11]))
PY
done
