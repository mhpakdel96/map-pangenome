#!/usr/bin/env python3
"""
56_lineage_gene_function.py

Functional description of lineage-specific genes, using the read-based
three-state gene presence calls (the final gene-content matrix of the paper).
Descriptive only: no enrichment test, because lineage-specific genes sit in a
few contiguous genomic blocks (large sequence polymorphisms) and are not
independent observations.

Steps
  1. Presence frequency of every gene in each lineage group, among genomes
     with a call. As in 49 (PART B), absences in genomes with K-10 read depth
     below --min-depth are treated as undetermined.
       groups: S/Type I, S/Type III, Bison, C basal grade, C clade D
  2. Lineage-specific genes for five contrasts (present in >= --hi of the
     called genomes on one side and <= --lo on the other; each side needs
     >= --min-called genomes with a call):
       S_only        present in S, absent from C (Bison included in C)
       S_absent      absent from S, present in C
       TypeIII_vs_I  differ between S Type III and S Type I
       Bison_vs_C    differ between Bison and the other C genomes
       D_vs_basal    differ between clade D and the C basal grade
  3. Genomic location: representative sequences mapped with minimap2 to
     K-10 (C-type reference) and Telford (CP033688.1, S Type I reference);
     genes on the same reference closer than --block-gap bp form one block
     (= one independent locus). Genes found in K-10 are annotated with the
     overlapping K-10 CDS (Prokka locus tag and product).
  4. eggNOG-mapper on the translated representative sequences (COG category,
     preferred name, description, KEGG KO/pathway/module, CAZy, PFAM), run
     locally if emapper.py is available, or parsed from a result file
     produced by the eggNOG-mapper web server (--annotations).
  5. KEGG pathway and module names from the KEGG REST API (cached; IDs are
     kept if there is no network), COG category names, and summary tables.

Inputs (defaults under the project root)
  08_reads/summary/presence_reads_ternary.tsv   genes x genomes, 1/0/NA
  08_reads/summary/genome_qc_reads.tsv          K-10 depth per genome
  05_phylogeny/typing/lineage_assignments.tsv   lineage / sublineage
  05_phylogeny/typing/split_assignments.tsv     clade D (tree_side)
  08_reads/pan_ref.fa                           Panaroo pan_genome_reference.fa
  04_pangenome/dna_presence/gene_table.tsv      Panaroo annotation (optional)
  03_genomes/fasta/GCA_000007865.1.fasta        K-10
  03_genomes/gff/GCA_000007865.1.gff            K-10 Prokka annotation
  Telford (CP033688.1)                          searched under 05_phylogeny/typing,
                                                downloaded from ENA if absent

Outputs (07_functional/lineage_genes/)
  group_frequencies.tsv       every gene: presence frequency and n called per group
  lineage_specific_genes.tsv  selected genes with location, block and eggNOG/KEGG
  blocks.tsv                  genomic blocks per contrast, with outgroup presence
  cog_summary.tsv             genes per COG category and contrast
  kegg_pathways.tsv           genes per KEGG pathway and contrast (with gene names)
  selected_genes.fna / .faa   sequences given to eggNOG-mapper
  report.txt                  everything printed

Usage
  # 1) selection, location, sequences (ml_env; minimap2 from snippy_env)
  python 56_lineage_gene_function.py --minimap2 "$(conda run -n snippy_env which minimap2)"
  # 2a) eggNOG-mapper locally (separate env, see below), then parse
  python 56_lineage_gene_function.py --emapper-data /path/to/eggnog_db --cpu 8
  # 2b) or upload 07_functional/lineage_genes/selected_genes.faa to
  #     http://eggnog-mapper.embl.de (proteins, default settings), download
  #     the *.emapper.annotations file and run
  python 56_lineage_gene_function.py --annotations out.emapper.annotations

eggNOG-mapper environment (once):
  conda create -n eggnog_env -c conda-forge -c bioconda eggnog-mapper=2.1.12 -y
  conda activate eggnog_env
  mkdir -p /cfs/earth/scratch/tabasran/db/eggnog
  download_eggnog_data.py -y --data_dir /cfs/earth/scratch/tabasran/db/eggnog
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/cfs/earth/scratch/tabasran/MAP")
GROUPS = ["S/Type I", "S/Type III", "Bison", "C basal grade", "C clade D"]

# COG functional categories (NCBI COG database definitions)
COG = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "-": "No COG category or no eggNOG hit",
}

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]
START = {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}   # bacterial table 11 starts


# ----------------------------------------------------------------------------- helpers
class Log:
    def __init__(self, path):
        self.f = open(path, "w")

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s)
        self.f.write(s + "\n")
        self.f.flush()


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip().upper())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(path, seqs, width=60):
    with open(path, "w") as fh:
        for n, s in seqs.items():
            fh.write(">" + n + "\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


def translate(nt):
    """translation table 11; first codon read as Met if it is a start codon;
    stops: terminal one removed, internal ones written as X (eggNOG cannot
    use '*')."""
    nt = nt[: len(nt) // 3 * 3]
    aa = [CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt), 3)]
    if aa and nt[:3] in START:
        aa[0] = "M"
    if aa and aa[-1] == "*":
        aa = aa[:-1]
    return "".join("X" if a == "*" else a for a in aa)


def safe_name(n):
    """eggNOG-mapper and minimap2 accept '~' but some web forms do not; keep a map."""
    return n.replace("~~~", "__").replace("/", "_").replace(" ", "_")


def find_tool(name):
    return shutil.which(name)


# ----------------------------------------------------------------------------- 1. groups and frequencies
def lineage_groups(lin_path, split_path):
    lin = pd.read_csv(lin_path, sep="\t", dtype=str).set_index("genome")
    sp = pd.read_csv(split_path, sep="\t", dtype=str).set_index("genome")
    g = pd.Series("C basal grade", index=lin.index)
    s = lin["lineage"] == "S"
    g[s & (lin["sublineage"] == "III")] = "S/Type III"
    g[s & (lin["sublineage"] != "III")] = "S/Type I"
    g[(~s) & (lin["sublineage"] == "B")] = "Bison"
    d = sp["tree_side"].reindex(lin.index).fillna("").str.startswith("D ")
    g[d & ~s & (lin["sublineage"] != "B")] = "C clade D"
    return g


def frequencies(X, groups):
    """X genes x genomes with 1/0/NaN; returns freq and n called per group."""
    out = {}
    for name, members in groups.items():
        cols = [c for c in X.columns if c in members]
        sub = X[cols]
        called = sub.notna().sum(axis=1)
        pres = (sub == 1).sum(axis=1)
        out[("freq", name)] = (pres / called.replace(0, np.nan)).round(4)
        out[("n_called", name)] = called
        out[("n_present", name)] = pres
    F = pd.DataFrame(out)
    F.columns = ["{}|{}".format(a, b) for a, b in F.columns]
    return F


# ----------------------------------------------------------------------------- 3. location
def map_genes(minimap2, seqs_fna, ref, tag, outdir, cpu, min_id, min_cov):
    paf = outdir / "map_{}.paf".format(tag)
    cmd = [minimap2, "-x", "asm20", "-c", "--secondary=no", "-t", str(cpu), str(ref), str(seqs_fna)]
    with open(paf, "w") as fh:
        subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL, check=True)
    rows = []
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        qlen, qs, qe, strand, tname, tlen, ts, te, nm, alen = (
            int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10]))
        rows.append(dict(gene=f[0], ref=tag, contig=tname, start=ts + 1, end=te, strand=strand,
                         identity=nm / max(alen, 1), qcov=(qe - qs) / max(qlen, 1), aln=alen))
    P = pd.DataFrame(rows)
    if P.empty:
        return P
    P = P.sort_values(["gene", "aln"], ascending=[True, False]).drop_duplicates("gene")
    return P[(P["identity"] >= min_id) & (P["qcov"] >= min_cov)]


def k10_cds(gff):
    cds = []
    for line in open(gff):
        if line.startswith("##FASTA"):
            break
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[2] != "CDS":
            continue
        at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        cds.append((f[0], int(f[3]), int(f[4]), at.get("locus_tag", ""), at.get("gene", ""),
                    at.get("product", "").replace("%2C", ",")))
    return pd.DataFrame(cds, columns=["contig", "start", "end", "locus_tag", "gene", "product"])


def overlap_cds(loc, C):
    """K-10 CDS with the largest overlap for each mapped gene."""
    out = []
    for _, r in loc.iterrows():
        if r["ref"] != "K-10":
            out.append(("", "", ""))
            continue
        c = C[(C["start"] <= r["end"]) & (C["end"] >= r["start"])]
        if c.empty:
            out.append(("", "", ""))
            continue
        ov = np.minimum(c["end"], r["end"]) - np.maximum(c["start"], r["start"])
        b = c.loc[ov.idxmax()]
        out.append((b["locus_tag"], b["gene"], b["product"]))
    loc = loc.copy()
    loc[["k10_locus_tag", "k10_gene", "k10_product"]] = pd.DataFrame(out, index=loc.index)
    return loc


def blocks(df, gap):
    """block id per gene: same contrast, same reference contig, < gap bp apart."""
    df = df.sort_values(["contrast", "ref", "contig", "start"], na_position="last").copy()
    bid, key0, end0, n, prev_contrast = [], None, None, 0, None
    for _, r in df.iterrows():
        if r["contrast"] != prev_contrast:
            n, key0, end0, prev_contrast = 0, None, None, r["contrast"]
        if pd.isna(r["start"]):
            n += 1
            bid.append("{}_U{:02d}".format(r["contrast"], n))     # unmapped gene: its own block
            key0 = None
            continue
        key = (r["ref"], r["contig"])
        if key != key0 or r["start"] - end0 > gap:
            n += 1
            key0, end0 = key, r["end"]
        else:
            end0 = max(end0, r["end"])
        bid.append("{}_B{:02d}".format(r["contrast"], n))
    df["block"] = bid
    return df


# ----------------------------------------------------------------------------- 4. eggNOG
def run_emapper(faa, outdir, data_dir, cpu, log):
    emapper = find_tool("emapper.py")
    if not emapper:
        log("  [i] emapper.py not on PATH: activate eggnog_env, or use the web server and --annotations")
        return None
    prefix = "lineage_genes"
    cmd = [emapper, "-i", str(faa), "--itype", "proteins", "-m", "diamond", "--cpu", str(cpu),
           "--data_dir", str(data_dir), "--output", prefix, "--output_dir", str(outdir), "--override",
           "--tax_scope", "auto", "--go_evidence", "non-electronic", "--pfam_realign", "none"]
    log("  running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    return outdir / (prefix + ".emapper.annotations")


def read_emapper(path):
    rows, header = [], None
    for line in open(path):
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            header = [h.lstrip("#") for h in f]
            continue
        if header:
            rows.append(f + [""] * (len(header) - len(f)))
    E = pd.DataFrame(rows, columns=header)
    keep = [c for c in ("query", "seed_ortholog", "evalue", "max_annot_lvl", "COG_category", "Description",
                        "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction",
                        "BRITE", "CAZy", "PFAMs") if c in E.columns]
    E = E[keep].rename(columns={"query": "seq_id"})
    E = E.replace("-", "")
    return E.add_prefix("eggnog_").rename(columns={"eggnog_seq_id": "seq_id"})


def kegg_names(outdir, log):
    cache = outdir / "kegg_names.tsv"
    if cache.exists():
        K = pd.read_csv(cache, sep="\t", dtype=str)
        return dict(zip(K["id"], K["name"]))
    names = {}
    for what in ("pathway", "module"):
        try:
            with urllib.request.urlopen("https://rest.kegg.jp/list/" + what, timeout=30) as r:
                for line in r.read().decode().splitlines():
                    i, n = line.split("\t", 1)
                    names[i.replace("path:", "").replace("md:", "")] = n
        except Exception as e:  # noqa: BLE001
            log("  [i] KEGG {} names not fetched ({}); IDs kept".format(what, e))
    if names:
        pd.DataFrame({"id": list(names), "name": list(names.values())}).to_csv(cache, sep="\t", index=False)
    return names


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--hi", type=float, default=0.95, help="presence on the 'present' side (default 0.95)")
    ap.add_argument("--lo", type=float, default=0.05, help="presence on the 'absent' side (default 0.05)")
    ap.add_argument("--min-called", type=int, default=5, help="genomes with a call needed on each side")
    ap.add_argument("--min-depth", type=float, default=30.0, help="absences below this K-10 depth -> NA")
    ap.add_argument("--block-gap", type=int, default=5000, help="max distance (bp) within one block")
    ap.add_argument("--min-id", type=float, default=0.90)
    ap.add_argument("--min-cov", type=float, default=0.80)
    ap.add_argument("--telford", default="", help="Telford CP033688.1 FASTA (searched/downloaded if empty)")
    ap.add_argument("--outgroup", default="", help="outgroup FASTA to polarise gain vs loss "
                                                   "(default 10_outgroup/MAH104.fa from script 51)")
    ap.add_argument("--annotations", default="", help="eggNOG-mapper *.emapper.annotations (web server)")
    ap.add_argument("--emapper-data", default="", help="eggNOG data_dir; runs emapper.py if given")
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--minimap2", default="", help="path to minimap2 if it is not on PATH")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    root = Path(args.root)
    out = Path(args.out) if args.out else root / "07_functional/lineage_genes"
    out.mkdir(parents=True, exist_ok=True)
    log = Log(out / "report.txt")
    p = {
        "ternary": root / "08_reads/summary/presence_reads_ternary.tsv",
        "readqc": root / "08_reads/summary/genome_qc_reads.tsv",
        "lineages": root / "05_phylogeny/typing/lineage_assignments.tsv",
        "split": root / "05_phylogeny/typing/split_assignments.tsv",
        "panref": root / "08_reads/pan_ref.fa",
        "genetable": root / "04_pangenome/dna_presence/gene_table.tsv",
        "k10": root / "03_genomes/fasta/GCA_000007865.1.fasta",
        "k10gff": root / "03_genomes/gff/GCA_000007865.1.gff",
    }
    for k, v in p.items():
        if not v.exists() and k != "genetable":
            sys.exit("[!] missing input {}: {}".format(k, v))

    # ---------------------------------------------------------------- 1
    log("1. GROUP FREQUENCIES (read-based three-state calls)")
    X = pd.read_csv(p["ternary"], sep="\t", index_col=0)
    q = pd.read_csv(p["readqc"], sep="\t", index_col=0)
    low = [g for g in X.columns if q["depth"].get(g, np.nan) < args.min_depth]
    X[low] = X[low].where(X[low] != 0, np.nan)
    log("  genes {} | genomes {} | absences set to NA in {} genomes with K-10 depth < {}x".format(
        X.shape[0], X.shape[1], len(low), args.min_depth))
    grp = lineage_groups(p["lineages"], p["split"])
    grp = grp.reindex(X.columns).dropna()
    X = X[grp.index]
    log("  genomes per group: " + ", ".join("{} {}".format(k, int((grp == k).sum())) for k in GROUPS))
    members = {k: set(grp.index[grp == k]) for k in GROUPS}
    members["S"] = members["S/Type I"] | members["S/Type III"]
    members["C"] = members["Bison"] | members["C basal grade"] | members["C clade D"]
    members["C non-Bison"] = members["C basal grade"] | members["C clade D"]
    F = frequencies(X, members)
    GT = None
    if p["genetable"].exists():
        GT = pd.read_csv(p["genetable"], sep="\t")
        GT = GT.set_index(GT.columns[0])
        keep = [c for c in ("annotation", "len_bp") if c in GT.columns]
        F = F.join(GT[keep])
    F.to_csv(out / "group_frequencies.tsv", sep="\t")

    # ---------------------------------------------------------------- 2
    log("\n2. LINEAGE-SPECIFIC GENES (present >= {} | absent <= {} | >= {} called per side)".format(
        args.hi, args.lo, args.min_called))
    contrasts = [("S_only", "S", "C"), ("S_absent", "C", "S"),
                 ("TypeIII_only", "S/Type III", "S/Type I"), ("TypeIII_absent", "S/Type I", "S/Type III"),
                 ("Bison_only", "Bison", "C non-Bison"), ("Bison_absent", "C non-Bison", "Bison"),
                 ("D_only", "C clade D", "C basal grade"), ("D_absent", "C basal grade", "C clade D")]
    sel = []
    for name, pres, absn in contrasts:
        fp, fa = F["freq|" + pres], F["freq|" + absn]
        ok = ((fp >= args.hi) & (fa <= args.lo) & (F["n_called|" + pres] >= args.min_called)
              & (F["n_called|" + absn] >= args.min_called))
        genes = F.index[ok.fillna(False)]
        log("  {:<15} present in {:<14} absent from {:<14} {:>4} genes".format(name, pres, absn, len(genes)))
        for gname in genes:
            sel.append({"gene": gname, "contrast": name, "present_in": pres, "absent_from": absn,
                        "freq_present_side": fp[gname], "freq_absent_side": fa[gname],
                        "n_called_present_side": int(F.at[gname, "n_called|" + pres]),
                        "n_called_absent_side": int(F.at[gname, "n_called|" + absn])})
    S = pd.DataFrame(sel)
    if S.empty:
        log("  no lineage-specific genes at these thresholds")
        return
    # nested contrasts: a gene that is S_only is not also reported as Type III / Bison etc.
    S = S.drop_duplicates("gene", keep="first")
    if GT is not None:
        S = S.join(GT[[c for c in ("annotation", "len_bp") if c in GT.columns]], on="gene")

    # ---------------------------------------------------------------- 3
    log("\n3. GENOMIC LOCATION")
    pan = read_fasta(p["panref"])
    miss = [g for g in S["gene"] if g not in pan]
    if miss:
        log("  [i] {} selected genes without a representative sequence: {}".format(len(miss), ", ".join(miss[:10])))
    ids = {g: safe_name(g) for g in S["gene"] if g in pan}
    S["seq_id"] = S["gene"].map(ids)
    fna = {ids[g]: pan[g] for g in ids}
    write_fasta(out / "selected_genes.fna", fna)
    faa = {k: translate(v) for k, v in fna.items()}
    write_fasta(out / "selected_genes.faa", faa)
    n_int = sum(1 for v in faa.values() if "X" in v)
    log("  sequences written: {} (protein sequences with internal stop or ambiguous codon: {})".format(len(faa), n_int))

    tel = Path(args.telford) if args.telford else None
    if not tel or not tel.exists():
        cands = sorted(glob.glob(str(root / "05_phylogeny/typing/**/Telford*.f*a"), recursive=True))
        cands = [c for c in cands if not c.endswith(".fai")]
        tel = Path(cands[0]) if cands else out / "Telford_CP033688.1.fasta"
    if not tel.exists():
        try:
            log("  downloading Telford CP033688.1 from ENA")
            urllib.request.urlretrieve("https://www.ebi.ac.uk/ena/browser/api/fasta/CP033688.1", tel)
        except Exception as e:  # noqa: BLE001
            log("  [!] Telford not available ({}); S-only genes will stay unmapped".format(e))
            tel = None
    log("  references: K-10 {} | Telford {}".format(p["k10"], tel))
    mm = args.minimap2 or find_tool("minimap2")
    if not mm:
        sys.exit("[!] minimap2 not on PATH (it is in the environment used for scripts 23/24/34)")
    L = map_genes(mm, out / "selected_genes.fna", p["k10"], "K-10", out, args.cpu, args.min_id, args.min_cov)
    if tel:
        Lt = map_genes(mm, out / "selected_genes.fna", tel, "Telford", out, args.cpu, args.min_id, args.min_cov)
        L = pd.concat([L, Lt[~Lt["gene"].isin(L["gene"])]]) if not L.empty else Lt
    L = overlap_cds(L, k10_cds(p["k10gff"])) if not L.empty else L
    S = S.merge(L.rename(columns={"gene": "seq_id"}), on="seq_id", how="left")
    log("  mapped: K-10 {} | Telford {} | unmapped {}".format(
        int((S["ref"] == "K-10").sum()), int((S["ref"] == "Telford").sum()), int(S["ref"].isna().sum())))
    # ---- outgroup: is the gene present outside MAP? (gain in one lineage vs loss in the other)
    og = Path(args.outgroup) if args.outgroup else root / "10_outgroup/MAH104.fa"
    if og.exists():
        O = map_genes(mm, out / "selected_genes.fna", og, "outgroup", out, args.cpu, args.min_id, args.min_cov)
        hit = set(O["gene"]) if not O.empty else set()
        S["in_outgroup"] = S["seq_id"].isin(hit)
        log("  outgroup {}: {} of {} genes present".format(og.name, int(S["in_outgroup"].sum()), len(S)))
        log("    present in the outgroup = ancestral, so lost in the lineage that lacks it;"
            " absent = gained after MAP split from M. avium")
        for c, d in S.groupby("contrast", sort=False):
            log("    {:<12} {:>3}/{:<3} present in the outgroup".format(c, int(d["in_outgroup"].sum()), len(d)))
    else:
        log("  [i] outgroup not found at {} (run 51_outgroup_state.sh first); gain/loss not polarised".format(og))
        S["in_outgroup"] = np.nan

    S = blocks(S, args.block_gap)
    B = (S.groupby(["contrast", "block"])
          .agg(ref=("ref", "first"), contig=("contig", "first"), start=("start", "min"), end=("end", "max"),
               n_genes=("gene", "size"), n_in_outgroup=("in_outgroup", "sum"),
               k10_from=("k10_locus_tag", lambda x: min([v for v in x if v] or [""])),
               k10_to=("k10_locus_tag", lambda x: max([v for v in x if v] or [""])),
               genes=("gene", lambda x: "; ".join(x)))
          .reset_index())
    B.to_csv(out / "blocks.tsv", sep="\t", index=False)
    for c in S["contrast"].unique():
        b = B[B["contrast"] == c]
        log("  {:<15} {:>4} genes in {:>3} blocks (largest block {} genes)".format(
            c, int(b["n_genes"].sum()), len(b), int(b["n_genes"].max())))

    # ---------------------------------------------------------------- 4
    log("\n4. eggNOG-mapper")
    ann = Path(args.annotations) if args.annotations else None
    if ann is None and args.emapper_data:
        ann = run_emapper(out / "selected_genes.faa", out, args.emapper_data, args.cpu, log)
    if ann is None:
        default = out / "lineage_genes.emapper.annotations"
        ann = default if default.exists() else None
    if ann and Path(ann).exists():
        E = read_emapper(ann)
        S = S.merge(E, on="seq_id", how="left")
        ec = [c for c in S.columns if c.startswith("eggnog_")]
        S[ec] = S[ec].fillna("")
        log("  annotations: {} | genes with an eggNOG hit {}/{}".format(
            ann, int(S["eggnog_seed_ortholog"].fillna("").ne("").sum()), len(S)))
    else:
        log("  [i] no eggNOG annotations yet. Run eggNOG-mapper on {} (see header) and re-run with "
            "--annotations".format(out / "selected_genes.faa"))

    # ---------------------------------------------------------------- 5
    if "eggnog_COG_category" in S.columns:
        log("\n5. SUMMARIES (descriptive)")
        rows = []
        for _, r in S.iterrows():
            cats = r["eggnog_COG_category"] or "-"
            for ch in cats:
                rows.append({"contrast": r["contrast"], "COG": ch, "COG_name": COG.get(ch, ch),
                             "gene": r["gene"], "block": r["block"]})
        C = pd.DataFrame(rows)
        CS = (C.groupby(["contrast", "COG", "COG_name"])
               .agg(genes=("gene", "nunique"), blocks=("block", "nunique"),
                    gene_list=("gene", lambda x: "; ".join(sorted(set(x)))))
               .reset_index().sort_values(["contrast", "genes"], ascending=[True, False]))
        CS.to_csv(out / "cog_summary.tsv", sep="\t", index=False)
        log("  COG categories (genes / independent blocks):")
        for c, d in CS.groupby("contrast", sort=False):
            log("   " + c + ": " + "; ".join("{} {} ({}/{})".format(a, b[:40], g, k) for a, b, g, k in
                                              zip(d["COG"], d["COG_name"], d["genes"], d["blocks"])))
        names = kegg_names(out, log)
        rows = []
        for _, r in S.iterrows():
            for pid in (r.get("eggnog_KEGG_Pathway") or "").split(","):
                if pid.startswith("map"):
                    rows.append({"contrast": r["contrast"], "pathway": pid, "name": names.get(pid, ""),
                                 "gene": r["gene"], "label": r.get("eggnog_Preferred_name") or r["gene"],
                                 "block": r["block"]})
        if rows:
            K = pd.DataFrame(rows)
            KS = (K.groupby(["contrast", "pathway", "name"])
                   .agg(genes=("gene", "nunique"), blocks=("block", "nunique"),
                        gene_list=("label", lambda x: "; ".join(sorted(set(x)))))
                   .reset_index().sort_values(["contrast", "genes"], ascending=[True, False]))
            # broad overview maps (map01100 etc.) are kept but flagged
            KS["overview_map"] = KS["pathway"].str.match(r"map01[1-2]\d\d")
            # maps 03200+ are organismal/disease maps of KEGG: a bacterial gene reaches them
            # only through a shared KO, so they are not interpretable here
            KS["eukaryotic_map"] = KS["pathway"].str.extract(r"map(\d+)")[0].astype(int) >= 3200
            KS.to_csv(out / "kegg_pathways.tsv", sep="\t", index=False)
            log("  KEGG pathways with >= 2 genes (excluding overview and organismal maps):")
            for c, d in KS[(~KS["overview_map"]) & (~KS["eukaryotic_map"])
                           & (KS["genes"] >= 2)].groupby("contrast", sort=False):
                log("   " + c + ": " + "; ".join("{} {} ({} genes, {} blocks: {})".format(
                    a, b, g, k, l) for a, b, g, k, l in zip(d["pathway"], d["name"], d["genes"], d["blocks"],
                                                          d["gene_list"])))
        if "eggnog_KEGG_Module" in S.columns:
            S["eggnog_KEGG_Module_name"] = S["eggnog_KEGG_Module"].fillna("").map(
                lambda s: "; ".join(names.get(m, m) for m in s.split(",") if m))

    cols = [c for c in ("contrast", "block", "gene", "seq_id", "annotation", "len_bp", "present_in", "absent_from",
                        "freq_present_side", "freq_absent_side", "n_called_present_side", "n_called_absent_side",
                        "ref", "contig", "start", "end", "strand", "identity", "qcov", "k10_locus_tag", "k10_gene",
                        "k10_product") if c in S.columns]
    cols += [c for c in S.columns if c.startswith("eggnog_")]
    S = S[cols].sort_values(["contrast", "block", "start"])
    S.to_csv(out / "lineage_specific_genes.tsv", sep="\t", index=False)
    log("\nwritten: " + str(out))


if __name__ == "__main__":
    main()
