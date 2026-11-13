#!/bin/bash
# =============================================================================
# 00_archive_and_restructure.sh
#
# MAP pan-genome project — archive all previous runs and create a clean
# directory structure. NOTHING IS DELETED: every existing output is MOVED
# into 00_archive/ so that earlier results remain inspectable.
#
# Run once, from the project root. Review the printed plan before confirming.
#
# Usage:  bash 00_archive_and_restructure.sh
# =============================================================================

set -euo pipefail

ROOT="/cfs/earth/scratch/tabasran/MAP"
GFF="${ROOT}/Final_gff"
ARCH="${ROOT}/00_archive"

cd "${ROOT}"

echo "==============================================================="
echo " MAP project restructure"
echo " Root:    ${ROOT}"
echo " Archive: ${ARCH}"
echo "==============================================================="
echo

# ---------------------------------------------------------------------------
# 0. Safety check
# ---------------------------------------------------------------------------
if [[ -d "${ARCH}" ]]; then
    echo "[!] ${ARCH} already exists. Aborting so nothing gets overwritten."
    echo "    Rename or remove it manually if you want to re-run this script."
    exit 1
fi

if [[ ! -d "${GFF}" ]]; then
    echo "[!] ${GFF} not found. Check ROOT."
    exit 1
fi

# ---------------------------------------------------------------------------
# 1. Create the new layout
# ---------------------------------------------------------------------------
mkdir -p "${ARCH}"/{roary_run1,roary_run2,panaroo_run1,panaroo_run2,scoary_old,trees_old,logs_old}
mkdir -p "${ROOT}"/01_metadata
mkdir -p "${ROOT}"/02_qc/{quast,checkm2,multiqc}
mkdir -p "${ROOT}"/03_genomes/{fasta,gff,excluded}
mkdir -p "${ROOT}"/04_pangenome
mkdir -p "${ROOT}"/05_phylogeny
mkdir -p "${ROOT}"/06_pangwas
mkdir -p "${ROOT}"/07_functional
mkdir -p "${ROOT}"/08_figures
mkdir -p "${ROOT}"/scripts

echo "[1/5] Directory skeleton created."

# ---------------------------------------------------------------------------
# 2. Archive previous pan-genome / phylogeny / GWAS outputs
# ---------------------------------------------------------------------------
move_if_exists () {   # $1 = source, $2 = destination dir
    if [[ -e "$1" ]]; then
        mv "$1" "$2"/
        echo "      archived: $(basename "$1")"
    fi
}

echo "[2/5] Archiving previous runs..."
move_if_exists "${GFF}/Old_Result_Roary"   "${ARCH}/roary_run1"
move_if_exists "${GFF}/roary.sbatch"       "${ARCH}/roary_run1"
move_if_exists "${GFF}/Output95_Cleaned"   "${ARCH}/roary_run2"
move_if_exists "${GFF}/roary_new.sbatch"   "${ARCH}/roary_run2"

move_if_exists "${GFF}/panaroo_out"        "${ARCH}/panaroo_run1"
move_if_exists "${GFF}/panaroo.sbatch"     "${ARCH}/panaroo_run1"
move_if_exists "${GFF}/panaroo_Rana"       "${ARCH}/panaroo_run2"
move_if_exists "${GFF}/panaroo_Rana.zip"   "${ARCH}/panaroo_run2"
move_if_exists "${GFF}/Panarro_Rana.sbatch" "${ARCH}/panaroo_run2"

move_if_exists "${GFF}/tree.nwk"           "${ARCH}/trees_old"
move_if_exists "${GFF}/nohup.out"          "${ARCH}/logs_old"

# Scoary logs (dated) and SLURM logs
for f in "${GFF}"/scoary_*.log; do move_if_exists "$f" "${ARCH}/scoary_old"; done
for f in "${GFF}"/Slurm-*;      do move_if_exists "$f" "${ARCH}/logs_old";  done
for f in "${GFF}"/slurm-*.out;  do move_if_exists "$f" "${ARCH}/logs_old";  done

# ---------------------------------------------------------------------------
# 3. Move the 23 contaminated GFFs into the new excluded/ folder
# ---------------------------------------------------------------------------
echo "[3/5] Moving contaminated genomes..."
if [[ -d "${GFF}/Contaminated_Genomes_GFF" ]]; then
    mkdir -p "${ROOT}/03_genomes/excluded/step1_gc_outliers"
    mv "${GFF}/Contaminated_Genomes_GFF"/*.gff \
       "${ROOT}/03_genomes/excluded/step1_gc_outliers"/ 2>/dev/null || true
    rmdir "${GFF}/Contaminated_Genomes_GFF" 2>/dev/null || true
    n=$(ls -1 "${ROOT}/03_genomes/excluded/step1_gc_outliers"/*.gff 2>/dev/null | wc -l)
    echo "      ${n} GC-outlier genomes moved to excluded/step1_gc_outliers/"
fi

# ---------------------------------------------------------------------------
# 4. Move the GFFs themselves into 03_genomes/gff/
#    NOTE: GCA_002208705 is present as BOTH .1 and .2 — the older .1 is
#    quarantined here so only one version per genome enters the pan-genome.
# ---------------------------------------------------------------------------
echo "[4/5] Moving GFF files..."
mv "${GFF}"/*.gff "${ROOT}/03_genomes/gff"/ 2>/dev/null || true

mkdir -p "${ROOT}/03_genomes/excluded/step2_duplicates"
if [[ -f "${ROOT}/03_genomes/gff/GCA_002208705.1.gff" && \
      -f "${ROOT}/03_genomes/gff/GCA_002208705.2.gff" ]]; then
    mv "${ROOT}/03_genomes/gff/GCA_002208705.1.gff" \
       "${ROOT}/03_genomes/excluded/step2_duplicates"/
    echo "      GCA_002208705.1 quarantined (superseded by .2)"
fi

# ---------------------------------------------------------------------------
# 5. Inventory
# ---------------------------------------------------------------------------
echo "[5/5] Building inventory..."
ls -1 "${ROOT}/03_genomes/gff"/*.gff 2>/dev/null \
    | xargs -n1 basename | sed 's/\.gff$//' | sort \
    > "${ROOT}/01_metadata/gff_inventory.txt"

N_GFF=$(wc -l < "${ROOT}/01_metadata/gff_inventory.txt")

# leftovers that are not .gff and were not archived
LEFT=$(find "${GFF}" -mindepth 1 -maxdepth 1 2>/dev/null | wc -l)

echo
echo "==============================================================="
echo " GFF files now in 03_genomes/gff/ : ${N_GFF}"
echo " Inventory written to             : 01_metadata/gff_inventory.txt"
echo " Items still left in Final_gff/   : ${LEFT}"
echo "==============================================================="
if [[ "${LEFT}" -gt 0 ]]; then
    echo " Remaining items (inspect and move manually):"
    find "${GFF}" -mindepth 1 -maxdepth 1 -printf '   %f\n'
fi
echo
echo "Expected: 1819 GFFs (1820 minus the duplicated GCA_002208705.1)."
echo "If the number differs, STOP and investigate before running QUAST."
