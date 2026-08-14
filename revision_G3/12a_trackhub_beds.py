#!/usr/bin/env python
"""Write the BED9 files that 12_build_trackhub.sh converts to bigBed.

The hub exists so a reviewer can open any locus in the UCSC browser and see exactly what
this paper measured, on the same assembly (hs1 = T2T-CHM13v2.0). Everything it shows is
derived from the two canonical inputs (plan section 3.0) or from the same gene-set
construction the Gene Ontology analysis uses, so the hub cannot disagree with the figures.

Tracks written
    TEs_<class>.bed        one per TE class, name = subfamily, score = divergence,
                           itemRgb = the project class palette
    TSS_10kb_windows.bed   the 38,704 TSS neighbourhoods the mapping used
    genes_TE_top.bed       the 1,436 genes with the most TEs in their windows
    genes_TE_bottom.bed    the 1,436 genes with the fewest
    IFNA_domain.bed        the 220 kb interferon-alpha domain

Two things worth knowing about the score column. RepeatMasker divergence in this table runs
0-480, not 0-1000, so it is written through unscaled and `trackDb.txt` says so rather than
pretending it is a UCSC-style 0-1000 score. And bigBed requires score <= 1000, which 480
satisfies, so no clamping is needed.

Usage
    python 12a_trackhub_beds.py <output-dir>
"""

import importlib
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import revision_lib as rl  # noqa: E402

# The class palette the figures use, as bigBed itemRgb triples. Keeping the hub and the
# figures on one palette is the point: a reviewer flipping between them sees the same colours.
CLASS_RGB = {
    "LINE": "204,102,11",     # #cc660b
    "LTR": "112,69,60",       # #70453c
    "SINE": "171,31,32",      # #ab1f20
    "DNA": "25,95,144",       # #195f90
    "Retroposon": "118,82,151",  # #765297  (SVA)
    "RC": "35,128,35",        # #238023    (Helitron)
}
# The file uses the RepeatMasker class names; the paper uses SVA and Helitron. The hub keeps
# the paper's names in the track labels and the source names in the data.
CLASS_LABEL = {"Retroposon": "SVA", "RC": "Helitron"}

IFNA = ("chr9", 21_150_692, 21_370_055, "IFNA_domain")
TOP_N = 1436  # same as CLASS_TOP_N in 06_go_rerun_fdr005.py

REPEATS = os.path.join(REPO, "T2T_repeat_masker_processed_sorted.bed")
GENES = os.path.join(REPO, "T2T_genes.bed")
GENOME = os.path.join(REPO, "chm13.genome")


def chromosome_sizes():
    sizes = pd.read_csv(GENOME, sep="\t", header=None, names=["chrom", "size"])
    return dict(zip(sizes["chrom"], sizes["size"]))


def clip_to_chromosomes(frame, sizes, label):
    """Clip interval ends to the chromosome length and starts to zero.

    `T2T_genes.bed` is built by extending 5 kb either side of each TSS without checking the
    chromosome boundary, so a handful of windows run past the end of their chromosome. That is
    harmless for `bedtools intersect`, which simply finds nothing out there, but `bedToBigBed`
    rejects it outright. Clipping is the honest fix: the window really does stop at the end of
    the chromosome. How many were clipped is reported rather than silently corrected.
    """
    limits = frame["chrom"].map(sizes)
    if limits.isna().any():
        unknown = sorted(frame.loc[limits.isna(), "chrom"].unique())
        raise SystemExit(f"{label}: chromosomes absent from {GENOME}: {unknown}")
    over_end = int((frame["end"] > limits).sum())
    under_start = int((frame["start"] < 0).sum())
    if over_end or under_start:
        print(f"    clipped {over_end} interval(s) at a chromosome end and "
              f"{under_start} below zero in {label}")
    out = frame.copy()
    out["end"] = out[["end"]].join(limits.rename("limit")).min(axis=1).astype(int)
    out["start"] = out["start"].clip(lower=0)
    return out[out["end"] > out["start"]]


def write_te_class_beds(out_dir):
    """One BED9 per class, streamed so the 162 MB input is never held in memory twice."""
    columns = ["chrom", "start", "end", "score", "subfamily", "family", "class"]
    handles = {}
    counts = {name: 0 for name in CLASS_RGB}
    try:
        for chunk in pd.read_csv(REPEATS, sep="\t", header=None, names=columns,
                                 chunksize=1_000_000, low_memory=False):
            for class_name, rows in chunk.groupby("class", sort=False):
                if class_name not in CLASS_RGB:
                    raise SystemExit(f"unexpected class in {REPEATS}: {class_name!r}")
                if class_name not in handles:
                    path = os.path.join(out_dir, f"TEs_{CLASS_LABEL.get(class_name, class_name)}.bed")
                    handles[class_name] = open(path, "w")
                frame = pd.DataFrame({
                    "chrom": rows["chrom"],
                    "start": rows["start"],
                    "end": rows["end"],
                    "name": rows["subfamily"],
                    "score": rows["score"].clip(0, 1000),
                    "strand": ".",
                    "thickStart": rows["start"],
                    "thickEnd": rows["end"],
                    "itemRgb": CLASS_RGB[class_name],
                })
                frame.to_csv(handles[class_name], sep="\t", header=False, index=False)
                counts[class_name] += len(rows)
    finally:
        for handle in handles.values():
            handle.close()
    for class_name, n in counts.items():
        print(f"  TEs_{CLASS_LABEL.get(class_name, class_name):10s} {n:>9,} elements")
    return counts


def write_windows_bed(out_dir, sizes):
    windows = pd.read_csv(GENES, sep="\t", header=None,
                          names=["chrom", "start", "end", "gene"])
    windows = clip_to_chromosomes(windows, sizes, "TSS_10kb_windows")
    frame = pd.DataFrame({
        "chrom": windows["chrom"], "start": windows["start"], "end": windows["end"],
        "name": windows["gene"], "score": 0, "strand": ".",
        "thickStart": windows["start"], "thickEnd": windows["end"],
        "itemRgb": "100,100,100",
    })
    path = os.path.join(out_dir, "TSS_10kb_windows.bed")
    frame.to_csv(path, sep="\t", header=False, index=False)
    print(f"  TSS_10kb_windows      {len(frame):>9,} windows")
    return len(frame)


def gene_sets():
    """The TE_top and TE_bottom gene sets, built exactly as 06_go_rerun_fdr005.py builds them.

    Importing that module rather than re-deriving the sets is deliberate: if the construction
    ever changes, the hub changes with it instead of quietly disagreeing with the GO analysis.
    """
    go_rerun = importlib.import_module("06_go_rerun_fdr005")
    table = go_rerun.load_gene_table()
    sets = {}
    for name, ascending in [("TE_top", False), ("TE_bottom", True)]:
        sets[name] = list(
            table.sort_values("TE_number", ascending=ascending)[["Gene_name"]]
            .drop_duplicates()["Gene_name"])[:go_rerun.CLASS_TOP_N]
    return sets


def write_gene_set_beds(out_dir, sizes):
    windows = pd.read_csv(GENES, sep="\t", header=None,
                          names=["chrom", "start", "end", "gene"])
    windows = clip_to_chromosomes(windows, sizes, "gene sets")
    written = {}
    for name, genes in gene_sets().items():
        # A gene with several TSS contributes several windows; all of them are shown, which
        # is the same multiple-TSS behaviour the enrichment analysis has.
        selected = windows[windows["gene"].isin(set(genes))]
        colour = "0,128,0" if name == "TE_top" else "200,0,0"
        frame = pd.DataFrame({
            "chrom": selected["chrom"], "start": selected["start"], "end": selected["end"],
            "name": selected["gene"], "score": 0, "strand": ".",
            "thickStart": selected["start"], "thickEnd": selected["end"],
            "itemRgb": colour,
        })
        path = os.path.join(out_dir, f"genes_{name}.bed")
        frame.to_csv(path, sep="\t", header=False, index=False)
        written[name] = (len(genes), len(frame))
        print(f"  genes_{name:15s} {len(frame):>9,} windows for {len(genes):,} genes")
    return written


def write_ifna_bed(out_dir):
    chrom, start, end, name = IFNA
    path = os.path.join(out_dir, "IFNA_domain.bed")
    with open(path, "w") as handle:
        handle.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t.\t{start}\t{end}\t0,0,0\n")
    print(f"  IFNA_domain           {end - start:>9,} bp")


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    out_dir = sys.argv[1]
    os.makedirs(out_dir, exist_ok=True)
    for path in (REPEATS, GENES, GENOME):
        if not os.path.exists(path):
            raise SystemExit(f"missing canonical input: {path}\n"
                             "See REPRODUCE.md for how to rebuild it.")
    sizes = chromosome_sizes()
    print("writing BED9 files")
    write_te_class_beds(out_dir)
    write_windows_bed(out_dir, sizes)
    write_gene_set_beds(out_dir, sizes)
    write_ifna_bed(out_dir)
    print(f"-> {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
