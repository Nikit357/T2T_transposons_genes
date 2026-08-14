# `revision_G3/external/` — third-party inputs, and exactly where each came from

Everything here was downloaded on **2026-08-04** for WP4 (reviewer major comment 4)
and is *not* produced by this repository. Nothing in here is ever modified in place.
If a file is missing, the commands below re-fetch it.

| File | Bytes | md5 | Used by |
|---|---|---|---|
| `lu2020/mmc2.xlsx` | 77,240 | — | `04a_lu2020_geneset_overlap.py` (their Table S1) |
| `lu2020/mmc3.xlsx` … `mmc7.xlsx` | 32,783 – 642,833 | — | not used; kept because the URLs rot |
| `HOM_MouseHumanSequence.rpt` | 15,110,281 | — | `04a` (mouse→human orthology) — **gitignored, re-download** |
| `hs1ToHg38.over.chain.gz` | 2,882,482 | `6652082a3d77877507e69425ce5dfe13` | `04b_newly_resolved_regions.py` |

## Lu et al. 2020 supplementary tables

Lu JY, Shao W, Chang L, et al. *Genomic Repeats Categorize Genes with Distinct
Functions for Orchestrated Regulation.* Cell Rep 2020;30(10):3296-3311.e5.
DOI [10.1016/j.celrep.2020.02.048](https://doi.org/10.1016/j.celrep.2020.02.048).
This is **reference 14** of our manuscript. PMID 32160538, PMCID PMC7195444,
Elsevier PII `S2211124720302096` (taken from the Crossref `alternative-id`, not
guessed — an earlier guess resolved to a different Cell Reports paper entirely).

```bash
PII=S2211124720302096
for i in 2 3 4 5 6 7; do
  curl -L -A "Mozilla/5.0" -o "lu2020/mmc$i.xlsx" \
    "https://ars.els-cdn.com/content/image/1-s2.0-$PII-mmc$i.xlsx"
done
```

`mmc1.pdf` (Document S1, figures) and `mmc8.pdf` (the article) were downloaded and
then deleted — they are 10 MB of non-data.

Routes that do **not** work, recorded so nobody repeats them:

- `pmc.ncbi.nlm.nih.gov/.../bin/NIHMS1574784-supplement-2.xlsx` — PMC now serves a
  JavaScript proof-of-work challenge instead of the file.
- `www.cell.com/cms/...` — 403 to any non-browser client.
- Europe PMC `/supplementaryFiles` — 404 for this PMCID.
- `ftp.ncbi.nlm.nih.gov/pub/pmc/oa_package/...` — 404 despite the OA service
  advertising that path.

### What `mmc2.xlsx` (their Table S1) actually contains

Three columns of **mouse** gene symbols, matching the counts in their Results text
exactly:

| Column | Genes |
|---|---|
| `L1-enriched genes` | 1,480 |
| `Low-complexity-repeat-enriched genes` | 2,439 |
| `SINE-enriched genes` | 2,041 |

Their fourth category — 383 satellite-repeat-enriched genes — is described in the
text but is **not** in Table S1, so it cannot be compared.

These are mm9 categories built on mouse SINE subfamilies (B1, B2, B4). That is why
`04a` has to map through orthology, and why the response letter has to say that a
comparison "on the same dataset" was never available. See the script docstring.

## MGI mouse–human homology

Mouse Genome Informatics, Jackson Laboratory. 24,592 human + 21,929 mouse entries
grouped into homology classes by `DB Class Key`.

```bash
curl -L -o HOM_MouseHumanSequence.rpt \
  "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
```

Gitignored (15 MB, freely re-downloadable). `04a` uses only mouse symbols whose
homology class contains exactly one human symbol — 18,784 of 20,181.

## UCSC T2T-CHM13 → GRCh38 liftOver chains

```bash
curl -L -o hs1ToHg38.over.chain.gz \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz"
```

13,198 chains, 856,823 aligned blocks on the chm13 side. `04b` uses these to bound
how much chm13v2.0 sequence an hg38-based study could not have seen at all. The
derived 208.8 Mb agrees well with the ~182–189 Mb of novel sequence reported by
Nurk et al. 2022, which is the independent check that the parse is right.

## Extra Python dependency

Reading `.xls`/`.xlsx` supplementary tables needs `openpyxl` (already present).
`xlrd` 2.0.2 was installed into `~/venvs/Retroelements_3_11` while chasing a
legacy `.xls` copy of the wrong paper; the final inputs are all `.xlsx`, so
**`xlrd` is not required by any script here** and can be removed if Daniil prefers a
clean environment.
