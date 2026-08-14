# Vector icons for the biological-process groups of Figure 10

Hand-authored SVG icons that replace the seven raster pictures (images picked off the internet)
sitting next to the process labels in **Figure 10**, Figma frame **`861:35`** — the integrated map
(the frame the revision renumbers from `Figure 9`; see `../PLACEMENT.md` §2, caveat C9: renumber
**by node ID only**, because `Figure 9`/`Figure 10` names collide with the subfamilies paper).

Replacing them matters beyond looks: vector art exports at any resolution, and none of it is
third-party artwork of unknown licence. Note that the seven process rasters are **not** all the
pixel content in the frame — `767:1210` (the central insertion cartoon, superseded by the central
panel below) and `768:1284` (the 422 × 421 ring-script image) are images too. Frame `861:35` only
becomes fully vector once all nine are gone; the ring script is out of scope here.

## Which icon replaces which raster

| Process label (text node) | Icon | Raster it replaces | Raster node | Raster size |
|---|---|---|---|---|
| RNA splicing and processing (`774:1312`) | `Fig10_icon_RNA_splicing.svg` | `Group 117` (two stacked images) | `784:87` — children `783:72`, `784:85` | 90 × 89 |
| DNA replication (`776:10`) | `Fig10_icon_DNA_replication.svg` | `image 1529892` | `777:47` | 53 × 107 |
| Metals metabolism (`776:12`) | `Fig10_icon_metals_metabolism.svg` | `image 1529896` | `781:69` | 59 × 71 |
| DNA recombination (`776:6`) | `Fig10_icon_DNA_recombination.svg` | `image 1529900` | `784:89` | 86 × 88 |
| Lipids metabolism (`776:8`) | `Fig10_icon_lipids_metabolism.svg` | `image 1529894` | `780:62` | 60 × 80 |
| Embryogenesis (`773:1286`) | `Fig10_icon_embryogenesis.svg` | `image 1529893` | `779:58` | 82 × 58 |
| Nervous system (`774:1288`) | `Fig10_icon_nervous_system.svg` | `Danya Nikitin PhD theme overview (1) 1` | `780:60` | 95 × 53 |

Two alternates are supplied where the canonical pictogram is a matter of taste — pick one and
delete the other, do not use both:

| Alternate | Instead of | Why it might be preferred |
|---|---|---|
| `Fig10_icon_embryogenesis_alt_blastocyst.svg` | the fetus silhouette | reads unambiguously at very small sizes; says "early development" rather than "fetus" |
| `Fig10_icon_nervous_system_alt_neuron.svg` | the brain outline | a neuron is the cell-level counterpart of the other six icons, which are all molecular |

## What each icon depicts

| Icon | Content |
|---|---|
| RNA splicing | pre-mRNA: two exons (solid blocks) with the intron looped out above; arrow; mature mRNA with the exons joined and the intron gone |
| DNA replication | a replication fork: parent duplex with base-pair rungs unwinding, each parent strand templating a new strand |
| Metals metabolism | a membrane interrupted by a two-helix channel, with metal ions above, inside the pore, and released below |
| DNA recombination | two homologous chromatids crossing over at a chiasma; the break in one strand marks them as separate molecules exchanging arms |
| Lipids metabolism | a phospholipid: polar head group with two hydrocarbon tails |
| Embryogenesis | an embryo in the curled fetal posture — head, curled trunk tapering to a tail, limb buds |
| Nervous system | brain in lateral view with gyri and brainstem |

## The central mechanism panel

`Fig10_center_TE_insertion_epigenetic_impact.svg` — 300 × 290 (`viewBox="25 0 300 290"`), full
colour. It states the paper's mechanism explicitly, which the current centre only implies: a TE
inserts **near a gene, inside the 10 kb window around its TSS**, brings its own epigenetic marks,
and changes both the gene's chromatin profile and its expression.

Three rows, following the reference's own top-to-bottom logic:

| Row | Content |
|---|---|
| 1 | the mobile element — a blue duplex carrying its own histone marks — and a blue arrow to the target site |
| 2 | the locus before insertion: a repressive mark and two DNA-methylation lollipops at the promoter, a small dark-red bent arrow = weak transcription |
| 3 | after insertion: the TE sits upstream inside the window, its active marks have spread to the promoter (dashed arrow), the methylation is gone, and a large green bent arrow = strong transcription |

Symbols and colours are taken from the existing centre, not invented, so the panel is continuous
with the rest of the figure:

| Symbol | Meaning | Colour | Where it came from |
|---|---|---|---|
| Blue duplex / blue rounded box | the TE (free, then inserted) | `#1a7fc1`, dark strand `#0e5f94` | the blue helix in raster `767:1210` |
| Orange boxes on a grey line | host gene exons on genomic DNA | `#e0822f`, outline `#8a4a13`, line `#333333` | the orange host helix in the same raster |
| Pentagon on a stalk | histone mark — blue/green active, red repressive | `#0046a1`, `#00a100`, `#d90202`, black outline | mark groups `767:1231`–`767:1256`, `767:1225` |
| Lollipop | DNA methylation | `#ffcc00`, black outline | groups `767:1213`, `767:1214` |
| Small bent arrow | weak transcription | `#7f0404` | group `767:1268` |
| Large bent arrow | strong transcription | `#2f7f04` | group `767:1267` |
| Bracket with a centre tick | the 10 kb window around the TSS | `#333333` | new — the window is the paper's core design and was not drawn before |

**What it supersedes.** The panel carries its own marks and arrows, so placing it replaces the
central raster **and** the twelve loose vector mark groups sitting on top of it:

```
767:1210  raster (insertion cartoon)
767:1213  767:1214              methylation lollipops
767:1225  767:1226              red pentagons
767:1231  767:1236  767:1241  767:1246  767:1251  767:1256   histone pentagons
767:1267  767:1268              strong / weak transcription arrows
```

Keep the text nodes: `767:1261` "TE insertions epigenetic impact", and reuse `767:1275`
"Weak expression" and `767:1277` "Strong expression" as the row 2 and row 3 labels.

**No text is inside the SVG** — same rule as every other panel, so the figure keeps one typeface
set in Figma. Anchor points for the labels worth adding, in panel coordinates:

| Label | Anchor (x, y) |
|---|---|
| TE with its own marks | 195, 30 |
| TSS | 160, 100 (left of the bent arrow, row 2) |
| Weak expression | 196, 110 |
| DNA methylation | 130, 168 |
| Strong expression | 234, 210 |
| 10 kb window around the TSS | 170, 288 (centred under the bracket) |

Recommended placement size: **~300 × 290 at 1:1**, centred on the footprint of raster `767:1210`
(which is 264 × 156 — the panel is taller because it adds the gene model, the TSS and the window
that the raster never showed, so the centre of the figure needs a little more vertical room).

A monochrome variant is not supplied: this panel belongs to the coloured centre of the figure, not
to the black process-icon set. Ask if one is wanted.

## Conventions

- **Canvas** `viewBox="0 0 64 64"`, `width`/`height` `64` — every icon is square and on the same
  grid, so all seven can be placed at one uniform size (unlike the rasters, whose aspect ratios
  ranged from 53 × 107 to 95 × 53).
- **Colour** pure black `#000000`, **no background element** — the background is transparent. The
  only non-black paint in the set is the white eye dot inside the embryo's black head.
- **Stroke** 3.5 units for primary outlines, 2.5–3 for interior detail, 4 for the two recombination
  chromatids; `stroke-linecap`/`stroke-linejoin` round throughout. At the recommended 58 px
  placement that is a ~3 px primary stroke, which holds up against the figure's heavy arrows and
  16 px labels without going muddy.
- **No `vector-effect`, no clip paths, no embedded raster, no external references** — Figma imports
  these as plain vector networks that scale cleanly.

## Recommended placement

Place each icon at a **uniform 58 × 58** and centre it on the footprint of the raster it replaces —
that is what `_preview_Fig10_with_icons.png` shows, and the label-to-icon distances in the current
layout all still work at that size. Nothing in `revision_G3/` writes to Figma; placement is manual,
same as every other panel in `../PLACEMENT.md`.

In Figma: **Import** the SVG (or paste its contents), then delete the raster node listed above.
Keep the icon's black fill as-is; do not recolour to the TE class palette — these mark biological
processes, not TE classes, and the palette is reserved for classes across both manuscripts.

## Reference renders (not deliverables)

| File | What it is |
|---|---|
| `_preview_icon_sheet.png` | all nine icons, 168 px detail render over the 58 px in-figure size |
| `_preview_Fig10_with_icons.png` | frame `861:35` with the seven rasters whited out and the icons composited at 58 px, for approval before touching the Figma file |
| `_preview_center_panel.png` | the central mechanism panel at 2× |

Regenerate either with `cairosvg` (installed into `~/venvs/Retroelements_3_11`); the icons
themselves are hand-authored and have no build step.
