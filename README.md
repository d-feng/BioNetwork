# BioNetwork

R pipelines for CRISPRi/CRISPRko screen integration and immune cell pathway visualisation.

---

## Motivation: beyond force-directed network layouts

Standard network visualisation tools (Cytoscape, igraph, STRING) use **force-directed layouts** where node position is determined by edge weight, not biology. These layouts are dense, non-reproducible across runs, and collapse to unreadable clusters at publication figure size. This pipeline addresses those limitations with a curated, biologically anchored design:

| Conventional approach | This pipeline |
|---|---|
| Force-directed / arbitrary node positions | Curated x/y coordinates follow **biological compartments** (membrane → cytoplasm → nucleus) |
| All genes rendered equally | **Hit-centric**: only significant CRISPRi hits are full nodes; non-hit connectors are minimal grey bridges |
| Single colour per node (up/down) | **5-state concordance encoding** across N screens |
| One screen at a time | **N-screen merge** with Fisher's combined p-value rewarding replication |
| No domain annotation | **SH2-domain halo** and other structural features overlaid directly on nodes |
| Separate figures per screen | **Dot-matrix** below each node shows per-assay direction in a single figure |

---

## Design principles

### 1. Biologically meaningful layout
Nodes are placed on a fixed canvas reflecting actual cell biology:
- **Top row** — membrane receptors (TCR complex, co-receptors, checkpoint)
- **Middle rows** — cytoplasmic signalling cascades (kinase → adaptor → effector)
- **Bottom rows** — nuclear transcription factors

Modules are spatially separated: TCR/Proximal → MAPK/NFAT → PI3K/AKT → NF-κB → JAK-STAT → Epi/Apoptosis → Mito/Metabolic.
A biologist can instantly orient within the figure without reading labels.

### 2. Hit-centric node rendering
Only genes passing the FDR/LFC threshold appear as full-size coloured nodes.
Non-hit "connector" nodes (e.g. RAS, IKKα) are drawn small and grey — present only to preserve pathway topology between two hit nodes.
Node **radius scales with −log₁₀(Fisher p)**, so the most replicated hits are visually largest.

### 3. Multi-screen concordance encoding
When two or more MAGeCK screens are merged, each node receives one of five colours:

| Colour | Meaning |
|---|---|
| Deep blue | Concordant GoF — all screens agree (depleted) |
| Deep red | Concordant LoF — all screens agree (enriched) |
| Per-assay blue/green/purple… | Single-assay GoF hit only |
| Per-assay red/orange/pink… | Single-assay LoF hit only |
| Gold | Discordant — screens disagree on direction |

### 4. Dot-matrix per-assay indicator
A row of small dots below each hit node encodes per-screen direction independently:
- Position (left → right) = screen order
- Dot colour = direction in that screen (blue = GoF, red = LoF, grey = not a hit)

This lets the reader see at a glance: *"this gene was a GoF hit in both screens"* vs *"only screen 1 called it"* — without needing a separate supplementary figure.

### 5. Fisher's combined p-value
Meta-significance is computed across all screens where a gene is a hit using Fisher's method:
`χ² = −2 Σ log(p_i)`, `df = 2k`. Genes replicated in multiple screens receive a smaller combined p and appear as larger nodes.

### 6. Protein domain annotation
Structural/functional annotations (e.g. SH2 domains) are overlaid as a halo ring directly on the node, keeping biological context in the same figure without a separate panel.

---

## Usage

```r
# ── Settings ─────────────────────────────────────────────────
MERGE_MODE <- TRUE
ASSAY_CONFIG <- tribble(
  ~path,                              ~name,               ~convention, ~toupper_genes,
  "path/to/screen1_gene_summary.txt", "Screen1_CRISPRi",   "CRISPRi",   FALSE,
  "path/to/screen2_gene_summary.txt", "Screen2_CRISPRi",   "CRISPRi",   FALSE
  # mouse screen: set toupper_genes = TRUE
)

# ── Run ──────────────────────────────────────────────────────
source("CRISPR_merge_utils.R")
# then source the full CRISPR_to_ImmuneCell_pathway.R
```

### Adding a new screen
Uncomment one row in `ASSAY_CONFIG`. No other changes needed.

### Convention options

| `convention` | Use case |
|---|---|
| `"CRISPRi"` | Knockdown; depleted = GoF (default) |
| `"CRISPRko"` | Knockout; same sign convention as CRISPRi |
| `"CRISPRko_flip"` | Knockout; lab uses reversed terminology |
| `"CRISPRa"` | Activation; signs flipped automatically |

---

## Files

| File | Description |
|---|---|
| `CRISPR_to_ImmuneCell_pathway.R` | Main pipeline: layout, ggplot assembly, legend, ggsave |
| `CRISPR_merge_utils.R` | Utilities: `load_mageck()`, `merge_assays()`, `node_colours_merged()`, `build_dotmatrix()`, `assay_summary()`, `plot_overlap()` |
