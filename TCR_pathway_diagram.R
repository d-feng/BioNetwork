# ============================================================
# TCR Signaling Pathway Diagram — CRISPRi hit overlay
# Carnevale et al. MAGeCK gene_summary
#
# Approach : ggplot2 + ggforce (manual coordinates, typed edges)
# Node shapes: kinase=hexagon, TF=diamond, adaptor=ellipse,
#              receptor/complex=rounded rect
# Edge types : activation, phosphorylation, inhibition (blunt),
#              indirect activation (dashed), binding
# CRISPRi    : blue=GoF(depleted), red=LoF(enriched), grey=not hit
#              node size scaled by -log10(FDR)
# ============================================================

library(ggplot2)
library(ggforce)    # geom_regon, geom_ellipse, geom_shape
library(ggrepel)    # label placement
library(ggnewscale) # new_scale_fill() — two fill scales in one plot
library(dplyr)
library(tibble)
library(tidyr)
library(scales)     # rescale()
library(grid)       # arrow()

# ── 0. Settings ──────────────────────────────────────────────
INPUT_FILE <- "C:/Users/difen/Downloads/Carnevale_CRISPRi_Tcell_mageck_gene_summary.txt"
FDR_CUTOFF <- 0.1
LFC_CUTOFF <- 0.2

# ── 1. Load CRISPRi data ──────────────────────────────────────
raw <- read.delim(INPUT_FILE, check.names = FALSE, stringsAsFactors = FALSE)

crispr_hits <- raw %>%
  transmute(
    gene    = id,
    hit_dir = dplyr::case_when(
      neg.fdr < FDR_CUTOFF & neg.lfc < -LFC_CUTOFF ~ "depleted",
      pos.fdr < FDR_CUTOFF & pos.lfc >  LFC_CUTOFF ~ "enriched",
      TRUE                                          ~ "none"
    ),
    best_fdr = pmin(neg.fdr, pos.fdr),
    nlp_fdr  = -log10(pmax(pmin(neg.fdr, pos.fdr), 1e-6))
  )

# ── 2. Backbone node table (biologically positioned) ─────────
# Canvas: x ∈ [0,10], y ∈ [0,11]
# y: 10.5=extracellular → 9=plasma membrane → 2=cytoplasm → 0=nucleus
# x: NFAT(1) | MAPK(2-3) | TCR core(4-5) | PI3K/AKT(5-6) | NF-kB(7-8)
nodes_raw <- tribble(
  ~gene,    ~label,         ~x,    ~y,    ~type,       ~compartment,
  # ── Extracellular ──────────────────────────────────────────
  "pMHCII", "pMHC-II",      5.0,  10.8,  "ligand",    "extracellular",
  # ── Plasma membrane ────────────────────────────────────────
  "LAG3",   "LAG3",         1.5,   9.1,  "receptor",  "membrane",
  "PDCD1",  "PD-1",         3.0,   9.1,  "receptor",  "membrane",
  "CD247",  "CD3ζ",         4.2,   9.1,  "receptor",  "membrane",
  "TCR",    "TCR",          5.0,   9.1,  "complex",   "membrane",
  "CD4",    "CD4",          5.8,   9.1,  "receptor",  "membrane",
  "CD28",   "CD28",         7.0,   9.1,  "receptor",  "membrane",
  "CTLA4",  "CTLA4",        8.2,   9.1,  "receptor",  "membrane",
  # ── Proximal cytoplasm ──────────────────────────────────────
  "FYN",    "FYN",          4.0,   8.0,  "kinase",    "cytoplasm",
  "LCK",    "LCK",          5.8,   8.0,  "kinase",    "cytoplasm",
  "ZAP70",  "ZAP70",        5.0,   7.0,  "kinase",    "cytoplasm",
  "LAT",    "LAT",          4.2,   6.1,  "adaptor",   "membrane",
  "LCP2",   "SLP-76",       3.2,   5.4,  "adaptor",   "cytoplasm",
  "GRB2",   "GRB2",         5.2,   5.8,  "adaptor",   "cytoplasm",
  "ITK",    "ITK",          2.8,   5.8,  "kinase",    "cytoplasm",
  "VAV1",   "VAV1",         2.3,   5.1,  "adaptor",   "cytoplasm",
  # ── PLCγ / second messengers ───────────────────────────────
  "PLCG1",  "PLCγ1",        3.5,   5.1,  "kinase",    "cytoplasm",
  "RASGRP1","RasGRP1",      3.6,   4.4,  "adaptor",   "cytoplasm",
  # ── PI3K / AKT ─────────────────────────────────────────────
  "PIK3CD", "PI3Kδ",        5.2,   4.8,  "kinase",    "cytoplasm",
  "PIK3R1", "p85α",         6.5,   5.8,  "adaptor",   "endosome",
  "PDK1",   "PDK1",         5.2,   3.8,  "kinase",    "cytoplasm",
  "AKT1",   "AKT",          5.2,   3.0,  "kinase",    "cytoplasm",
  "MTOR",   "mTOR",         5.8,   2.2,  "complex",   "cytoplasm",
  # ── MAPK branch ────────────────────────────────────────────
  "KRAS",   "RAS",          3.5,   3.7,  "kinase",    "membrane",
  "RAF1",   "RAF1",         3.5,   3.0,  "kinase",    "cytoplasm",
  "MAP2K1", "MEK1/2",       3.5,   2.2,  "kinase",    "cytoplasm",
  "MAPK3",  "ERK1/2",       3.5,   1.4,  "kinase",    "cytoplasm",
  # ── PKCθ / NF-κB branch ─────────────────────────────────────
  "PRKCQ",  "PKCθ",         5.8,   4.2,  "kinase",    "cytoplasm",
  "CARD11", "CARMA1",       6.5,   3.5,  "complex",   "cytoplasm",
  "BCL10",  "BCL10",        7.1,   3.0,  "adaptor",   "cytoplasm",
  "MALT1",  "MALT1",        7.7,   3.0,  "adaptor",   "cytoplasm",
  "IKBKB",  "IKKβ",         7.1,   2.2,  "kinase",    "cytoplasm",
  "CHUK",   "IKKα",         7.8,   2.2,  "kinase",    "cytoplasm",
  "NFKBIA", "IκBα",         7.4,   1.4,  "adaptor",   "cytoplasm",
  # ── Calcineurin / NFAT ──────────────────────────────────────
  "PPP3CA", "Calcineurin",  1.8,   3.0,  "complex",   "cytoplasm",
  # ── Nucleus TFs ─────────────────────────────────────────────
  "NFATC1", "NFAT",         1.5,   0.5,  "tf",        "nucleus",
  "MAPK1",  "ERK",          3.0,   0.5,  "kinase",    "nucleus",
  "FOXO1",  "FOXO1",        4.5,   0.5,  "tf",        "nucleus",
  "RELA",   "p65 (RELA)",   6.5,   0.5,  "tf",        "nucleus",
  "NFKB1",  "p50 (NFKB1)", 7.5,   0.5,  "tf",        "nucleus"
)

# ── 3. Backbone edge table ────────────────────────────────────
edges_raw <- tribble(
  ~from,     ~to,       ~type,
  # pMHC → TCR
  "pMHCII",  "TCR",     "binding",
  # TCR complex
  "TCR",     "CD247",   "binding",
  "TCR",     "CD4",     "binding",
  # Proximal kinases
  "CD4",     "LCK",     "activation",
  "CD247",   "FYN",     "activation",
  "LCK",     "ZAP70",   "phosphorylation",
  "FYN",     "ZAP70",   "phosphorylation",
  # ZAP70 → adaptors
  "ZAP70",   "LAT",     "phosphorylation",
  "ZAP70",   "LCP2",    "phosphorylation",
  # LAT scaffold
  "LAT",     "PLCG1",   "activation",
  "LAT",     "GRB2",    "binding",
  "LAT",     "ITK",     "activation",
  # SLP-76 scaffold
  "LCP2",    "VAV1",    "activation",
  "LCP2",    "RASGRP1", "activation",
  "LCP2",    "ITK",     "binding",
  # ITK → PLCγ1
  "ITK",     "PLCG1",   "phosphorylation",
  # PLCγ1 branches
  "PLCG1",   "RASGRP1", "activation",
  "PLCG1",   "PRKCQ",   "activation",
  "PLCG1",   "PPP3CA",  "indirect_act",
  # MAPK cascade
  "RASGRP1", "KRAS",    "activation",
  "KRAS",    "RAF1",    "activation",
  "RAF1",    "MAP2K1",  "phosphorylation",
  "MAP2K1",  "MAPK3",   "phosphorylation",
  "MAPK3",   "MAPK1",   "activation",
  "MAPK1",   "MAPK1_n", "translocation",  # cytoplasm ERK → nuclear ERK
  # PI3K/AKT
  "CD28",    "PIK3CD",  "activation",
  "VAV1",    "PIK3CD",  "activation",
  "PIK3R1",  "PIK3CD",  "binding",
  "PIK3CD",  "PDK1",    "indirect_act",
  "PDK1",    "AKT1",    "phosphorylation",
  "AKT1",    "MTOR",    "activation",
  "AKT1",    "FOXO1",   "phosphorylation",
  # NF-κB
  "PRKCQ",   "CARD11",  "phosphorylation",
  "CARD11",  "BCL10",   "binding",
  "BCL10",   "MALT1",   "binding",
  "MALT1",   "IKBKB",   "activation",
  "IKBKB",   "NFKBIA",  "phosphorylation",
  "CHUK",    "NFKBIA",  "phosphorylation",
  "NFKBIA",  "RELA",    "inhibition",
  "NFKBIA",  "NFKB1",   "inhibition",
  # Calcineurin → NFAT
  "PPP3CA",  "NFATC1",  "activation",
  # Checkpoint inhibitors
  "CTLA4",   "LCK",     "inhibition",
  "PDCD1",   "ZAP70",   "inhibition",
  "LAG3",    "LCK",     "inhibition",
  # GRB2 → RAS (via SOS)
  "GRB2",    "KRAS",    "indirect_act",
  # Nuclear activation (terminal arrows)
  "NFATC1",  "gene_exp","activation",
  "MAPK1",   "gene_exp","activation",
  "FOXO1",   "gene_exp","inhibition",
  "RELA",    "gene_exp","activation",
  "NFKB1",   "gene_exp","activation"
)

# Add a "gene expression" terminal node
nodes_raw <- nodes_raw %>%
  dplyr::add_row(gene="gene_exp", label="Gene\nExpression",
                 x=4.5, y=-0.4, type="output", compartment="nucleus")

# ── 4. Merge CRISPRi data ─────────────────────────────────────
nodes <- nodes_raw %>%
  dplyr::left_join(crispr_hits, by = "gene") %>%
  dplyr::mutate(
    hit_dir  = tidyr::replace_na(hit_dir, "none"),
    nlp_fdr  = tidyr::replace_na(nlp_fdr, 0),
    fill_col = dplyr::case_when(
      hit_dir == "depleted" ~ "#2166AC",   # blue: essential/GoF
      hit_dir == "enriched" ~ "#D6604D",   # red:  suppressor/LoF
      TRUE                  ~ "#D9D9D9"    # grey: not significant
    ),
    border_col = dplyr::case_when(
      hit_dir == "depleted" ~ "#08306B",
      hit_dir == "enriched" ~ "#67000D",
      TRUE                  ~ "#888888"
    ),
    # Node radius in data units (hexagon circumradius / rect half-size)
    node_r = dplyr::case_when(
      hit_dir == "none" ~ 0.28,
      TRUE              ~ scales::rescale(pmin(nlp_fdr, 5),
                                          to    = c(0.32, 0.52),
                                          from  = c(0, 5))
    ),
    label_face = ifelse(hit_dir == "none", "plain", "bold"),
    label_col  = ifelse(hit_dir == "none", "grey40", "grey10")
  )

# ── 5. Build edge coordinate table ────────────────────────────
node_xy <- nodes %>% dplyr::select(gene, x, y)

edges <- edges_raw %>%
  dplyr::left_join(node_xy, by = c("from" = "gene")) %>%
  dplyr::rename(x1 = x, y1 = y) %>%
  dplyr::left_join(node_xy, by = c("to" = "gene")) %>%
  dplyr::rename(x2 = x, y2 = y) %>%
  dplyr::filter(!is.na(x1), !is.na(x2)) %>%
  # Shorten edges so arrowheads don't overlap node centres
  dplyr::mutate(
    dx   = x2 - x1,
    dy   = y2 - y1,
    dist = sqrt(dx^2 + dy^2),
    # Trim 0.32 data units from each end (≈ node radius)
    trim = 0.30,
    x1s  = x1 + dx / dist * trim,
    y1s  = y1 + dy / dist * trim,
    x2s  = x2 - dx / dist * trim,
    y2s  = y2 - dy / dist * trim
  )

# Inhibition blunt-end perpendicular tick
inh_edges <- edges %>%
  dplyr::filter(type == "inhibition") %>%
  dplyr::mutate(
    px  = -dy / dist * 0.20,
    py  =  dx / dist * 0.20,
    bx1 = x2s - px, by1 = y2s - py,
    bx2 = x2s + px, by2 = y2s + py
  )

# ── 6. Compartment boxes ──────────────────────────────────────
compartments <- tribble(
  ~label,            ~xmin,  ~xmax,  ~ymin,  ~ymax,   ~fill,
  "Extracellular",    0.5,    9.5,   10.3,   11.2,   "#EBF5FB",
  "Plasma membrane",  0.5,    9.5,    8.5,   10.3,   "#FEF9E7",
  "Cytoplasm",        0.5,    9.5,    0.9,    8.5,   "#EAFAF1",
  "Endosome",         6.0,    9.0,    5.2,    6.5,   "#F5EEF8",
  "Nucleus",          0.5,    9.5,   -0.7,    0.9,   "#EBF5FB"
)

# ── 7. Arrow and colour specs ─────────────────────────────────
arr       <- arrow(length = unit(0.10, "inches"), type = "closed")
arr_open  <- arrow(length = unit(0.09, "inches"), type = "open")
arr_small <- arrow(length = unit(0.07, "inches"), type = "closed")

col_act   <- "grey25"
col_phos  <- "#C07020"
col_inh   <- "#C0392B"
col_ind   <- "grey55"
col_bind  <- "#5B7FA6"
col_trans <- "#8E44AD"

# ── 8. Node shape helper functions ───────────────────────────
# Rounded rectangle polygon (4 corners) for geom_shape
make_rect_poly <- function(df, hw = NULL, hh = 0.22) {
  df %>%
    dplyr::rowwise() %>%
    dplyr::reframe(
      gene      = gene,
      fill_col  = fill_col,
      border_col= border_col,
      label     = label,
      label_face= label_face,
      label_col = label_col,
      node_r    = node_r,
      px = c(x - (if(is.null(hw)) node_r * 1.3 else hw),
             x + (if(is.null(hw)) node_r * 1.3 else hw),
             x + (if(is.null(hw)) node_r * 1.3 else hw),
             x - (if(is.null(hw)) node_r * 1.3 else hw)),
      py = c(y - hh, y - hh, y + hh, y + hh)
    )
}

receptor_poly <- make_rect_poly(dplyr::filter(nodes, type %in% c("receptor","ligand","output")))
complex_poly  <- make_rect_poly(dplyr::filter(nodes, type == "complex"))

# ── 9. Pathway group annotations (text boxes) ─────────────────
pathway_labels <- tribble(
  ~x,   ~y,    ~label,
  1.6,   6.5,  "Ca²⁺ / NFAT",
  2.8,   4.1,  "RAS–MAPK",
  5.0,   2.6,  "PI3K–AKT–mTOR",
  7.0,   4.0,  "NF-κB",
  4.8,   7.5,  "TCR proximal\nsignaling"
)

# ── 10. Assemble plot ─────────────────────────────────────────
p <- ggplot() +

  # ── Compartment background boxes ──
  # fill = identity (pre-specified hex colour strings)
  geom_rect(
    data = compartments,
    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill),
    colour = "grey70", linewidth = 0.4, alpha = 0.7
  ) +
  scale_fill_identity() +   # compartment rect fill: use colour strings directly

  # Switch to a fresh fill scale for all subsequent node layers
  new_scale_fill() +

  geom_text(
    data = compartments,
    aes(x = xmin + 0.12, y = ymax - 0.07, label = label),
    hjust = 0, vjust = 1, size = 3.0, colour = "grey45",
    fontface = "bold.italic"
  ) +

  # ── Pathway group labels ──
  geom_label(
    data = pathway_labels,
    aes(x = x, y = y, label = label),
    fill = alpha("white", 0.6), colour = "grey30",
    size = 2.6, fontface = "bold.italic",
    linewidth = 0.2, label.r = unit(0.15, "lines")
  ) +

  # ════ EDGES (drawn before nodes) ════════════════════════════

  # Activation — solid dark arrow
  geom_segment(
    data    = dplyr::filter(edges, type == "activation"),
    aes(x=x1s, y=y1s, xend=x2s, yend=y2s),
    arrow     = arr,
    colour    = col_act,
    linewidth = 0.55,
    lineend   = "round"
  ) +

  # Phosphorylation — solid orange arrow
  geom_segment(
    data    = dplyr::filter(edges, type == "phosphorylation"),
    aes(x=x1s, y=y1s, xend=x2s, yend=y2s),
    arrow     = arr,
    colour    = col_phos,
    linewidth = 0.65,
    lineend   = "round"
  ) +

  # Indirect activation — dashed grey arrow
  geom_segment(
    data    = dplyr::filter(edges, type == "indirect_act"),
    aes(x=x1s, y=y1s, xend=x2s, yend=y2s),
    arrow    = arr,
    colour   = col_ind,
    linewidth = 0.45,
    linetype  = "dashed"
  ) +

  # Binding — plain line, no arrowhead
  geom_segment(
    data    = dplyr::filter(edges, type == "binding"),
    aes(x=x1s, y=y1s, xend=x2s, yend=y2s),
    colour    = col_bind,
    linewidth = 0.45
  ) +

  # Translocation — purple dashed arrow
  geom_curve(
    data    = dplyr::filter(edges, type == "translocation"),
    aes(x=x1s, y=y1s, xend=x2s, yend=y2s),
    arrow      = arr_small,
    colour     = col_trans,
    linewidth  = 0.45,
    linetype   = "dashed",
    curvature  = 0.3
  ) +

  # Inhibition — red line (no arrowhead; blunt tick added below)
  geom_segment(
    data    = dplyr::filter(edges, type == "inhibition"),
    aes(x=x1s, y=y1s, xend=x2s, yend=y2s),
    colour    = col_inh,
    linewidth = 0.55,
    lineend   = "butt"
  ) +
  # Blunt-end perpendicular tick
  geom_segment(
    data      = inh_edges,
    aes(x=bx1, y=by1, xend=bx2, yend=by2),
    colour    = col_inh,
    linewidth = 1.0
  ) +

  # ════ CHECKPOINT inhibition curves (curved so they don't cross) ════
  geom_curve(
    data      = dplyr::filter(edges, type == "inhibition",
                              from %in% c("CTLA4","LAG3")),
    aes(x=x1s, y=y1s, xend=x2s, yend=y2s),
    colour    = col_inh,
    linewidth = 0.55,
    curvature = 0.35,
    lineend   = "butt"
  ) +

  # ════ NODES ══════════════════════════════════════════════════

  # Receptors & ligands — rounded rectangle
  geom_shape(
    data    = receptor_poly,
    aes(x=px, y=py, group=gene, fill=fill_col, colour=border_col),
    radius  = unit(4, "pt"),
    linewidth = 0.45
  ) +

  # Complexes — rounded rectangle, dashed border
  geom_shape(
    data     = complex_poly,
    aes(x=px, y=py, group=gene, fill=fill_col, colour=border_col),
    radius   = unit(4, "pt"),
    linewidth = 0.5,
    linetype  = "dashed"
  ) +

  # Kinases — hexagons (sides is an aesthetic in ggforce 0.5+)
  geom_regon(
    data    = dplyr::filter(nodes, type == "kinase"),
    aes(x0=x, y0=y, r=node_r * 0.90, angle=pi/6,
        sides=6, fill=fill_col, colour=border_col),
    linewidth = 0.45
  ) +

  # Transcription factors — diamonds (4-sided rotated 45°)
  geom_regon(
    data    = dplyr::filter(nodes, type == "tf"),
    aes(x0=x, y0=y, r=node_r * 1.05, angle=pi/4,
        sides=4, fill=fill_col, colour=border_col),
    linewidth = 0.45
  ) +

  # Adaptors — ellipses
  geom_ellipse(
    data    = dplyr::filter(nodes, type == "adaptor"),
    aes(x0=x, y0=y, a=node_r * 1.30, b=node_r * 0.65,
        angle=0, fill=fill_col, colour=border_col),
    linewidth = 0.45
  ) +

  # All node layers use pre-computed fill_col strings → identity scale
  scale_fill_identity() +
  scale_colour_identity() +

  # ════ LABELS ══════════════════════════════════════════════════
  geom_label_repel(
    data         = dplyr::filter(nodes, type != "output"),
    aes(x=x, y=y, label=label,
        fontface = label_face,
        colour   = label_col),
    fill          = alpha("white", 0.80),
    size          = 2.4,
    label.size    = 0.15,              # border thickness (ggrepel uses label.size, not linewidth)
    label.r       = unit(0.12, "lines"),
    box.padding   = 0.18,
    point.padding = 0.05,
    segment.colour= "grey60",
    segment.size  = 0.25,
    seed          = 99,
    max.overlaps  = Inf,
    show.legend   = FALSE
  ) +

  # Gene expression output label
  geom_label(
    data = dplyr::filter(nodes, type == "output"),
    aes(x=x, y=y, label=label),
    fill = "#F8C471", colour = "grey20",
    size = 3, fontface = "bold",
    linewidth = 0.4,                   # replaces deprecated label.size
    label.r = unit(0.2, "lines")
  ) +

  # ════ LEGEND (manual) ════════════════════════════════════════
  # CRISPRi hit legend — small coloured shapes bottom-right
  annotate("rect",   xmin=8.8, xmax=9.3, ymin=-0.65, ymax=-0.35,
           fill="#2166AC", colour="#08306B", linewidth=0.4) +
  annotate("text",   x=9.4, y=-0.50,
           label="Depleted (GoF, essential)", hjust=0, size=2.4) +
  annotate("rect",   xmin=8.8, xmax=9.3, ymin=-1.05, ymax=-0.75,
           fill="#D6604D", colour="#67000D", linewidth=0.4) +
  annotate("text",   x=9.4, y=-0.90,
           label="Enriched (LoF, suppressor)", hjust=0, size=2.4) +
  annotate("rect",   xmin=8.8, xmax=9.3, ymin=-1.45, ymax=-1.15,
           fill="#D9D9D9", colour="#888888", linewidth=0.4) +
  annotate("text",   x=9.4, y=-1.30,
           label="Not significant", hjust=0, size=2.4) +
  # Edge type legend
  annotate("segment", x=0.6, xend=1.3, y=-0.50, yend=-0.50,
           colour=col_act,  linewidth=0.7,
           arrow=arrow(length=unit(0.07,"in"), type="closed")) +
  annotate("text", x=1.4, y=-0.50, label="Activation", hjust=0, size=2.4) +
  annotate("segment", x=0.6, xend=1.3, y=-0.85, yend=-0.85,
           colour=col_phos, linewidth=0.7,
           arrow=arrow(length=unit(0.07,"in"), type="closed")) +
  annotate("text", x=1.4, y=-0.85, label="Phosphorylation", hjust=0, size=2.4) +
  annotate("segment", x=0.6, xend=1.3, y=-1.20, yend=-1.20,
           colour=col_inh,  linewidth=0.7) +
  annotate("segment", x=1.3, xend=1.3, y=-1.35, yend=-1.05,
           colour=col_inh,  linewidth=1.0) +
  annotate("text", x=1.4, y=-1.20, label="Inhibition", hjust=0, size=2.4) +
  annotate("segment", x=3.2, xend=3.9, y=-0.50, yend=-0.50,
           colour=col_ind,  linewidth=0.7, linetype="dashed",
           arrow=arrow(length=unit(0.07,"in"), type="closed")) +
  annotate("text", x=4.0, y=-0.50, label="Indirect activation", hjust=0, size=2.4) +
  annotate("segment", x=3.2, xend=3.9, y=-0.85, yend=-0.85,
           colour=col_bind, linewidth=0.7) +
  annotate("text", x=4.0, y=-0.85, label="Binding", hjust=0, size=2.4) +
  annotate("segment", x=3.2, xend=3.9, y=-1.20, yend=-1.20,
           colour=col_trans, linewidth=0.7, linetype="dashed",
           arrow=arrow(length=unit(0.07,"in"), type="closed")) +
  annotate("text", x=4.0, y=-1.20, label="Translocation", hjust=0, size=2.4) +
  # Shape legend
  annotate("text", x=6.0, y=-0.50, label="■ = Receptor/Complex  ⬡ = Kinase",
           hjust=0, size=2.4, colour="grey40") +
  annotate("text", x=6.0, y=-0.85, label="◆ = Transcription factor  ○ = Adaptor",
           hjust=0, size=2.4, colour="grey40") +
  annotate("text", x=6.0, y=-1.20, label="Size ∝ −log₁₀(FDR) of CRISPRi hit",
           hjust=0, size=2.4, colour="grey40") +

  # ════ SCALES & THEME ═════════════════════════════════════════
  # Note: scale_fill_identity() and scale_colour_identity() are already
  # declared inline above (after node layers) — do not repeat here.
  coord_fixed(ratio = 1, xlim = c(0.3, 13.5), ylim = c(-1.6, 11.4)) +
  theme_void(base_size = 11) +
  theme(
    plot.margin    = margin(8, 8, 8, 8),
    plot.title     = element_text(face = "bold", size = 14, hjust = 0.5,
                                  margin = margin(b = 4)),
    plot.subtitle  = element_text(size  = 8.5,  hjust = 0.5, colour = "grey40",
                                  margin = margin(b = 6))
  ) +
  labs(
    title    = "TCR Signaling Pathway — Carnevale et al. CRISPRi Screen",
    subtitle = paste0(
      "Node fill: blue = GoF/depleted (FDR<", FDR_CUTOFF, ", LFC<-", LFC_CUTOFF, ")   ",
      "red = LoF/enriched (FDR<", FDR_CUTOFF, ", LFC>", LFC_CUTOFF, ")   ",
      "grey = not significant\n",
      "Node size ∝ −log₁₀(FDR)   Shape: hexagon=kinase, diamond=TF, oval=adaptor, rect=receptor/complex"
    )
  )

print(p)

ggsave("TCR_pathway_CRISPRi.pdf", p, width = 17, height = 14, units = "in")
ggsave("TCR_pathway_CRISPRi.png", p, width = 17, height = 14, units = "in", dpi = 300)
cat("Saved: TCR_pathway_CRISPRi.pdf / .png\n")

# ── Summary of hits on the backbone ──────────────────────────
cat("\n── CRISPRi hits mapped onto backbone ───────────────────\n")
nodes %>%
  dplyr::filter(hit_dir != "none", type != "output") %>%
  dplyr::select(gene, label, hit_dir, nlp_fdr, best_fdr, compartment) %>%
  dplyr::arrange(dplyr::desc(nlp_fdr)) %>%
  print(n = 30)
