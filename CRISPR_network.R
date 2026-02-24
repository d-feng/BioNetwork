# ============================================================
# CRISPRi Network Analysis — Carnevale et al. T Cell Screen
# Input : Carnevale_CRISPRi_Tcell_mageck_gene_summary.txt
# New   : pathway grouping (clusterProfiler/ReactomePA),
#         hub / master-regulator scoring, annotated network plot
# ============================================================

# ── Required packages ────────────────────────────────────────
# Install if missing:
# install.packages(c("ggforce", "patchwork", "ggnewscale"))
# BiocManager::install(c("STRINGdb","clusterProfiler","ReactomePA","org.Hs.eg.db"))

# Load Bioconductor packages FIRST (they mask some dplyr verbs),
# then tidyverse LAST so dplyr wins back select/filter/rename etc.
library(STRINGdb)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)   # gene symbol ↔ Entrez mapping

library(tidyverse)      # MUST come after Bioc packages → dplyr::select wins
library(igraph)
library(ggraph)
library(ggrepel)
library(ggforce)        # geom_mark_hull (pathway hulls)
library(patchwork)      # combine plots
library(ggnewscale)     # new_scale_colour() — allows two colour scales in one plot

# Belt-and-suspenders: pin every tidyverse verb that Bioc packages can mask.
# AnnotationDbi / clusterProfiler / org.Hs.eg.db overwrite many dplyr & tidyr names.
select       <- dplyr::select
filter       <- dplyr::filter
rename       <- dplyr::rename
mutate       <- dplyr::mutate
transmute    <- dplyr::transmute
arrange      <- dplyr::arrange
summarise    <- dplyr::summarise
summarize    <- dplyr::summarize
group_by     <- dplyr::group_by
ungroup      <- dplyr::ungroup
slice_head   <- dplyr::slice_head
distinct     <- dplyr::distinct
left_join    <- dplyr::left_join
bind_rows    <- dplyr::bind_rows
pull         <- dplyr::pull
count        <- dplyr::count
desc         <- dplyr::desc
tibble       <- tibble::tibble
separate_rows <- tidyr::separate_rows
replace_na   <- tidyr::replace_na
fct_reorder  <- forcats::fct_reorder
map_dbl      <- purrr::map_dbl
str_trunc    <- stringr::str_trunc

# ── 0. Settings ──────────────────────────────────────────────
INPUT_FILE       <- "C:/Users/difen/Downloads/Carnevale_CRISPRi_Tcell_mageck_gene_summary.txt"
SPECIES          <- 9606       # Homo sapiens
STRING_SCORE     <- 400        # 400 = medium confidence, 700 = high
TOP_N_GENES      <- 80         # wider net so pathway enrichment has power
HIT_MODE         <- "both"     # "neg" | "pos" | "both"
CELL_TYPE        <- "CRISPRi T cells (Carnevale et al.)"
EXPORT_CYTOSCAPE <- FALSE      # set TRUE if Cytoscape is running locally

# Pathway annotation source: "Reactome" or "KEGG"
PATHWAY_SOURCE   <- "Reactome"
N_TOP_PATHWAYS   <- 6          # how many pathway hulls to draw on the network
PATHWAY_FDR      <- 0.05       # ORA FDR cutoff

# Hub gene thresholds (used in regulator table)
HUB_TOP_N        <- 15         # highlight top N hub genes in network

# ── 1. Load MAGeCK output ────────────────────────────────────
raw <- read_tsv(INPUT_FILE, show_col_types = FALSE)
cat(sprintf("Loaded %d genes from %s\n", nrow(raw), basename(INPUT_FILE)))

# ── 2. Select hits ────────────────────────────────────────────
select_hits <- function(data, mode, top_n) {
  half <- ceiling(top_n / 2)

  neg_hits <- data %>%
    filter(neg.fdr < 0.1, neg.lfc < -0.3) %>%
    arrange(neg.fdr, neg.lfc) %>%
    slice_head(n = if (mode == "both") half else top_n) %>%
    transmute(
      Gene       = id,
      LFC        = neg.lfc,
      FDR        = neg.fdr,
      score      = neg.score,
      good_sgrna = neg.goodsgrna,
      direction  = "Depleted (gain of function)"   # knockdown hurts → gene drives T cells
    )

  pos_hits <- data %>%
    filter(pos.fdr < 0.1, pos.lfc > 0.3) %>%
    arrange(pos.fdr, desc(pos.lfc)) %>%
    slice_head(n = if (mode == "both") half else top_n) %>%
    transmute(
      Gene       = id,
      LFC        = pos.lfc,
      FDR        = pos.fdr,
      score      = pos.score,
      good_sgrna = pos.goodsgrna,
      direction  = "Enriched (loss of function)"   # knockdown helps → gene suppresses T cells
    )

  switch(mode,
    neg  = neg_hits,
    pos  = pos_hits,
    both = bind_rows(neg_hits, pos_hits) %>% distinct(Gene, .keep_all = TRUE)
  )
}

hits <- select_hits(raw, HIT_MODE, TOP_N_GENES)
cat(sprintf("Hits selected: %d genes (mode = '%s')\n", nrow(hits), HIT_MODE))
if (nrow(hits) < 3) stop("Too few hits — relax thresholds or change HIT_MODE.")

# ── 3. STRINGdb query ─────────────────────────────────────────
string_db <- STRINGdb$new(
  version         = "12.0",
  species         = SPECIES,
  score_threshold = STRING_SCORE,
  input_directory = ""
)

hits_mapped_all <- string_db$map(
  my_data_frame              = as.data.frame(hits),
  my_data_frame_id_col_names = "Gene",
  removeUnmappedRows         = FALSE   # keep unmapped so we can report them
)

# ── Categorise every hit gene ────────────────────────────────
genes_unmapped <- hits_mapped_all %>%
  filter(is.na(STRING_id)) %>%
  dplyr::pull(Gene)

# Keep only mapped rows; deduplicate ambiguous mappings → best match per gene
hits_mapped <- hits_mapped_all %>%
  filter(!is.na(STRING_id)) %>%
  dplyr::group_by(Gene) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

cat(sprintf("Mapped to STRING : %d / %d genes\n", nrow(hits_mapped), nrow(hits)))
if (length(genes_unmapped) > 0)
  cat(sprintf("NOT mapped (lost): %d genes — %s\n",
              length(genes_unmapped), paste(genes_unmapped, collapse = ", ")))

interactions <- string_db$get_interactions(hits_mapped$STRING_id)
cat(sprintf("STRING interactions: %d\n", nrow(interactions)))
if (nrow(interactions) == 0) stop("No interactions. Lower STRING_SCORE or expand TOP_N_GENES.")

# ── 4. Build igraph ───────────────────────────────────────────
id2gene <- setNames(hits_mapped$Gene, hits_mapped$STRING_id)

edges <- interactions %>%
  transmute(
    from   = id2gene[from],
    to     = id2gene[to],
    weight = combined_score / 1000
  ) %>%
  filter(!is.na(from), !is.na(to)) %>%
  filter(from != to) %>%
  dplyr::distinct(from, to, .keep_all = TRUE)

# ALL mapped genes as vertices (including those with no interactions)
nodes <- hits_mapped %>%
  dplyr::select(Gene, LFC, FDR, score, good_sgrna, direction) %>%
  dplyr::distinct(Gene, .keep_all = TRUE) %>%
  mutate(neg_log10_fdr = -log10(pmax(FDR, 1e-10)))

g <- graph_from_data_frame(d = edges, directed = FALSE, vertices = nodes)

# Identify isolated nodes (mapped but no interactions with other hits)
isolated_genes <- V(g)$name[degree(g) == 0]
if (length(isolated_genes) > 0)
  cat(sprintf("Mapped but isolated (no PPI edges): %d genes — %s\n",
              length(isolated_genes), paste(isolated_genes, collapse = ", ")))

# ── Keep isolated nodes in graph but flag them ───────────────
# Do NOT delete them — they are still significant CRISPR hits,
# just not connected to other hits in STRING at this threshold.
# They appear as standalone nodes on the plot, labelled in italics.
V(g)$is_isolated <- degree(g) == 0

cat(sprintf(
  "\nGene accounting:\n  Selected: %d  |  Unmapped: %d  |  Mapped+connected: %d  |  Mapped+isolated: %d\n",
  nrow(hits), length(genes_unmapped),
  sum(degree(g) > 0), sum(degree(g) == 0)
))

# ── 5. Network centrality → hub / regulator scoring ──────────
# Centrality is only meaningful for connected nodes.
# Isolated nodes receive 0 for all metrics and are excluded from hub ranking.
V(g)$degree      <- degree(g)
V(g)$betweenness <- betweenness(g, weights = 1 / E(g)$weight, normalized = TRUE)
V(g)$eigenvector <- eigen_centrality(g, weights = E(g)$weight, scale = TRUE)$vector
V(g)$pagerank    <- page_rank(g, weights = E(g)$weight)$vector

# Z-score each metric and average → composite hub score
zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
V(g)$hub_score <- rowMeans(cbind(
  zscore(V(g)$degree),
  zscore(V(g)$betweenness),
  zscore(V(g)$eigenvector),
  zscore(V(g)$pagerank)
), na.rm = TRUE)

# Only rank connected nodes as hubs (isolated nodes score 0 by definition)
hub_rank <- ifelse(V(g)$is_isolated, Inf, rank(-V(g)$hub_score))
V(g)$is_hub <- hub_rank <= HUB_TOP_N

cat(sprintf("Final network: %d nodes, %d edges\n", vcount(g), ecount(g)))

# ── 6. Pathway enrichment (ORA) ───────────────────────────────
# Convert gene symbols → Entrez IDs
gene_list <- V(g)$name

entrez_map <- bitr(
  geneID   = gene_list,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
cat(sprintf("Entrez mapped: %d / %d genes\n", nrow(entrez_map), length(gene_list)))

run_ora <- function(entrez_ids, source = "Reactome") {
  if (source == "Reactome") {
    enrichPathway(
      gene          = entrez_ids,
      organism      = "human",
      pvalueCutoff  = PATHWAY_FDR,
      qvalueCutoff  = 0.2,
      readable      = TRUE   # show gene symbols in result
    )
  } else {
    enrichKEGG(
      gene          = entrez_ids,
      organism      = "hsa",
      pvalueCutoff  = PATHWAY_FDR,
      qvalueCutoff  = 0.2
    )
  }
}

ora_result <- run_ora(entrez_map$ENTREZID, PATHWAY_SOURCE)

if (is.null(ora_result) || nrow(as.data.frame(ora_result)) == 0) {
  warning("No enriched pathways found — pathway hull annotation will be skipped.")
  pathway_df   <- tibble()
  gene2pathway <- tibble(Gene = character(), pathway_label = character())
} else {
  pathway_df <- as.data.frame(ora_result) %>%
    arrange(p.adjust) %>%
    slice_head(n = N_TOP_PATHWAYS) %>%
    mutate(
      # Shorten long pathway names for the legend
      pathway_label = str_trunc(Description, 45),
      pathway_rank  = row_number()
    )

  cat(sprintf("Top %d enriched pathways:\n", nrow(pathway_df)))
  pathway_df %>% dplyr::select(pathway_label, p.adjust, Count) %>% print()

  # Map each gene to its best (lowest p.adjust) pathway
  gene2pathway <- pathway_df %>%
    dplyr::select(pathway_label, geneID) %>%
    separate_rows(geneID, sep = "/") %>%
    dplyr::rename(Gene = geneID) %>%
    group_by(Gene) %>%
    slice_head(n = 1) %>%          # one pathway per gene
    ungroup()
}

# Attach pathway label to graph vertices
node_df <- tibble(Gene = V(g)$name) %>%
  left_join(gene2pathway, by = "Gene") %>%
  replace_na(list(pathway_label = "Other"))

V(g)$pathway <- node_df$pathway_label

# ── 7. Subcellular compartment annotation (GO CC) ────────────
# Uses GOALL which traverses the GO DAG — a gene annotated to
# "mitochondrial inner membrane" will correctly match "mitochondria"
compartment_go <- list(
  "Cell surface"   = c("GO:0009986", "GO:0005886"),  # cell surface, plasma membrane
  "Nucleus"        = c("GO:0005634", "GO:0000228"),  # nucleus, nuclear chromosome
  "Mitochondria"   = c("GO:0005739", "GO:0005743"),  # mitochondrion, inner membrane
  "Endosome"       = c("GO:0005768", "GO:0010008"),  # endosome, early endosome
  "Lysosome"       = c("GO:0005764"),
  "ER"             = c("GO:0005783", "GO:0005788"),  # ER, ER lumen
  "Golgi"          = c("GO:0005794", "GO:0005795"),
  "Cytoplasm"      = c("GO:0005737", "GO:0005829"),  # cytoplasm, cytosol
  "Cytoskeleton"   = c("GO:0005856"),
  "Extracellular"  = c("GO:0005576", "GO:0005615")   # extracellular region, space
)

# Priority order: if a gene fits multiple, first match wins
compartment_priority <- names(compartment_go)

# Query GO CC for all network genes using GOALL (ancestor-inclusive)
gene_symbols_net <- V(g)$name
entrez_cc <- tryCatch(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = gene_symbols_net,
    columns = c("SYMBOL", "GOALL", "ONTOLOGYALL"),
    keytype = "SYMBOL"
  ) %>%
    filter(ONTOLOGYALL == "CC") %>%
    dplyr::select(SYMBOL, GOALL) %>%
    dplyr::distinct(),
  error = function(e) {
    warning("GO CC query failed: ", e$message)
    tibble(SYMBOL = character(), GOALL = character())
  }
)

# Map each gene to its highest-priority compartment
assign_compartment <- function(symbol, go_df, priority, go_map) {
  gene_go <- go_df$GOALL[go_df$SYMBOL == symbol]
  for (comp in priority) {
    if (any(gene_go %in% go_map[[comp]])) return(comp)
  }
  "Other"
}

compartment_vec <- sapply(gene_symbols_net, assign_compartment,
                           go_df      = entrez_cc,
                           priority   = compartment_priority,
                           go_map     = compartment_go)
V(g)$compartment <- compartment_vec

cat(sprintf("Compartment annotation: %d / %d genes annotated\n",
            sum(V(g)$compartment != "Other"), vcount(g)))
print(sort(table(V(g)$compartment), decreasing = TRUE))

# ── 7b. Protein complex annotation (CORUM) ───────────────────
# Downloads humanComplexes.txt once and caches locally
corum_file <- "humanComplexes.txt"
corum_url  <- "https://mips.helmholtz-muenchen.de/corum/download/releases/current/humanComplexes.txt.zip"

if (!file.exists(corum_file)) {
  message("Downloading CORUM database (~2 MB)...")
  tmp <- tempfile(fileext = ".zip")
  tryCatch({
    download.file(corum_url, tmp, quiet = TRUE)
    unzip(tmp, exdir = ".")
    message("CORUM downloaded and saved to: ", corum_file)
  }, error = function(e) {
    warning("CORUM download failed: ", e$message,
            "\nComplex annotation will be skipped.")
    corum_file <<- NULL
  })
}

if (!is.null(corum_file) && file.exists(corum_file)) {
  corum_raw <- read.delim(corum_file, sep = "\t", stringsAsFactors = FALSE,
                          check.names = FALSE)

  # Column may be named "subunits(Gene name)" — normalise
  subunit_col <- grep("subunits.*[Gg]ene", colnames(corum_raw), value = TRUE)[1]

  corum_long <- do.call(rbind, lapply(seq_len(nrow(corum_raw)), function(i) {
    genes <- trimws(strsplit(corum_raw[[subunit_col]][i], ";")[[1]])
    genes <- genes[nzchar(genes)]
    if (!length(genes)) return(NULL)
    data.frame(Gene         = genes,
               complex_name = corum_raw$ComplexName[i],
               stringsAsFactors = FALSE)
  }))

  # Collapse all complexes per gene into one string
  complex_annot <- corum_long %>%
    filter(Gene %in% gene_symbols_net) %>%
    group_by(Gene) %>%
    summarise(complexes = paste(unique(complex_name), collapse = " | "),
              n_complexes = n_distinct(complex_name),
              .groups = "drop")

  # Attach to graph
  complex_df <- tibble(Gene = gene_symbols_net) %>%
    left_join(complex_annot, by = "Gene") %>%
    replace_na(list(complexes = "none", n_complexes = 0L))

  V(g)$complexes   <- complex_df$complexes
  V(g)$n_complexes <- complex_df$n_complexes

  cat(sprintf("Complex annotation: %d / %d genes in ≥1 CORUM complex\n",
              sum(V(g)$n_complexes > 0), vcount(g)))
} else {
  V(g)$complexes   <- "none"
  V(g)$n_complexes <- 0L
}

# ── 8. Community detection ────────────────────────────────────
set.seed(123)
communities   <- cluster_louvain(g, weights = E(g)$weight)
V(g)$community <- as.character(membership(communities))
n_communities  <- max(membership(communities))

# ── 8. Network plot with pathway hulls ───────────────────────
# Hull border colour palette (discrete) — kept separate from node fill (continuous)
pathway_levels   <- unique(V(g)$pathway)
named_pathways   <- setdiff(pathway_levels, "Other")
n_named          <- length(named_pathways)

hull_border_pal <- setNames(
  scales::hue_pal()(max(n_named, 1)),
  named_pathways
)

# ── Build layout: FR for connected nodes, neat grid for isolated ─
set.seed(42)
layout_coords <- create_layout(g, layout = "fr")
layout_coords$pathway      <- V(g)$pathway
layout_coords$is_hub       <- V(g)$is_hub
layout_coords$hub_score    <- V(g)$hub_score
layout_coords$is_isolated  <- V(g)$is_isolated
layout_coords$compartment  <- V(g)$compartment
layout_coords$complexes    <- V(g)$complexes
layout_coords$n_complexes  <- V(g)$n_complexes

# Rescale connected-node layout to a fixed x range [0, 10]
connected_idx <- which(!layout_coords$is_isolated)
isolated_idx  <- which(layout_coords$is_isolated)

if (length(connected_idx) > 0) {
  cx <- layout_coords$x[connected_idx]
  cy <- layout_coords$y[connected_idx]
  # Normalise to [0, 10] × [0, 10]
  layout_coords$x[connected_idx] <- (cx - min(cx)) / (max(cx) - min(cx) + 1e-9) * 10
  layout_coords$y[connected_idx] <- (cy - min(cy)) / (max(cy) - min(cy) + 1e-9) * 10
}

# Place isolated nodes in a tidy column to the right (x = 11.8+)
# sorted by LFC so direction is visually grouped
if (length(isolated_idx) > 0) {
  iso_lfc   <- layout_coords$LFC[isolated_idx]
  iso_order <- order(iso_lfc)                       # sort by LFC (neg → pos)
  n_iso     <- length(isolated_idx)
  cols      <- max(1, ceiling(n_iso / 20))          # wrap into columns of 20
  rows_per_col <- ceiling(n_iso / cols)

  col_idx <- ((seq_along(iso_order) - 1) %/% rows_per_col)
  row_idx <- ((seq_along(iso_order) - 1) %%  rows_per_col)

  layout_coords$x[isolated_idx[iso_order]] <- 12.0 + col_idx * 1.4
  layout_coords$y[isolated_idx[iso_order]] <- 10 - row_idx * (10 / rows_per_col)
}

# x position of the dividing line between network and isolated strip
divider_x <- 11.3

# Hull data: only annotated genes (exclude "Other") and not isolated
hull_data <- layout_coords %>% filter(pathway != "Other", !is_isolated)

p_network <- ggraph(layout_coords) +

  # ── Divider between main network and isolated-gene strip ──
  geom_vline(
    xintercept = divider_x,
    linetype   = "dashed",
    colour     = "grey70",
    linewidth  = 0.4
  ) +
  annotate(
    "text",
    x      = divider_x + 0.15,
    y      = 10.3,
    label  = "No PPI partners\namong hits",
    hjust  = 0,
    vjust  = 1,
    size   = 2.5,
    colour = "grey50",
    fontface = "italic"
  ) +

  # ── Pathway hulls: use 'colour' aesthetic (border) to avoid fill conflict ──
  # fill is set to a fixed semi-transparent white so only the border carries colour
  geom_mark_hull(
    data         = hull_data,
    aes(x = x, y = y, colour = pathway, label = pathway),
    fill         = alpha("white", 0.08),   # fixed light wash — not mapped
    linewidth    = 0.8,
    label.fill   = alpha("white", 0.75),
    label.colour = "grey15",
    label.size   = 2.5,
    expand       = unit(6, "mm"),
    show.legend  = TRUE
  ) +
  scale_colour_manual(
    values = hull_border_pal,
    name   = paste(PATHWAY_SOURCE, "pathway")
  ) +

  # Allow a second independent colour scale for node borders
  new_scale_colour() +

  # ── Edges ──
  geom_edge_link(
    aes(width = weight, alpha = weight),
    colour      = "grey55",
    show.legend = FALSE
  ) +
  scale_edge_width(range = c(0.3, 2.2)) +
  scale_edge_alpha(range = c(0.12, 0.70)) +

  # ── Compartment colour palette for node borders ──────────────
  # Defined here so it's available for both geom and scale
  # ── Nodes: fill = LFC, border colour = compartment, size = FDR, shape = direction ──
  geom_node_point(
    aes(fill   = LFC,
        size   = neg_log10_fdr,
        shape  = direction,
        colour = compartment),   # border = subcellular compartment
    stroke = 1.2
  ) +
  scale_fill_gradient2(
    low      = "#2D87BB",   # blue  = neg LFC → depleted → gain of function
    mid      = "grey92",
    high     = "#E84545",   # red   = pos LFC → enriched → loss of function
    midpoint = 0,
    name     = "LFC (CRISPRi)"
  ) +
  scale_colour_manual(
    values = c(
      "Cell surface"  = "#E69F00",   # orange
      "Nucleus"       = "#9B5DE5",   # purple
      "Mitochondria"  = "#E84545",   # red
      "Endosome"      = "#56B4E9",   # sky blue
      "Lysosome"      = "#00B4D8",   # cyan
      "ER"            = "#2A9D8F",   # teal
      "Golgi"         = "#F4A261",   # peach
      "Cytoplasm"     = "#90BE6D",   # green
      "Cytoskeleton"  = "#F9C74F",   # yellow
      "Extracellular" = "#C77DFF",   # lavender
      "Other"         = "grey80"
    ),
    name   = "Compartment",
    guide  = guide_legend(override.aes = list(fill = "grey70", size = 3))
  ) +
  scale_size_continuous(range = c(3, 10), name = "−log₁₀(FDR)") +
  scale_shape_manual(
    values = c(
      "Depleted (gain of function)" = 21,
      "Enriched (loss of function)" = 24
    ),
    name = "CRISPRi hit"
  ) +

  # ── Hub gene labels (larger, bold) ──
  geom_node_text(
    data = function(x) filter(x, is_hub & !is_isolated),
    aes(label = name),
    repel        = TRUE,
    size         = 3.2,
    fontface     = "bold",
    colour       = "grey5",
    max.overlaps = 30,
    box.padding  = 0.35
  ) +
  # ── Non-hub connected labels ──
  geom_node_text(
    data = function(x) filter(x, !is_hub & !is_isolated),
    aes(label = name),
    repel        = TRUE,
    size         = 2.3,
    colour       = "grey30",
    max.overlaps = 20,
    box.padding  = 0.25
  ) +
  # ── Isolated node labels (italic, dashed border) ──
  geom_node_text(
    data = function(x) filter(x, is_isolated),
    aes(label = name),
    repel        = TRUE,
    size         = 2.2,
    fontface     = "italic",
    colour       = "grey50",
    max.overlaps = 20,
    box.padding  = 0.25
  ) +

  # ── Hub gene ring highlight — extra outer ring in black ──────
  # Uses a slightly larger fixed size to create a visible outer ring
  # colour is already mapped, so we draw this as a separate point layer
  geom_node_point(
    data = function(x) filter(x, is_hub),
    aes(size = neg_log10_fdr),
    fill        = NA,
    colour      = "black",
    shape       = 21,
    stroke      = 2.2,       # thicker than node stroke (1.2) to show outside it
    show.legend = FALSE
  ) +

  labs(
    title    = paste("CRISPRi Gene Network —", CELL_TYPE),
    subtitle = sprintf(
      "%d genes  |  %d interactions (STRING ≥ %d)  |  %d Louvain communities  |  top %d hubs outlined",
      vcount(g), ecount(g), STRING_SCORE, n_communities, HUB_TOP_N
    ),
    caption  = paste0(
      "● Depleted (neg LFC) = knockdown hurts T cells → gain-of-function gene   ",
      "▲ Enriched (pos LFC) = knockdown helps T cells → checkpoint/suppressor gene\n",
      "Node fill = LFC   Node border colour = subcellular compartment (GO CC)   ",
      "Coloured hulls = ", PATHWAY_SOURCE, " pathway   ",
      "Bold labels + thick black ring = top hub/regulator genes"
    )
  ) +

  theme_graph(base_family = "sans") +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.45, "cm"),
    legend.text      = element_text(size = 7.5),
    legend.title     = element_text(size = 8, face = "bold"),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 8.5, colour = "grey40"),
    plot.caption     = element_text(size = 7,   colour = "grey50"),
    plot.margin      = margin(10, 10, 10, 10)
  )

# ── 9. Hub / regulator table plot ────────────────────────────
hub_tbl <- tibble(
    Gene        = V(g)$name,
    direction   = V(g)$direction,
    LFC         = V(g)$LFC,
    FDR         = V(g)$FDR,
    degree      = V(g)$degree,
    betweenness = V(g)$betweenness,
    eigenvector = V(g)$eigenvector,
    hub_score   = V(g)$hub_score,
    pathway     = V(g)$pathway
  ) %>%
  arrange(desc(hub_score)) %>%
  slice_head(n = HUB_TOP_N)

p_hub <- hub_tbl %>%
  mutate(
    Gene      = fct_reorder(Gene, hub_score),
    dir_col   = ifelse(LFC < 0, "#2D87BB", "#E84545")
  ) %>%
  ggplot(aes(x = hub_score, y = Gene)) +
  geom_segment(
    aes(x = 0, xend = hub_score, yend = Gene, colour = direction),
    linewidth = 1.2
  ) +
  geom_point(aes(colour = direction, size = -log10(FDR + 1e-10))) +
  geom_text(
    aes(label = pathway),
    hjust  = -0.1,
    size   = 2.4,
    colour = "grey40"
  ) +
  scale_colour_manual(
    values = c(
      "Depleted (gain of function)" = "#2D87BB",
      "Enriched (loss of function)" = "#E84545"
    ),
    guide = "none"
  ) +
  scale_size_continuous(range = c(2, 6), name = "−log₁₀(FDR)") +
  labs(
    title    = paste("Top", HUB_TOP_N, "Hub / Master Regulator Genes"),
    subtitle = "Composite score = mean z-score of degree, betweenness, eigenvector, PageRank",
    x        = "Composite hub score (z)",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey40"),
    axis.text.y   = element_text(size = 9)
  )

# ── 10. Pathway ORA dot plot ──────────────────────────────────
if (nrow(pathway_df) > 0) {
  p_pathway <- pathway_df %>%
    mutate(
      pathway_label = fct_reorder(pathway_label, -p.adjust),
      gene_ratio_num = map_dbl(GeneRatio, ~ eval(parse(text = .x)))
    ) %>%
    ggplot(aes(x = gene_ratio_num, y = pathway_label)) +
    geom_point(aes(size = Count, colour = p.adjust)) +
    scale_colour_gradient(
      low  = "#E84545",
      high = "grey70",
      name = "FDR",
      trans = "log10"
    ) +
    scale_size_continuous(range = c(3, 9), name = "Gene count") +
    labs(
      title    = paste(PATHWAY_SOURCE, "pathway enrichment (ORA)"),
      subtitle = paste0("Input: ", length(gene_list), " network genes"),
      x        = "Gene ratio",
      y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title  = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 9)
    )
} else {
  p_pathway <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No enriched pathways at FDR < 0.05\nTry relaxing PATHWAY_FDR or expanding TOP_N_GENES",
             size = 4, colour = "grey50") +
    theme_void()
}

# ── 11. Assemble & save ───────────────────────────────────────
print(p_network)

side_panel <- p_hub / p_pathway + plot_layout(heights = c(1, 1))
combined   <- p_network | side_panel
combined   <- combined + plot_layout(widths = c(3, 1.4))

ggsave("CRISPR_Carnevale_network.pdf",   p_network, width = 15, height = 12)
ggsave("CRISPR_Carnevale_network.png",   p_network, width = 15, height = 12, dpi = 300)
ggsave("CRISPR_Carnevale_combined.pdf",  combined,  width = 22, height = 12)
ggsave("CRISPR_Carnevale_combined.png",  combined,  width = 22, height = 12, dpi = 300)
cat("Saved: CRISPR_Carnevale_network.pdf/.png and CRISPR_Carnevale_combined.pdf/.png\n")

# ── 12. Optional: export to Cytoscape ────────────────────────
if (EXPORT_CYTOSCAPE) {
  if (!requireNamespace("RCy3", quietly = TRUE)) stop("Install RCy3: BiocManager::install('RCy3')")
  library(RCy3)
  cytoscapePing()
  createNetworkFromIgraph(g, title = paste("CRISPRi", CELL_TYPE), collection = "CRISPR Networks")

  style_name <- "CRISPRi_style"
  createVisualStyle(style_name)
  setNodeColorMapping("LFC",        c(min(V(g)$LFC), 0, max(V(g)$LFC)),
                      c("#2D87BB","#FFFFFF","#E84545"), "c", style_name)
  setNodeSizeMapping("neg_log10_fdr", range(V(g)$neg_log10_fdr), c(20, 80), "c", style_name)
  setVisualStyle(style_name)
  cat("Network exported to Cytoscape\n")
}

# ── 13. Save tables ───────────────────────────────────────────
node_tbl <- tibble(
  Gene        = V(g)$name,
  direction   = V(g)$direction,
  LFC         = V(g)$LFC,
  FDR         = V(g)$FDR,
  good_sgrna  = V(g)$good_sgrna,
  pathway     = V(g)$pathway,
  community   = V(g)$community,
  degree      = V(g)$degree,
  betweenness = round(V(g)$betweenness, 4),
  eigenvector = round(V(g)$eigenvector, 4),
  pagerank    = round(V(g)$pagerank,    6),
  hub_score   = round(V(g)$hub_score,   4),
  is_hub        = V(g)$is_hub,
  is_isolated   = V(g)$is_isolated,
  in_ppi        = !V(g)$is_isolated,
  compartment   = V(g)$compartment,      # subcellular location (GO CC)
  complexes     = V(g)$complexes,        # CORUM protein complex membership
  n_complexes   = V(g)$n_complexes
) %>% arrange(desc(hub_score))

# Unmapped genes — save separately so they are not silently lost
if (length(genes_unmapped) > 0) {
  unmapped_tbl <- hits %>%
    filter(Gene %in% genes_unmapped) %>%
    mutate(reason = "Not found in STRING database")
  write_csv(unmapped_tbl, "CRISPR_Carnevale_unmapped_genes.csv")
  cat(sprintf("Saved: CRISPR_Carnevale_unmapped_genes.csv (%d genes)\n",
              nrow(unmapped_tbl)))
}

edge_tbl <- as_data_frame(g, what = "edges")

write_csv(node_tbl, "CRISPR_Carnevale_network_nodes.csv")
write_csv(edge_tbl, "CRISPR_Carnevale_network_edges.csv")

if (nrow(pathway_df) > 0)
  write_csv(pathway_df, "CRISPR_Carnevale_pathways.csv")

cat("Saved: CRISPR_Carnevale_network_nodes.csv / _edges.csv / _pathways.csv\n")

# ── 14. Console: top regulators ──────────────────────────────
cat("\n══ Top", HUB_TOP_N, "hub / regulator genes ══════════════════════════\n")
node_tbl %>%
  dplyr::select(Gene, direction, LFC, FDR, pathway, degree, hub_score, is_hub) %>%
  slice_head(n = HUB_TOP_N) %>%
  print(n = HUB_TOP_N)
