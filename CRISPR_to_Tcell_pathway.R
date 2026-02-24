# ============================================================
# Hit-Centric CRISPRi T Cell Pathway Diagram
# Carnevale et al. — MAGeCK gene-level results
#
# Design:
#   1. Load ALL significant hits (neg.fdr<0.10 & neg.lfc<-0.20,
#      or pos.fdr<0.10 & pos.lfc>0.20)
#   2. Map every hit to a curated x/y position on a T cell
#      pathway canvas using a lookup table.
#   3. Non-hit "connector" nodes are placed ONLY to bridge
#      two hits; they are small, grey, and barely labelled.
#   4. Remaining hits that have no curated position are shown
#      in a colour-coded overflow strip beneath the canvas.
#   5. Reactome ORA for subtitle annotation only.
#
# CRISPRi convention (Carnevale):
#   depleted (neg.lfc < 0) = GoF / essential  → BLUE
#   enriched (pos.lfc > 0) = LoF / suppressor → RED
# ============================================================

# ── 0.  Packages ─────────────────────────────────────────────
needed_bioc <- c("clusterProfiler", "ReactomePA", "org.Hs.eg.db")
needed_cran <- c("ggplot2", "ggforce", "ggrepel", "ggnewscale",
                 "dplyr", "tibble", "tidyr", "scales", "grid",
                 "forcats", "stringr", "purrr")
new_b <- needed_bioc[!sapply(needed_bioc, requireNamespace, quietly = TRUE)]
if (length(new_b)) BiocManager::install(new_b)
new_c <- needed_cran[!sapply(needed_cran, requireNamespace, quietly = TRUE)]
if (length(new_c)) install.packages(new_c)

# Load Bioc first (they mask dplyr); tidyverse last
library(clusterProfiler); library(ReactomePA); library(org.Hs.eg.db)
library(ggplot2); library(ggforce); library(ggrepel); library(ggnewscale)
library(dplyr); library(tibble); library(tidyr); library(scales)
library(grid); library(forcats); library(stringr); library(purrr)

# Pin every dplyr/tidyr verb against Bioc masking
select     <- dplyr::select;   filter     <- dplyr::filter
rename     <- dplyr::rename;   mutate     <- dplyr::mutate
arrange    <- dplyr::arrange;  summarise  <- dplyr::summarise
group_by   <- dplyr::group_by; ungroup    <- dplyr::ungroup
left_join  <- dplyr::left_join; bind_rows <- dplyr::bind_rows
distinct   <- dplyr::distinct; pull       <- dplyr::pull
slice_head <- dplyr::slice_head; transmute <- dplyr::transmute
case_when  <- dplyr::case_when; anti_join  <- dplyr::anti_join
replace_na <- tidyr::replace_na
fct_reorder <- forcats::fct_reorder
str_trunc   <- stringr::str_trunc
map_chr     <- purrr::map_chr

# ── 1.  Settings ─────────────────────────────────────────────
FDR_CUT <- 0.10
LFC_CUT <- 0.20
OUT     <- "C:/Users/difen/Rcode/Tcell_CRISPRi_pathway"
TSTAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")   # e.g. 20250222_143512

# ── 1b.  Multi-screen configuration ──────────────────────────
# Add one row per screen. First row = primary / reference assay.
# To run single-screen mode: keep only one row, or set
#   MERGE_MODE <- FALSE.
#
# Conventions
#   "CRISPRi"       knockdown; neg.lfc<0 (depleted) → GoF/blue
#   "CRISPRko"      knockout;  same sign convention as CRISPRi
#   "CRISPRko_flip" knockout;  lab reverses neg/pos labels
#   "CRISPRa"       activation; neg/pos flipped vs CRISPRi
#
# Merge options
#   MERGE_METHOD  "union"     gene hit in ANY screen → on diagram
#                 "intersect" gene hit in ALL screens → on diagram
#   MERGE_RULE    "majority"  direction with most votes wins
#                 "strict"    only fully concordant genes get direction
#                 "first"     direction from row-1 assay always wins
#
MERGE_MODE   <- TRUE         # FALSE = single-screen; TRUE = merge all rows
MERGE_METHOD <- "union"
MERGE_RULE   <- "majority"

ASSAY_CONFIG <- tribble(
  ~path,                                                                                ~name,               ~convention,  ~toupper_genes,
  "C:/Users/difen/Downloads/Carnevale_CRISPRi_Tcell_mageck_gene_summary.txt",          "Carnevale_CRISPRi", "CRISPRi",    FALSE,
  "C:/Users/difen/Downloads/Shifrut_CRISPRi_CD8T_mageck_gene_summary.txt",             "Shifrut_CRISPRi",   "CRISPRi",    FALSE  # human CD8 T cell
  # ── Add further screens by uncommenting and filling these rows: ─────────────────────
  # "C:/Users/difen/Downloads/Belk_CRISPRi_CD8T_mageck_gene_summary.txt",              "Belk_CRISPRi",      "CRISPRi",    TRUE   # mouse – needs toupper
  # "C:/Users/difen/Downloads/Screen3_gene_summary.txt",                               "Screen3",           "CRISPRi",    FALSE
)

# ── 2.  Load & classify hits ─────────────────────────────────
source("C:/Users/difen/Rcode/CRISPR_merge_utils.R")

# Build assay list from config table
assay_list <- purrr::pmap(
  ASSAY_CONFIG,
  function(path, name, convention, ...) {
    load_mageck(path       = path,
                name       = name,
                convention = convention,
                fdr_cut    = FDR_CUT,
                lfc_cut    = LFC_CUT)
  }
)

if (MERGE_MODE && length(assay_list) > 1) {
  # ── Multi-screen path ─────────────────────────────────────
  mg       <- merge_assays(assay_list,
                           method         = MERGE_METHOD,
                           consensus_rule = MERGE_RULE)
  assay_summary(mg)          # concordance report to console
  dat      <- mg$dat
  hits     <- mg$hits
  hits_dep <- mg$hits_dep
  hits_enr <- mg$hits_enr
  MERGE_LABEL <- sprintf("%s  [%s, n=%d screens]",
                         paste(mg$assay_names, collapse = " + "),
                         MERGE_METHOD, length(assay_list))
  cat(sprintf("MERGE_MODE active: %s\n", MERGE_LABEL))
} else {
  # ── Single-screen path ───────────────────────────────────
  a1       <- assay_list[[1]]
  dat      <- a1$dat
  hits     <- a1$hits
  hits_dep <- a1$hits_dep
  hits_enr <- a1$hits_enr
  mg       <- NULL
  MERGE_LABEL <- a1$name
  if (length(assay_list) > 1)
    message("Note: MERGE_MODE is FALSE — only '", a1$name,
            "' used. Set MERGE_MODE <- TRUE to activate merge.")
}
cat(sprintf("Hits: %d depleted (GoF), %d enriched (LoF)\n",
            nrow(hits_dep), nrow(hits_enr)))

# ── 3.  ORA (Reactome) — for subtitle only ───────────────────
emap <- tryCatch(
  bitr(hits$gene, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db),
  error = function(e) data.frame(SYMBOL = character(), ENTREZID = character()))

ora <- tryCatch(
  enrichPathway(gene = emap$ENTREZID, organism = "human",
                pvalueCutoff = 0.10, qvalueCutoff = 0.25, readable = TRUE),
  error = function(e) NULL)

if (!is.null(ora) && nrow(as.data.frame(ora)) > 0) {
  ora_df  <- as.data.frame(ora) %>% arrange(p.adjust) %>% slice_head(n = 20)
  ora_sub <- paste(str_trunc(head(ora_df$Description, 3), 50), collapse = " · ")
  cat(sprintf("ORA: %d pathways enriched\n", nrow(ora_df)))
} else {
  ora_df  <- data.frame()
  ora_sub <- "ORA: no significant pathways"
}

# ── 4.  Hit position lookup table ────────────────────────────
# Every gene here gets a node ON the pathway canvas.
# Columns: gene, x, y, shape (kinase/tf/adaptor/receptor/complex),
#          module (for colour band), label (display name)
#
# Canvas coordinate system:
#   x  0 –  3.4   TCR / proximal
#   x  3.6 – 7.0  MAPK / NFAT
#   x  7.2 – 10.5 PI3K / AKT / mTOR
#   x 10.7 – 14.0 NF-κB / transcription
#   x 14.2 – 17.5 JAK-STAT
#   x 17.7 – 21.5 Epigenetic / apoptosis
#
#   y  9.4 – 9.9   plasma membrane
#   y  7.0 – 9.3   proximal cytoplasm
#   y  4.5 – 6.9   mid cytoplasm
#   y  2.0 – 4.4   distal cytoplasm
#   y  0.3 – 1.9   nucleus / TFs

HIT_POS <- tribble(
  ~gene,      ~x,    ~y,   ~shape,     ~module,          ~label,

  # ── TCR / PROXIMAL ──────────────────────────────────────────────────────────
  # Membrane receptors
  "CD247",    0.6,   9.7, "receptor", "TCR/Proximal",   "CD3ζ",
  "CD3D",     1.2,   9.7, "receptor", "TCR/Proximal",   "CD3δ",
  "CD3E",     1.8,   9.7, "receptor", "TCR/Proximal",   "CD3ε",
  "CD28",     2.5,   9.5, "receptor", "TCR/Proximal",   "CD28",
  "PDCD1",    0.3,   9.4, "receptor", "TCR/Proximal",   "PD-1",
  "LAG3",     0.8,   9.1, "receptor", "TCR/Proximal",   "LAG3",
  "CD5",      2.3,   9.1, "receptor", "TCR/Proximal",   "CD5",
  "CTLA4",    3.0,   9.4, "receptor", "TCR/Proximal",   "CTLA4",
  "CD8A",     3.2,   9.7, "receptor", "TCR/Proximal",   "CD8A",
  # Proximal kinases
  "LCK",      1.7,   8.2, "kinase",   "TCR/Proximal",   "LCK",
  "FYN",      0.7,   8.2, "kinase",   "TCR/Proximal",   "FYN",
  "ZAP70",    1.2,   7.3, "kinase",   "TCR/Proximal",   "ZAP70",
  "ITK",      0.4,   6.5, "kinase",   "TCR/Proximal",   "ITK",
  # Adaptors / scaffold
  "LAT",      1.6,   6.8, "adaptor",  "TCR/Proximal",   "LAT",
  "LCP2",     0.7,   6.0, "adaptor",  "TCR/Proximal",   "SLP-76",
  "GRAP2",    2.3,   7.0, "adaptor",  "TCR/Proximal",   "GRAP2",
  "VAV1",     0.4,   5.2, "adaptor",  "TCR/Proximal",   "VAV1",
  "SH2D1A",   1.1,   5.2, "adaptor",  "TCR/Proximal",   "SAP",
  "RAC2",     0.4,   4.5, "kinase",   "TCR/Proximal",   "RAC2",
  "WAS",      1.1,   4.5, "adaptor",  "TCR/Proximal",   "WASp",
  "RHOH",     1.8,   4.5, "adaptor",  "TCR/Proximal",   "RhoH",
  # LoF suppressors — proximal
  "CBLB",     2.2,   7.5, "adaptor",  "TCR/Proximal",   "CBLB",
  "UBASH3A",  0.2,   7.0, "adaptor",  "TCR/Proximal",   "UBASH3A",
  "PTPN6",    2.8,   7.8, "adaptor",  "TCR/Proximal",   "SHP-1",

  # ── MAPK / RAS / NFAT ───────────────────────────────────────────────────────
  "PLCG1",    3.8,   6.8, "kinase",   "MAPK/NFAT",      "PLCγ1",
  "RASGRP1",  3.8,   5.9, "adaptor",  "MAPK/NFAT",      "RasGRP1",
  "KRAS",     3.8,   5.0, "kinase",   "MAPK/NFAT",      "RAS",
  "RAF1",     3.8,   4.2, "kinase",   "MAPK/NFAT",      "RAF",
  "MAP2K1",   3.8,   3.4, "kinase",   "MAPK/NFAT",      "MEK1/2",
  "MAPK3",    3.8,   2.5, "kinase",   "MAPK/NFAT",      "ERK1/2",
  "MAPK1",    3.8,   0.9, "tf",       "MAPK/NFAT",      "ERK(n)",
  "PPP3CA",   5.8,   5.5, "complex",  "MAPK/NFAT",      "Calcineurin",
  "NFATC1",   5.5,   0.9, "tf",       "MAPK/NFAT",      "NFAT",
  # LoF suppressor — MAPK
  "RASA2",    5.2,   5.0, "adaptor",  "MAPK/NFAT",      "RASA2",
  "DUSP4",    5.2,   2.5, "adaptor",  "MAPK/NFAT",      "DUSP4",

  # ── PI3K / AKT / mTOR ───────────────────────────────────────────────────────
  "PIK3CD",   7.5,   6.5, "kinase",   "PI3K/AKT",       "PI3Kδ",
  "PDPK1",    7.5,   5.7, "kinase",   "PI3K/AKT",       "PDK1",
  "AKT1",     7.5,   4.9, "kinase",   "PI3K/AKT",       "AKT",
  "MTOR",     7.5,   4.0, "complex",  "PI3K/AKT",       "mTOR",
  "FOXO1",    7.5,   0.9, "tf",       "PI3K/AKT",       "FOXO1",
  "CCND2",    8.5,   3.2, "adaptor",  "PI3K/AKT",       "CycD2",
  "CCND3",    9.2,   3.2, "adaptor",  "PI3K/AKT",       "CycD3",
  "CDK6",     9.8,   3.2, "kinase",   "PI3K/AKT",       "CDK6",
  # LoF
  "DGKZ",     8.8,   6.5, "kinase",   "PI3K/AKT",       "DGKζ",
  "DGKA",     9.5,   6.5, "kinase",   "PI3K/AKT",       "DGKα",
  "CDKN1B",   9.8,   0.9, "tf",       "PI3K/AKT",       "p27",

  # ── NF-κB / TRANSCRIPTION ───────────────────────────────────────────────────
  "PRKCQ",   10.9,   6.5, "kinase",   "NF-κB",          "PKCθ",
  "CARD11",  11.5,   5.7, "complex",  "NF-κB",          "CARMA1",
  "BCL10",   12.1,   4.9, "adaptor",  "NF-κB",          "BCL10",
  "MALT1",   12.8,   4.9, "adaptor",  "NF-κB",          "MALT1",
  "IKBKB",   11.8,   4.0, "kinase",   "NF-κB",          "IKKβ",
  "RELA",    11.2,   0.9, "tf",       "NF-κB",          "p65",
  "NFKB1",   12.0,   0.9, "tf",       "NF-κB",          "p50",
  "IRF4",    10.8,   0.9, "tf",       "NF-κB",          "IRF4",
  "GATA3",   12.9,   0.9, "tf",       "NF-κB",          "GATA3",
  "JUNB",    13.7,   0.9, "tf",       "NF-κB",          "JUNB",
  "YBX1",    13.7,   1.7, "tf",       "NF-κB",          "YBX1",
  "HIVEP2",  10.8,   1.7, "tf",       "NF-κB",          "HIVEP2",
  "IKZF2",   11.6,   1.7, "tf",       "NF-κB",          "IKZF2",
  "NR4A1",   12.4,   1.7, "tf",       "NF-κB",          "NR4A1",
  "TOB1",    13.1,   1.7, "tf",       "NF-κB",          "TOB1",
  # LoF
  "TNFAIP3", 13.5,   4.9, "adaptor",  "NF-κB",          "A20",
  "TCEB2",   13.5,   4.0, "adaptor",  "NF-κB",          "Elongin-B",
  "ARIH2",   13.5,   3.2, "adaptor",  "NF-κB",          "ARIH2",
  "SUPT4H1", 10.8,   4.9, "adaptor",  "NF-κB",          "SUPT4H1",

  # ── JAK-STAT / CYTOKINE ─────────────────────────────────────────────────────
  "JAK3",    14.5,   8.3, "kinase",   "JAK-STAT",       "JAK3",
  "JAK1",    15.2,   8.3, "kinase",   "JAK-STAT",       "JAK1",
  "JAK2",    17.0,   8.3, "kinase",   "JAK-STAT",       "JAK2",
  "STAT5A",  14.5,   7.3, "tf",       "JAK-STAT",       "STAT5A",
  "STAT5B",  15.2,   7.3, "tf",       "JAK-STAT",       "STAT5B",
  "STAT1",   17.0,   7.3, "tf",       "JAK-STAT",       "STAT1",
  # LoF
  "SOCS1",   14.0,   6.3, "adaptor",  "JAK-STAT",       "SOCS1",
  "IRF2BP2", 16.2,   0.9, "tf",       "JAK-STAT",       "IRF2BP2",

  # ── EPIGENETIC / CHROMATIN / APOPTOSIS ──────────────────────────────────────
  "EZH2",    18.2,   2.5, "complex",  "Epi/Apoptosis",  "EZH2",
  "DNMT3A",  19.5,   2.5, "adaptor",  "Epi/Apoptosis",  "DNMT3A",
  "ARID1A",  20.2,   2.5, "adaptor",  "Epi/Apoptosis",  "ARID1A",
  "SMARCB1", 20.9,   2.5, "adaptor",  "Epi/Apoptosis",  "SMARCB1",
  "MYC",     18.5,   0.9, "tf",       "Epi/Apoptosis",  "MYC",
  "MEF2D",   17.8,   0.9, "tf",       "Epi/Apoptosis",  "MEF2D",
  "TP53",    21.5,   0.9, "tf",       "Epi/Apoptosis",  "p53",
  "EOMES",   19.3,   0.9, "tf",       "Epi/Apoptosis",  "EOMES",
  "TBX21",   20.1,   0.9, "tf",       "Epi/Apoptosis",  "T-bet",
  "TCF7",    20.9,   0.9, "tf",       "Epi/Apoptosis",  "TCF1",
  "PRDM1",   19.3,   1.7, "tf",       "Epi/Apoptosis",  "BLIMP1",
  "BATF",    20.1,   1.7, "tf",       "Epi/Apoptosis",  "BATF",
  "BACH2",   20.9,   1.7, "tf",       "Epi/Apoptosis",  "BACH2",
  "BCL2",    18.2,   4.8, "adaptor",  "Epi/Apoptosis",  "BCL2",
  "BCL2L1",  19.0,   5.5, "adaptor",  "Epi/Apoptosis",  "BCL-XL",
  "MCL1",    17.7,   5.5, "adaptor",  "Epi/Apoptosis",  "MCL1",
  "BAX",     19.8,   5.5, "adaptor",  "Epi/Apoptosis",  "BAX",
  # LoF
  "TMEM222", 21.5,   2.5, "adaptor",  "Epi/Apoptosis",  "TMEM222",
  "MEF2D",   17.8,   0.9, "tf",       "Epi/Apoptosis",  "MEF2D",

  # ── MITOCHONDRIAL / METABOLIC ────────────────────────────────────────────────
  # GoF (depleted – blue)
  "PHB",     22.6,   6.2, "complex",  "Mito/Metabolic", "PHB",
  "PHB2",    23.4,   6.2, "complex",  "Mito/Metabolic", "PHB2",
  "SDHB",    24.2,   6.2, "complex",  "Mito/Metabolic", "SDHB",
  "NFS1",    25.0,   6.2, "adaptor",  "Mito/Metabolic", "NFS1",
  "RPTOR",   22.6,   5.0, "complex",  "Mito/Metabolic", "RPTOR",
  "VHL",     23.4,   5.0, "adaptor",  "Mito/Metabolic", "VHL",
  "GCLC",    24.2,   5.0, "adaptor",  "Mito/Metabolic", "GCLC",
  "GCLM",    25.0,   5.0, "adaptor",  "Mito/Metabolic", "GCLM",
  "AK2",     24.6,   7.2, "kinase",   "Mito/Metabolic", "AK2",
  "PRMT5",   22.6,   2.5, "adaptor",  "Mito/Metabolic", "PRMT5",
  "HUWE1",   23.4,   2.5, "adaptor",  "Mito/Metabolic", "HUWE1",
  "KMT2D",   24.2,   2.5, "adaptor",  "Mito/Metabolic", "KMT2D",
  # LoF (enriched – red)
  "GNA13",   22.6,   0.9, "adaptor",  "Mito/Metabolic", "GNA13",
  "NDUFB10", 23.4,   6.9, "adaptor",  "Mito/Metabolic", "NDUFB10",
  "GLRX",    24.2,   0.9, "adaptor",  "Mito/Metabolic", "GLRX",
  "FIBP",    25.0,   0.9, "adaptor",  "Mito/Metabolic", "FIBP",
  "AGO1",    22.6,   1.7, "adaptor",  "Mito/Metabolic", "AGO1",
  "BIRC6",   23.4,   1.7, "adaptor",  "Mito/Metabolic", "BIRC6",
  "PCBP2",   24.2,   1.7, "adaptor",  "Mito/Metabolic", "PCBP2",
  "RPRD1B",  25.0,   1.7, "adaptor",  "Mito/Metabolic", "RPRD1B",
  "C2orf68", 25.0,   3.6, "adaptor",  "Mito/Metabolic", "C2orf68"
)

# Remove any duplicate gene rows (keep first)
HIT_POS <- HIT_POS %>% distinct(gene, .keep_all = TRUE)

# ── 5.  Connector nodes (non-hit grey bridges) ───────────────
# Appear ONLY to connect two hit nodes where biology requires it.
CONNECTORS <- tribble(
  ~gene,     ~x,    ~y,   ~shape,     ~module,         ~label,
  "pMHCII",  1.5,  10.5, "ligand",   "TCR/Proximal",  "pMHC",
  "TCR",     1.5,   9.9, "complex",  "TCR/Proximal",  "TCRαβ",
  "CD4",     2.1,   9.9, "receptor", "TCR/Proximal",  "CD4",
  "GRB2",    2.4,   6.3, "adaptor",  "TCR/Proximal",  "GRB2",
  "PIK3R1",  8.2,   7.2, "adaptor",  "PI3K/AKT",      "p85α",
  "PTEN",    9.0,   6.5, "adaptor",  "PI3K/AKT",      "PTEN",
  "CHUK",    12.8,  4.0, "kinase",   "NF-κB",         "IKKα",
  "NFKBIA",  12.3,  3.2, "adaptor",  "NF-κB",         "IκBα",
  "IL2RB",   14.5,  9.5, "receptor", "JAK-STAT",      "IL-2Rβ",
  "IL2RG",   15.2,  9.5, "receptor", "JAK-STAT",      "γc",
  "IFNGR1",  17.0,  9.5, "receptor", "JAK-STAT",      "IFNγR1",
  "STAT5n",  15.0,  0.4, "tf",       "JAK-STAT",      "STAT5(n)",
  "EED",     18.9,  3.2, "adaptor",  "Epi/Apoptosis", "EED",
  "SUZ12",   17.7,  3.2, "adaptor",  "Epi/Apoptosis", "SUZ12",
  "gene_exp",11.0, -0.5, "output",   "shared",        "Gene\nExpression"
)

# ── 5b. SH2-domain annotation ────────────────────────────────
# Proteins with at least one confirmed SH2 domain on this canvas.
# Source: UniProt domain annotations (InterPro SH2 family).
SH2_GENES <- c(
  # Proximal kinases
  "LCK",    # 1 SH2
  "FYN",    # 1 SH2
  "ZAP70",  # 2 SH2 (tandem)
  "ITK",    # 1 SH2
  # Adaptors / scaffold
  "VAV1",   # 1 SH2
  "SH2D1A", # SH2-only (SAP)
  "GRAP2",  # 1 SH2 (GADS)
  # Effectors
  "PLCG1",  # split SH2
  "PTPN6",  # 2 SH2 (SHP-1)
  # JAK-STAT
  "STAT5A", # 1 SH2
  "STAT5B", # 1 SH2
  "STAT1",  # 1 SH2
  "SOCS1",  # 1 SH2 (negative regulator)
  # Connector nodes
  "GRB2",   # 1 SH2
  "PIK3R1"  # 2 SH2 (p85α regulatory subunit)
)

# ── 6.  Build unified node table ─────────────────────────────
# Hit nodes first; connectors only if not already a hit
all_nodes <- bind_rows(
  HIT_POS %>% mutate(is_hit = TRUE),
  CONNECTORS %>%
    filter(!gene %in% HIT_POS$gene) %>%
    mutate(is_hit = FALSE)
)

# Join CRISPRi data
nodes <- all_nodes %>%
  left_join(dat %>% select(gene, hit_dir, best_lfc, best_fdr, nlp),
            by = "gene") %>%
  mutate(
    hit_dir    = replace_na(hit_dir, "none"),
    nlp        = replace_na(nlp, 0),
    # Default single-assay colours (overridden below for merged mode)
    fill_col   = case_when(
      hit_dir == "depleted" ~ "#2166AC",   # blue  – GoF
      hit_dir == "enriched" ~ "#D6604D",   # red   – LoF
      TRUE                  ~ "#D5D8DC"),  # grey  – connector
    border_col = case_when(
      hit_dir == "depleted" ~ "#08306B",
      hit_dir == "enriched" ~ "#67000D",
      TRUE                  ~ "#909497"),
    # Node radius: hits scale with −log10(FDR); connectors fixed tiny
    node_r     = case_when(
      hit_dir != "none" ~ rescale(pmin(nlp, 6), to = c(0.27, 0.50), from = c(0, 6)),
      TRUE              ~ 0.17),
    label_size = if_else(hit_dir != "none", 2.1, 1.55),
    label_face = if_else(hit_dir != "none", "bold", "plain"),
    label_col  = if_else(hit_dir != "none", "#111111", "#777777")
  )

# In MERGE_MODE override fill/border to encode concordance
if (!is.null(mg)) {
  merged_cols <- node_colours_merged(mg$hits, mg$assay_names)
  nodes <- nodes %>%
    left_join(merged_cols, by = "gene", suffix = c("", ".mg")) %>%
    mutate(
      fill_col   = if_else(!is.na(fill_col.mg),   fill_col.mg,   fill_col),
      border_col = if_else(!is.na(border_col.mg), border_col.mg, border_col)
    ) %>%
    select(-fill_col.mg, -border_col.mg)
}

# Tag nodes that carry an SH2 domain
nodes <- nodes %>% mutate(has_sh2 = gene %in% SH2_GENES)

# Report coverage
placed <- nodes %>% filter(hit_dir != "none")
cat(sprintf("\nHits placed on pathway : %d / %d (%.0f%%)\n",
            nrow(placed), nrow(hits),
            100 * nrow(placed) / nrow(hits)))

# Overflow: hits NOT in HIT_POS
overflow_hits <- hits %>% filter(!gene %in% HIT_POS$gene)
cat(sprintf("Overflow strip         : %d hits\n", nrow(overflow_hits)))

# ── 7.  Edges ────────────────────────────────────────────────
# All genes referenced must appear in nodes (hit or connector).
# is_cross = TRUE → drawn as a curved arc (cross-module link)
EDGES <- tribble(
  ~from,      ~to,        ~etype,             ~is_cross,
  # pMHC → TCR complex
  "pMHCII",  "TCR",      "binding",           FALSE,
  "TCR",     "CD247",    "binding",           FALSE,
  "TCR",     "CD3D",     "binding",           FALSE,
  "TCR",     "CD3E",     "binding",           FALSE,
  "CD4",     "LCK",      "activation",        FALSE,
  "CD28",    "LCK",      "activation",        FALSE,
  # Checkpoint inhibition
  "PDCD1",   "ZAP70",    "inhibition",        FALSE,
  "LAG3",    "LCK",      "inhibition",        FALSE,
  "CTLA4",   "LCK",      "inhibition",        FALSE,
  "CD5",     "LCK",      "inhibition",        FALSE,
  "UBASH3A", "ZAP70",    "inhibition",        FALSE,
  "CBLB",    "ZAP70",    "inhibition",        FALSE,
  "CBLB",    "LAT",      "inhibition",        FALSE,
  "PTPN6",   "LCK",      "inhibition",        FALSE,
  # Proximal signalling
  "LCK",     "ZAP70",    "phosphorylation",   FALSE,
  "FYN",     "ZAP70",    "phosphorylation",   FALSE,
  "ZAP70",   "LAT",      "phosphorylation",   FALSE,
  "ZAP70",   "LCP2",     "phosphorylation",   FALSE,
  "LAT",     "PLCG1",    "activation",        FALSE,
  "LAT",     "GRB2",     "binding",           FALSE,
  "LAT",     "ITK",      "activation",        FALSE,
  "LAT",     "GRAP2",    "binding",           FALSE,
  "LCP2",    "VAV1",     "activation",        FALSE,
  "LCP2",    "ITK",      "binding",           FALSE,
  "LCP2",    "RASGRP1",  "activation",        FALSE,
  "ITK",     "PLCG1",    "phosphorylation",   FALSE,
  "SH2D1A",  "LCK",      "activation",        FALSE,
  "VAV1",    "RAC2",     "activation",        FALSE,
  "RAC2",    "WAS",      "activation",        FALSE,
  # MAPK cascade
  "PLCG1",   "RASGRP1",  "activation",        FALSE,
  "PLCG1",   "PPP3CA",   "indirect_act",      FALSE,
  "RASGRP1", "KRAS",     "activation",        FALSE,
  "RASA2",   "KRAS",     "inhibition",        FALSE,
  "KRAS",    "RAF1",     "activation",        FALSE,
  "RAF1",    "MAP2K1",   "phosphorylation",   FALSE,
  "MAP2K1",  "MAPK3",    "phosphorylation",   FALSE,
  "MAPK3",   "MAPK1",    "translocation",     FALSE,
  "DUSP4",   "MAPK3",    "inhibition",        FALSE,
  "PPP3CA",  "NFATC1",   "activation",        FALSE,
  "NFATC1",  "gene_exp", "activation",        FALSE,
  "MAPK1",   "gene_exp", "activation",        FALSE,
  # PI3K / AKT
  "CD28",    "PIK3CD",   "activation",        TRUE,
  "VAV1",    "PIK3CD",   "activation",        TRUE,
  "CBLB",    "PIK3CD",   "inhibition",        TRUE,
  "DGKZ",    "PIK3CD",   "inhibition",        FALSE,
  "DGKA",    "PIK3CD",   "inhibition",        FALSE,
  "PIK3R1",  "PIK3CD",   "binding",           FALSE,
  "PTEN",    "PIK3CD",   "inhibition",        FALSE,
  "PIK3CD",  "PDPK1",    "indirect_act",      FALSE,
  "PDPK1",   "AKT1",     "phosphorylation",   FALSE,
  "AKT1",    "MTOR",     "activation",        FALSE,
  "AKT1",    "FOXO1",    "phosphorylation",   FALSE,
  "AKT1",    "BCL2",     "activation",        TRUE,
  "AKT1",    "CDKN1B",   "phosphorylation",   FALSE,
  "MTOR",    "CCND2",    "indirect_act",      FALSE,
  "MTOR",    "MYC",      "indirect_act",      TRUE,
  "CCND2",   "CDK6",     "binding",           FALSE,
  "CCND3",   "CDK6",     "binding",           FALSE,
  "FOXO1",   "gene_exp", "inhibition",        FALSE,
  "CDKN1B",  "gene_exp", "inhibition",        FALSE,
  # NF-κB
  "PLCG1",   "PRKCQ",    "activation",        TRUE,
  "PRKCQ",   "CARD11",   "phosphorylation",   FALSE,
  "CARD11",  "BCL10",    "binding",           FALSE,
  "BCL10",   "MALT1",    "binding",           FALSE,
  "MALT1",   "IKBKB",    "activation",        FALSE,
  "TNFAIP3", "MALT1",    "inhibition",        FALSE,
  "TNFAIP3", "IKBKB",    "inhibition",        FALSE,
  "IKBKB",   "NFKBIA",   "phosphorylation",   FALSE,
  "CHUK",    "NFKBIA",   "phosphorylation",   FALSE,
  "NFKBIA",  "RELA",     "inhibition",        FALSE,
  "NFKBIA",  "NFKB1",    "inhibition",        FALSE,
  "TCEB2",   "IKBKB",    "inhibition",        FALSE,
  "ARIH2",   "NFKBIA",   "indirect_act",      FALSE,
  "RELA",    "gene_exp", "activation",        FALSE,
  "NFKB1",   "gene_exp", "activation",        FALSE,
  "IRF4",    "gene_exp", "activation",        FALSE,
  "GATA3",   "gene_exp", "activation",        FALSE,
  "JUNB",    "gene_exp", "activation",        FALSE,
  "YBX1",    "gene_exp", "activation",        FALSE,
  "HIVEP2",  "gene_exp", "activation",        FALSE,
  "IKZF2",   "gene_exp", "inhibition",        FALSE,
  "NR4A1",   "gene_exp", "activation",        FALSE,
  "TOB1",    "gene_exp", "inhibition",        FALSE,
  "SUPT4H1", "gene_exp", "activation",        FALSE,
  # JAK-STAT
  "IL2RB",   "JAK1",     "activation",        FALSE,
  "IL2RG",   "JAK3",     "activation",        FALSE,
  "IFNGR1",  "JAK2",     "activation",        FALSE,
  "JAK3",    "STAT5A",   "phosphorylation",   FALSE,
  "JAK1",    "STAT5B",   "phosphorylation",   FALSE,
  "JAK2",    "STAT1",    "phosphorylation",   FALSE,
  "SOCS1",   "JAK1",     "inhibition",        FALSE,
  "SOCS1",   "JAK3",     "inhibition",        FALSE,
  "SOCS1",   "AKT1",     "inhibition",        TRUE,
  "STAT5A",  "STAT5n",   "translocation",     FALSE,
  "STAT5B",  "STAT5n",   "translocation",     FALSE,
  "STAT5n",  "gene_exp", "activation",        FALSE,
  "STAT1",   "gene_exp", "activation",        FALSE,
  "IRF2BP2", "gene_exp", "inhibition",        FALSE,
  # Epigenetic / apoptosis
  "EED",     "EZH2",     "binding",           FALSE,
  "SUZ12",   "EZH2",     "binding",           FALSE,
  "ARID1A",  "EZH2",     "inhibition",        FALSE,
  "SMARCB1", "EZH2",     "inhibition",        FALSE,
  "EZH2",    "gene_exp", "inhibition",        FALSE,
  "DNMT3A",  "gene_exp", "inhibition",        FALSE,
  "MYC",     "gene_exp", "activation",        FALSE,
  "MEF2D",   "gene_exp", "activation",        FALSE,
  "EOMES",   "gene_exp", "activation",        FALSE,
  "TBX21",   "gene_exp", "activation",        FALSE,
  "TCF7",    "gene_exp", "activation",        FALSE,
  "PRDM1",   "gene_exp", "inhibition",        FALSE,
  "BATF",    "gene_exp", "activation",        FALSE,
  "BACH2",   "gene_exp", "inhibition",        FALSE,
  "BCL2",    "BAX",      "inhibition",        FALSE,
  "BCL2L1",  "BAX",      "inhibition",        FALSE,
  "MCL1",    "BAX",      "inhibition",        FALSE,
  "TP53",    "BAX",      "activation",        FALSE,
  "TP53",    "BCL2",     "inhibition",        FALSE,
  # Mito/Metabolic
  "MTOR",    "RPTOR",    "binding",           TRUE,
  "RPTOR",   "PHB",      "indirect_act",      FALSE,
  "RPTOR",   "SDHB",     "indirect_act",      FALSE,
  "VHL",     "SDHB",     "binding",           FALSE,
  "NDUFB10", "AK2",      "binding",           FALSE,
  "PHB",     "PHB2",     "binding",           FALSE,
  "GCLC",    "GCLM",     "binding",           FALSE,
  "PRMT5",   "KMT2D",    "indirect_act",      FALSE,
  "HUWE1",   "VHL",      "indirect_act",      FALSE,
  "BIRC6",   "HUWE1",    "inhibition",        FALSE,
  "AGO1",    "PCBP2",    "binding",           FALSE,
  "GNA13",   "RPTOR",    "indirect_act",      TRUE,
  "AKT1",    "RPTOR",    "activation",        TRUE
)

# Keep only edges where both endpoints exist in nodes
node_genes <- unique(nodes$gene)
edges_filt <- EDGES %>%
  filter(from %in% node_genes, to %in% node_genes)

# Add coordinates
pos <- nodes %>% select(gene, x, y)
edges_coord <- edges_filt %>%
  left_join(pos, by = c("from" = "gene")) %>% rename(x1 = x, y1 = y) %>%
  left_join(pos, by = c("to"   = "gene")) %>% rename(x2 = x, y2 = y) %>%
  filter(!is.na(x1), !is.na(x2), !(x1 == x2 & y1 == y2))

# Edge style lookup
ESTYLE <- tribble(
  ~etype,             ~ecol,      ~elty,
  "activation",       "#1A5276",  "solid",
  "phosphorylation",  "#784212",  "solid",
  "inhibition",       "#7B241C",  "solid",
  "indirect_act",     "#1A5276",  "dashed",
  "binding",          "#5D6D7E",  "dashed",
  "translocation",    "#6C3483",  "dotted"
)
edges_coord <- edges_coord %>%
  left_join(ESTYLE, by = "etype") %>%
  mutate(ecol = replace_na(ecol, "#888888"),
         elty = replace_na(elty, "solid"))

# ── 8.  Overflow strip layout ─────────────────────────────────
strip_group <- function(g) {
  # Genes now placed in Mito/Metabolic module – exclude from strip
  mito_placed <- c("PHB","PHB2","SDHB","NFS1","VHL","RPTOR","PRMT5","HUWE1",
                   "KMT2D","GCLC","GCLM","AK2","GNA13","NDUFB10","GLRX",
                   "FIBP","AGO1","BIRC6","PCBP2","RPRD1B","C2orf68")
  case_when(
    g %in% mito_placed                             ~ NA_character_,
    grepl("^RP[SL]|^FAU|^MRPL|^MRPS",         g) ~ "Ribosome",
    grepl("^EIF|^EEF|TRMT|GSPT",               g) ~ "Translation",
    grepl("^POLR|^TAF|NELF|^GTF|^MED|SUPT",   g) ~ "Txn machinery",
    grepl("^NOP|^DDX|^DHX|^SNRP|^SF3|^LSM",   g) ~ "RNA processing",
    g %in% c("CCND1","CCND2","CDK4","RB1","E2F1",
              "CDKN1A","CDKN2A")                  ~ "Cell cycle",
    g %in% c("DAD1","DAPK1","XIAP","BIRC2",
              "BIRC3")                             ~ "Apoptosis",
    g %in% c("GNAI2","GNAQ","GNB1")               ~ "G-protein",
    g %in% c("AGO2","DICER1")                      ~ "RNA regulators",
    g %in% c("TMEM222")                            ~ "Other",
    TRUE ~ "Signalling"
  )
}

# filter out NA groups (mito-placed genes won't overflow)

if (nrow(overflow_hits) > 0) {
  overflow <- overflow_hits %>%
    mutate(grp = strip_group(gene)) %>%
    filter(!is.na(grp)) %>%                          # drop mito-placed genes
    arrange(grp, hit_dir, desc(abs(best_lfc)))

  # Pack groups left-to-right
  grp_tbl <- overflow %>%
    group_by(grp) %>% summarise(n = n(), .groups = "drop") %>%
    mutate(x_start = cumsum(c(0, head(n, -1))) * 0.90 + 0.4)

  overflow <- overflow %>%
    left_join(grp_tbl, by = "grp") %>%
    group_by(grp) %>%
    mutate(xi = (row_number() - 1) * 0.90 + x_start) %>%
    ungroup() %>%
    mutate(fill_col   = if_else(hit_dir == "depleted", "#2166AC", "#D6604D"),
           border_col = if_else(hit_dir == "depleted", "#08306B", "#67000D"))

  # Group label x positions
  grp_lab <- overflow %>%
    group_by(grp) %>% summarise(lx = mean(xi), .groups = "drop")
}

# ── 9.  Module colour reference ───────────────────────────────
MOD_COL <- c(
  "TCR/Proximal"  = "#2471A3",
  "MAPK/NFAT"     = "#1A7A4A",
  "PI3K/AKT"      = "#7D3C98",
  "NF-κB"         = "#BA4A00",
  "JAK-STAT"      = "#922B21",
  "Epi/Apoptosis" = "#0E7A6E",
  "Mito/Metabolic"= "#B7770D",
  "shared"        = "#555555"
)

module_headers <- tribble(
  ~name,             ~xmid,
  "TCR/Proximal",     1.7,
  "MAPK/NFAT",        4.9,
  "PI3K/AKT",         8.7,
  "NF-κB",           12.3,
  "JAK-STAT",        15.7,
  "Epi/Apoptosis",   19.8,
  "Mito/Metabolic",  23.8
) %>% mutate(col = MOD_COL[name])

# ── 10.  Compartment boxes ────────────────────────────────────
XMIN <- -0.4; XMAX <- 26.0
comp_boxes <- tribble(
  ~lab,              ~xmin_ov, ~xmax_ov, ~ymin,  ~ymax,  ~bg,
  "Extracellular",   XMIN,      21.6,    10.2,   11.0,  "#D6EAF8",
  "Plasma Membrane", XMIN,      21.6,     9.0,   10.15, "#D5F5E3",
  "Cytoplasm",       XMIN,      21.6,     1.2,    8.95, "#FDFEFE",
  "Nucleus",         XMIN,      21.6,    -0.55,   1.15, "#FEF9E7",
  "Mito/Metabolic",  21.8,      XMAX,    -0.55,   7.8,  "#FEF5E4"
) %>% mutate(xmin = xmin_ov, xmax = xmax_ov) %>%
  select(-xmin_ov, -xmax_ov)

mod_dividers <- tibble(
  xd   = c(3.4, 7.1, 10.6, 14.1, 17.6, 21.6),
  ymin = -0.5, ymax = 10.1)

# ── 11.  Draw helpers ─────────────────────────────────────────
arr_open  <- arrow(type = "open", length = unit(3.8, "pt"), angle = 22)
arr_blunt <- arrow(type = "open", length = unit(3.0, "pt"), angle = 90, ends = "last")

draw_segs <- function(df, col_v, lty_v, arr_v, lwd = 0.38, al = 0.80) {
  sub <- df %>% filter(ecol == col_v, elty == lty_v, !is_cross)
  if (nrow(sub) == 0) return(NULL)
  arr <- switch(arr_v, "open" = arr_open, "blunt" = arr_blunt, NULL)
  geom_segment(data = sub, aes(x = x1, y = y1, xend = x2, yend = y2),
    colour = col_v, linetype = lty_v, linewidth = lwd,
    alpha = al, arrow = arr, lineend = "round", show.legend = FALSE)
}

draw_arcs <- function(df, col_v, lty_v, arr_v, lwd = 0.32, curv = 0.22) {
  sub <- df %>% filter(ecol == col_v, elty == lty_v, is_cross)
  if (nrow(sub) == 0) return(NULL)
  arr <- switch(arr_v, "open" = arr_open, "blunt" = arr_blunt, NULL)
  geom_curve(data = sub, aes(x = x1, y = y1, xend = x2, yend = y2),
    colour = col_v, linetype = lty_v, linewidth = lwd,
    curvature = curv, alpha = 0.50, arrow = arr,
    lineend = "round", show.legend = FALSE)
}

draw_nodes <- function(df_sub, lwd_val) {
  rn <- df_sub %>% filter(shape %in% c("receptor","complex","output","ligand"))
  kn <- df_sub %>% filter(shape == "kinase")
  an <- df_sub %>% filter(shape == "adaptor")
  tn <- df_sub %>% filter(shape == "tf")
  list(
    if (nrow(rn) > 0)
      geom_tile(data = rn,
        aes(x = x, y = y, width = node_r * 2.3, height = node_r,
            fill = fill_col, colour = border_col),
        linewidth = lwd_val, show.legend = FALSE),
    if (nrow(kn) > 0)
      geom_regon(data = kn,
        aes(x0 = x, y0 = y, r = node_r, angle = pi/6,
            sides = 6, fill = fill_col, colour = border_col),
        linewidth = lwd_val, show.legend = FALSE),
    if (nrow(an) > 0)
      geom_ellipse(data = an,
        aes(x0 = x, y0 = y, a = node_r * 1.6, b = node_r * 0.72,
            angle = 0, fill = fill_col, colour = border_col),
        linewidth = lwd_val, show.legend = FALSE),
    if (nrow(tn) > 0)
      geom_regon(data = tn,
        aes(x0 = x, y0 = y, r = node_r * 1.22, angle = pi/4,
            sides = 4, fill = fill_col, colour = border_col),
        linewidth = lwd_val, show.legend = FALSE)
  )
}

# ── 12.  Assemble ggplot ──────────────────────────────────────
p_pathway <- ggplot() +

  # Compartment fills
  geom_rect(data = comp_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = lab),
    colour = "#CCCCCC", linewidth = 0.20, alpha = 0.38) +
  scale_fill_manual(name = "Compartment",
    values = setNames(comp_boxes$bg, comp_boxes$lab)) +
  new_scale_fill() +

  # Module dividers
  geom_segment(data = mod_dividers,
    aes(x = xd, xend = xd, y = ymin, yend = ymax),
    colour = "#C8C8C8", linewidth = 0.25, linetype = "dashed") +

  # Module headers
  annotate("rect",
    xmin = module_headers$xmid - 1.3,
    xmax = module_headers$xmid + 1.3,
    ymin = 10.55, ymax = 10.95,
    fill   = alpha(module_headers$col, 0.12),
    colour = module_headers$col, linewidth = 0.45) +
  annotate("text",
    x = module_headers$xmid, y = 10.75,
    label = module_headers$name,
    colour = module_headers$col, size = 2.55, fontface = "bold") +

  # ── Edges: intra-module segments ──
  draw_segs(edges_coord, "#1A5276", "solid",  "open")  +
  draw_segs(edges_coord, "#784212", "solid",  "open")  +
  draw_segs(edges_coord, "#7B241C", "solid",  "blunt") +
  draw_segs(edges_coord, "#1A5276", "dashed", "open")  +
  draw_segs(edges_coord, "#5D6D7E", "dashed", "none")  +
  draw_segs(edges_coord, "#6C3483", "dotted", "open")  +

  # ── Edges: cross-module arcs ──
  draw_arcs(edges_coord, "#1A5276", "solid",  "open")  +
  draw_arcs(edges_coord, "#784212", "solid",  "open")  +
  draw_arcs(edges_coord, "#7B241C", "solid",  "blunt") +
  draw_arcs(edges_coord, "#1A5276", "dashed", "open")  +

  # ── Connector (non-hit) nodes — small, grey, drawn first ──
  draw_nodes(nodes %>% filter(!is_hit), lwd_val = 0.35) +
  scale_fill_identity() + scale_colour_identity() +
  new_scale_fill() + new_scale_colour() +

  # ── Hit nodes — larger, coloured, drawn on top ──
  draw_nodes(nodes %>% filter(is_hit), lwd_val = 0.80) +
  scale_fill_identity() + scale_colour_identity() +
  new_scale_fill() + new_scale_colour() +

  # ── SH2-domain halo: orange ring around nodes with ≥1 SH2 domain ──
  ggforce::geom_circle(
    data = nodes %>% filter(has_sh2),
    aes(x0 = x, y0 = y, r = node_r + 0.10),
    fill = NA, colour = "#E67E22", linewidth = 0.65,
    linetype = "solid", show.legend = FALSE) +

  # ── Labels: hit nodes (bold, repel) ──
  geom_label_repel(
    data = nodes %>% filter(is_hit),
    aes(x = x, y = y, label = label, fontface = "bold"),
    size = 2.0, label.padding = unit(0.06, "lines"), label.size = 0.10,
    fill = alpha("white", 0.85), colour = "#111111",
    box.padding = 0.10, point.padding = 0.08,
    segment.linewidth = 0.18, segment.colour = "#BBBBBB",
    max.overlaps = 100, seed = 42, show.legend = FALSE) +

  # ── Labels: connector nodes (small, plain text, no box) ──
  geom_text(
    data = nodes %>% filter(!is_hit),
    aes(x = x, y = y, label = label),
    size = 1.50, colour = "#888888", fontface = "plain",
    vjust = -1.1, lineheight = 0.80) +

  # ── Dot-matrix: per-screen hit indicators (MERGE_MODE only) ──
  # One small dot per assay, centred below each hit node.
  # Dot fill = direction in that screen; position (L→R) = screen order.
  { if (!is.null(mg) && mg$n_assays >= 2) {
      dot_df <- build_dotmatrix(mg, nodes %>% filter(is_hit))
      list(
        geom_point(data = dot_df,
          aes(x = x_dot, y = y_dot, fill = dot_fill, colour = dot_border),
          shape = 21, size = 1.5, stroke = 0.35, show.legend = FALSE),
        scale_fill_identity(),
        scale_colour_identity(),
        new_scale_fill(),
        new_scale_colour()
      )
    } else NULL } +

  # ── Overflow strip ────────────────────────────────────────────
  { if (nrow(overflow_hits) > 0) {
      STRIP_Y  <- -1.85
      LABEL_Y  <- -2.35
      list(
        # Divider line
        annotate("segment",
          x = XMIN, xend = min(overflow$xi) + max(overflow$xi) + 1,
          y = STRIP_Y + 0.55, yend = STRIP_Y + 0.55,
          colour = "#CCCCCC", linewidth = 0.3),
        # Header
        annotate("text", x = XMIN + 0.1, y = STRIP_Y + 0.40,
          label = sprintf(
            "Additional hits not mapped to pathway backbone  (n=%d shown below)",
            nrow(overflow_hits)),
          size = 1.85, colour = "#444444", fontface = "bold", hjust = 0),
        # Group labels
        annotate("text",
          x = grp_lab$lx, y = LABEL_Y,
          label = grp_lab$grp, size = 1.65,
          colour = "#555555", fontface = "bold.italic", hjust = 0.5),
        # Tiles
        geom_tile(data = overflow,
          aes(x = xi, y = STRIP_Y, fill = fill_col, colour = border_col),
          width = 0.84, height = 0.38, linewidth = 0.40, show.legend = FALSE),
        scale_fill_identity(),
        scale_colour_identity(),
        new_scale_fill(), new_scale_colour(),
        # Gene names
        geom_text(data = overflow,
          aes(x = xi, y = STRIP_Y, label = gene),
          size = 1.30, colour = "white", fontface = "bold")
      )
    } else NULL } +

  # ── Legend: edge types + node key ──
  { yl <- -1.15; xl <- 0.0
    leg <- list(
      annotate("segment", x=xl,    xend=xl+0.6, y=yl,    yend=yl,
        colour="#1A5276", linewidth=0.6, arrow=arr_open),
      annotate("text",    x=xl+0.7,y=yl, label="Activation",     size=1.85, hjust=0),
      annotate("segment", x=xl,    xend=xl+0.6, y=yl-0.28, yend=yl-0.28,
        colour="#784212", linewidth=0.6, arrow=arr_open),
      annotate("text",    x=xl+0.7,y=yl-0.28,label="Phosphorylation",size=1.85,hjust=0),
      annotate("segment", x=xl+3.2,xend=xl+3.8,y=yl,    yend=yl,
        colour="#7B241C", linewidth=0.6, arrow=arr_blunt),
      annotate("text",    x=xl+3.9,y=yl, label="Inhibition",     size=1.85, hjust=0),
      annotate("segment", x=xl+3.2,xend=xl+3.8,y=yl-0.28,yend=yl-0.28,
        colour="#1A5276", linewidth=0.6, linetype="dashed"),
      annotate("text",    x=xl+3.9,y=yl-0.28,label="Indirect",   size=1.85, hjust=0),
      annotate("segment", x=xl+6.5,xend=xl+7.1,y=yl,    yend=yl,
        colour="#5D6D7E", linewidth=0.6, linetype="dashed"),
      annotate("text",    x=xl+7.2,y=yl, label="Binding",        size=1.85, hjust=0),
      annotate("segment", x=xl+6.5,xend=xl+7.1,y=yl-0.28,yend=yl-0.28,
        colour="#6C3483", linewidth=0.6, linetype="dotted", arrow=arr_open),
      annotate("text",    x=xl+7.2,y=yl-0.28,label="Translocation",size=1.85,hjust=0),
      # CRISPRi colour legend
      annotate("point",x=10.5,y=yl,      colour="#08306B",fill="#2166AC",size=3,shape=22,stroke=0.8),
      annotate("text", x=10.8,y=yl,
        label=sprintf("Depleted / GoF (n=%d)", nrow(hits_dep)),size=1.85,hjust=0),
      annotate("point",x=10.5,y=yl-0.28, colour="#67000D",fill="#D6604D",size=3,shape=22,stroke=0.8),
      annotate("text", x=10.8,y=yl-0.28,
        label=sprintf("Enriched / LoF (n=%d)", nrow(hits_enr)),size=1.85,hjust=0),
      annotate("point",x=15.5,y=yl,      colour="#909497",fill="#D5D8DC",size=2.2,shape=22,stroke=0.5),
      annotate("text", x=15.8,y=yl,      label="Backbone connector",size=1.85,hjust=0),
      # SH2-domain halo legend
      annotate("point",x=20.5,y=yl,      colour="#E67E22",fill="#FFFFFF",size=3.5,shape=21,stroke=1.5),
      annotate("text", x=20.8,y=yl,
        label=sprintf("SH2 domain (%d)", sum(nodes$has_sh2)),size=1.85,hjust=0),
      annotate("text", x=15.5,y=yl-0.28, size=1.70, hjust=0,
        colour="#666666", fontface="italic",
        label="⬡ Kinase  ◇ TF  ○ Adaptor  ▭ Receptor/Complex | size ∝ −log₁₀(Fisher p)")
    )
    # Multi-screen colour + dot-matrix legend (MERGE_MODE only)
    if (!is.null(mg) && mg$n_assays >= 2) {
      # Pull the same palette used by node_colours_merged()
      n_a  <- mg$n_assays
      apal <- ASSAY_PALETTE[((seq_len(n_a) - 1L) %% nrow(ASSAY_PALETTE)) + 1L, ]

      merge_leg <- list(
        # ── Header row ──────────────────────────────────────────────
        annotate("text", x=10.5, y=yl-0.56,
          label="Multi-screen node colours:", size=1.85,
          colour="#444444", hjust=0, fontface="bold"),
        # ── Both-screens (concordant) ────────────────────────────────
        annotate("point",x=10.5,y=yl-0.84,colour="#08306B",fill="#08519C",size=3,shape=22,stroke=0.8),
        annotate("text", x=10.8,y=yl-0.84,
          label=sprintf("Both screens GoF (%d)",
            sum(mg$hits$concordance=="concordant_depleted")),size=1.85,hjust=0),
        annotate("point",x=15.5,y=yl-0.84,colour="#67000D",fill="#CB181D",size=3,shape=22,stroke=0.8),
        annotate("text", x=15.8,y=yl-0.84,
          label=sprintf("Both screens LoF (%d)",
            sum(mg$hits$concordance=="concordant_enriched")),size=1.85,hjust=0),
        # ── Discordant ───────────────────────────────────────────────
        annotate("point",x=20.5,y=yl-0.84,colour="#9A7D0A",fill="#D4AC0D",size=3,shape=22,stroke=0.8),
        annotate("text", x=20.8,y=yl-0.84,
          label=sprintf("Discordant (%d)",
            sum(mg$hits$concordance=="discordant")),size=1.85,hjust=0)
      )

      # ── Per-assay rows (one row per screen) ─────────────────────────
      for (i in seq_len(n_a)) {
        aname   <- mg$assay_names[i]
        yrow    <- yl - 1.12 - (i - 1L) * 0.28
        n_gof_i <- sum(mg$hits[[paste0("hit_dir__", aname)]] == "depleted", na.rm = TRUE)
        n_lof_i <- sum(mg$hits[[paste0("hit_dir__", aname)]] == "enriched", na.rm = TRUE)
        merge_leg <- c(merge_leg, list(
          # Screen name label
          annotate("text", x=10.5, y=yrow,
            label=sprintf("%s only:", aname), size=1.75,
            colour="#333333", hjust=0, fontface="italic"),
          # GoF patch
          annotate("point",x=14.8,y=yrow,
            colour=apal$gof_border[i], fill=apal$gof_fill[i],
            size=3, shape=22, stroke=0.8),
          annotate("text", x=15.1,y=yrow,
            label=sprintf("GoF (%d)", n_gof_i), size=1.75, hjust=0),
          # LoF patch
          annotate("point",x=18.0,y=yrow,
            colour=apal$lof_border[i], fill=apal$lof_fill[i],
            size=3, shape=22, stroke=0.8),
          annotate("text", x=18.3,y=yrow,
            label=sprintf("LoF (%d)", n_lof_i), size=1.75, hjust=0)
        ))
      }

      # ── Dot-matrix key ───────────────────────────────────────────────
      dot_leg_y <- yl - 1.12 - n_a * 0.28 - 0.15
      merge_leg <- c(merge_leg,
        dotmatrix_legend_items(mg, xl = 0.0, yl = dot_leg_y))

      leg <- c(leg, merge_leg)
    }
    leg
  } +

  coord_cartesian(
    xlim = c(XMIN, XMAX),
    ylim = c(
      # Extra vertical space: 0.28 per assay row in the legend
      if (!is.null(mg) && mg$n_assays >= 2)
        -1.15 - 1.12 - mg$n_assays * 0.28 - 0.70   # legend bottom
      else if (nrow(overflow_hits) > 0) -2.55
      else                              -1.55,
      11.15),
    clip = "off") +
  labs(
    title = sprintf("Integrated T Cell Signalling — CRISPRi Screen (%s)", MERGE_LABEL),
    subtitle = sprintf(
      "%d GoF (blue) · %d LoF (red) · %d on pathway · %d in strip%s\nTop ORA: %s",
      nrow(hits_dep), nrow(hits_enr),
      nrow(placed),
      if (nrow(overflow_hits) > 0) nrow(overflow) else 0L,
      if (!is.null(mg) && mg$n_assays >= 2)
        sprintf(" · %d concordant · %d single-screen · %d discordant",
                sum(grepl("^concordant", mg$hits$concordance)),
                sum(mg$hits$concordance == "single_assay"),
                sum(mg$hits$concordance == "discordant"))
      else "",
      str_trunc(ora_sub, 80)),
    caption = sprintf(
      "Depleted = GoF · Enriched = LoF%s   |   Generated: %s",
      if (!is.null(mg) && mg$n_assays >= 2)
        sprintf(" · Node size ∝ −log₁₀(Fisher combined p) across %d screens",
                mg$n_assays)
      else " · Node size ∝ −log₁₀FDR",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.title    = element_text(size = 11.5, face = "bold",
                                 hjust = 0.5, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 7.0, colour = "#444444",
                                 hjust = 0.5, margin = margin(b = 1)),
    plot.caption  = element_text(size = 6.2, colour = "#888888",
                                 hjust = 0.5, margin = margin(t = 2)),
    plot.margin   = margin(8, 6, 14, 6),
    legend.position = "right",
    legend.title  = element_text(size = 8, face = "bold"),
    legend.text   = element_text(size = 7)
  )

# ── 13.  Export ───────────────────────────────────────────────
cat("\nSaving...\n")
OUT_TS <- paste0(OUT, "_", TSTAMP)   # timestamped base path
ggsave(paste0(OUT_TS, ".pdf"), p_pathway, width = 32, height = 16,
       units = "in", device = cairo_pdf)
ggsave(paste0(OUT_TS, ".png"), p_pathway, width = 32, height = 16,
       units = "in", dpi = 170)
cat(sprintf("  %s.pdf\n  %s.png\n", OUT_TS, OUT_TS))

# Overlap bar chart (MERGE_MODE only)
if (!is.null(mg)) {
  p_overlap <- plot_overlap(mg)
  ggsave(paste0(OUT_TS, "_overlap.pdf"), p_overlap,
         width = 7, height = 5, units = "in", device = cairo_pdf)
  cat(sprintf("  %s_overlap.pdf\n", OUT_TS))
}

# ORA dotplot
if (nrow(ora_df) >= 3) {
  ora_top <- ora_df %>%
    slice_head(n = 20) %>%
    mutate(Description = str_trunc(Description, 58),
           Description = fct_reorder(Description, -p.adjust))
  p_ora <- ggplot(ora_top,
    aes(x = GeneRatio, y = Description, size = Count, colour = p.adjust)) +
    geom_point() +
    scale_colour_gradient(low = "#C0392B", high = "#AEB6BF", name = "FDR") +
    scale_size_continuous(range = c(2, 8), name = "Gene count") +
    labs(title = "Reactome ORA — all CRISPRi hits",
         x = "Gene Ratio", y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8))
  ggsave(paste0(OUT_TS, "_ORA.pdf"), p_ora,
         width = 10, height = 7, units = "in", device = cairo_pdf)
  cat(sprintf("  %s_ORA.pdf\n", OUT_TS))
}

# Full hit table
hit_tbl <- bind_rows(
  placed %>%
    left_join(nodes %>% distinct(gene, module, shape), by = "gene") %>%
    mutate(location = "pathway"),
  if (nrow(overflow_hits) > 0)
    overflow %>%
      select(gene, hit_dir, best_lfc, best_fdr, nlp, grp) %>%
      rename(module = grp) %>%
      mutate(location = "strip", shape = NA_character_)
  else NULL
) %>%
  arrange(location, hit_dir, best_fdr) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
write.csv(hit_tbl, paste0(OUT_TS, "_hits_table.csv"), row.names = FALSE)
cat(sprintf("  %s_hits_table.csv  (%d rows)\n", OUT_TS, nrow(hit_tbl)))

cat("\n── Summary ──\n")
cat(sprintf("On-pathway : %d / %d  (%.0f%%)\n",
            nrow(placed), nrow(hits), 100 * nrow(placed) / nrow(hits)))
cat(sprintf("Strip      : %d\n", nrow(overflow_hits)))
cat(sprintf("Total shown: %d / %d  (%.0f%%)\n",
            nrow(placed) + nrow(overflow_hits), nrow(hits),
            100 * (nrow(placed) + nrow(overflow_hits)) / nrow(hits)))
