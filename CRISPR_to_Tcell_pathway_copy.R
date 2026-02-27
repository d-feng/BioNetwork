# ============================================================
# CRISPR → Integrated T Cell Pathway Diagram  (v4 – hit-centric)
# ============================================================
# Design philosophy:
#   • EVERY significant CRISPRi hit is shown on the diagram
#   • Non-hit backbone nodes appear ONLY as small grey connectors
#     when needed to link two hits
#   • Hits are placed by their biological role (module + y-band)
#   • Node size ∝ significance (−log10 FDR); colour = direction
#   • Hits with no backbone position → overflow panel (sorted by LFC)
#
# Layout (x-bands):
#   0 – 3.5   TCR / proximal signalling
#   3.5 – 7   MAPK / RAS / NFAT / PLCγ
#   7 – 10.5  PI3K / AKT / mTOR / cell cycle
#  10.5 – 14  NF-κB / ubiquitin / transcription
#  14 – 17.5  JAK-STAT / cytokine
#  17.5 – 21  Epigenetic / chromatin / apoptosis
#
# y-bands (shared across all modules):
#   9.5 – 10  plasma membrane receptors
#   6 – 9     cytoplasm – proximal signalling
#   3 – 6     cytoplasm – distal / effectors
#   1 – 3     cytoplasm – transcriptional regulators
#   0 – 1     nucleus – TFs
# ============================================================

# ── 0. Packages ──────────────────────────────────────────────
needed_bioc <- c("clusterProfiler","ReactomePA","org.Hs.eg.db")
needed_cran <- c("ggplot2","ggforce","ggrepel","ggnewscale",
                 "dplyr","tibble","tidyr","scales","grid",
                 "patchwork","forcats","stringr","purrr","RColorBrewer")
new_b <- needed_bioc[!sapply(needed_bioc, requireNamespace, quietly=TRUE)]
if (length(new_b)) BiocManager::install(new_b)
new_c <- needed_cran[!sapply(needed_cran, requireNamespace, quietly=TRUE)]
if (length(new_c)) install.packages(new_c)

library(clusterProfiler); library(ReactomePA); library(org.Hs.eg.db)
library(ggplot2); library(ggforce); library(ggrepel); library(ggnewscale)
library(dplyr); library(tibble); library(tidyr); library(scales)
library(grid); library(patchwork); library(forcats); library(stringr); library(purrr)

# Pin verbs
select        <- dplyr::select;   filter      <- dplyr::filter
rename        <- dplyr::rename;   mutate      <- dplyr::mutate
arrange       <- dplyr::arrange;  summarise   <- dplyr::summarise
group_by      <- dplyr::group_by; ungroup     <- dplyr::ungroup
left_join     <- dplyr::left_join; bind_rows  <- dplyr::bind_rows
distinct      <- dplyr::distinct; pull        <- dplyr::pull
slice_head    <- dplyr::slice_head; count     <- dplyr::count
transmute     <- dplyr::transmute; case_when  <- dplyr::case_when
anti_join     <- dplyr::anti_join; inner_join <- dplyr::inner_join
replace_na    <- tidyr::replace_na; separate_rows <- tidyr::separate_rows
fct_reorder   <- forcats::fct_reorder; str_trunc  <- stringr::str_trunc
map_chr       <- purrr::map_chr

# ── 1. Settings ──────────────────────────────────────────────
INPUT_FILE  <- "C:/Users/difen/Downloads/Carnevale_CRISPRi_Tcell_mageck_gene_summary.txt"
FDR_CUTOFF  <- 0.10
LFC_CUTOFF  <- 0.20
ORA_FDR     <- 0.10
OUT_PREFIX  <- "C:/Users/difen/Rcode/Tcell_CRISPRi_pathway"

# ── 2. MAGeCK data ───────────────────────────────────────────
raw <- read.delim(INPUT_FILE, check.names=FALSE, stringsAsFactors=FALSE)
cat(sprintf("Loaded %d genes\n", nrow(raw)))

crispr_all <- raw %>%
  transmute(
    gene    = id,
    neg_lfc = neg.lfc, neg_fdr = neg.fdr,
    pos_lfc = pos.lfc, pos_fdr = pos.fdr,
    hit_dir = case_when(
      neg.fdr < FDR_CUTOFF & neg.lfc < -LFC_CUTOFF ~ "depleted",
      pos.fdr < FDR_CUTOFF & pos.lfc >  LFC_CUTOFF ~ "enriched",
      TRUE ~ "none"),
    best_lfc = case_when(
      hit_dir=="depleted" ~ neg.lfc,
      hit_dir=="enriched" ~ pos.lfc,
      TRUE ~ 0),
    best_fdr  = pmin(neg.fdr, pos.fdr),
    nlp_fdr   = -log10(pmax(pmin(neg.fdr, pos.fdr), 1e-6))
  )

hits <- crispr_all %>% filter(hit_dir != "none") %>%
  arrange(hit_dir, best_fdr, best_lfc)
hits_dep <- hits %>% filter(hit_dir=="depleted")
hits_enr <- hits %>% filter(hit_dir=="enriched")
cat(sprintf("Hits: %d depleted (GoF/blue), %d enriched (LoF/red)\n",
            nrow(hits_dep), nrow(hits_enr)))

# ── 3. Reactome ORA ──────────────────────────────────────────
entrez_map <- tryCatch(
  bitr(hits$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db),
  error=function(e) data.frame(SYMBOL=character(), ENTREZID=character()))

ora <- tryCatch(
  enrichPathway(gene=entrez_map$ENTREZID, organism="human",
                pvalueCutoff=ORA_FDR, qvalueCutoff=0.25, readable=TRUE),
  error=function(e) NULL)

if (!is.null(ora) && nrow(as.data.frame(ora)) > 0) {
  ora_df <- as.data.frame(ora) %>% arrange(p.adjust) %>% slice_head(n=20)
  cat(sprintf("ORA: %d enriched Reactome pathways\n", nrow(ora_df)))
  ora_subtitle <- paste(head(str_trunc(ora_df$Description, 50), 3), collapse=" · ")
} else {
  ora_df <- data.frame()
  ora_subtitle <- "ORA: no enriched pathways at FDR<0.10"
}

# ── 4. Hit placement map ──────────────────────────────────────
# Each hit is assigned: module (x-band), role_y (y within band),
# label, node type.
# Non-assigned hits go to overflow strip.
#
# x-band centres: TCR=1.7, MAPK=5.2, PI3K=8.8, NFkB=12.2, JAK=15.8, Epi=19.2
# y values: membrane=9.7, proximal=7.5–8.5, mid=5–7, distal=3–5, nucleus_tf=0.5–2

hit_positions <- tribble(
  ~gene,      ~x,    ~y,    ~module,         ~node_type,  ~label,
  # ── TCR / Proximal (x 0.5–3.3) ─────────────────────────────────────────────
  # Membrane receptors
  "CD247",    1.0,   9.7,  "TCR Proximal",  "receptor",  "CD3ζ",
  "CD3D",     1.5,   9.7,  "TCR Proximal",  "receptor",  "CD3δ",
  "CD3E",     2.0,   9.7,  "TCR Proximal",  "receptor",  "CD3ε",
  "CD28",     2.7,   9.3,  "TCR Proximal",  "receptor",  "CD28",
  "PDCD1",    0.4,   9.3,  "TCR Proximal",  "receptor",  "PD-1",
  "LAG3",     0.9,   9.0,  "TCR Proximal",  "receptor",  "LAG3",
  "CD5",      2.5,   9.0,  "TCR Proximal",  "receptor",  "CD5",
  "CTLA4",    3.1,   9.3,  "TCR Proximal",  "receptor",  "CTLA4",
  "CD8A",     3.2,   9.7,  "TCR Proximal",  "receptor",  "CD8A",
  # Proximal kinases
  "LCK",      1.8,   8.3,  "TCR Proximal",  "kinase",    "LCK",
  "FYN",      0.8,   8.3,  "TCR Proximal",  "kinase",    "FYN",
  "ZAP70",    1.3,   7.5,  "TCR Proximal",  "kinase",    "ZAP70",
  "ITK",      0.5,   6.8,  "TCR Proximal",  "kinase",    "ITK",
  # Adaptors / scaffold
  "LCP2",     0.7,   6.0,  "TCR Proximal",  "adaptor",   "SLP-76",
  "LAT",      1.5,   6.5,  "TCR Proximal",  "adaptor",   "LAT",
  "GRAP2",    2.2,   6.8,  "TCR Proximal",  "adaptor",   "GRAP2",
  "GRB2",     2.5,   6.2,  "TCR Proximal",  "adaptor",   "GRB2",
  "VAV1",     0.4,   5.4,  "TCR Proximal",  "adaptor",   "VAV1",
  "SH2D1A",   1.0,   5.4,  "TCR Proximal",  "adaptor",   "SAP",
  # Negative regulators (LoF hits) – proximal
  "CBLB",     2.0,   7.5,  "TCR Proximal",  "adaptor",   "CBLB",
  "UBASH3A",  0.2,   7.1,  "TCR Proximal",  "adaptor",   "UBASH3A",
  "PTPN6",    2.8,   7.8,  "TCR Proximal",  "adaptor",   "SHP-1",
  # ── MAPK / RAS / NFAT (x 3.6–6.8) ─────────────────────────────────────────
  "PLCG1",    3.8,   6.5,  "MAPK/NFAT",    "kinase",    "PLCγ1",
  "RASGRP1",  3.8,   5.7,  "MAPK/NFAT",    "adaptor",   "RasGRP1",
  "KRAS",     3.8,   5.0,  "MAPK/NFAT",    "kinase",    "RAS",
  "RAF1",     3.8,   4.2,  "MAPK/NFAT",    "kinase",    "RAF",
  "MAP2K1",   3.8,   3.5,  "MAPK/NFAT",    "kinase",    "MEK1/2",
  "MAPK3",    3.8,   2.7,  "MAPK/NFAT",    "kinase",    "ERK1/2",
  "MAPK1",    3.8,   0.8,  "MAPK/NFAT",    "tf",        "ERK(nuc)",
  "PPP3CA",   5.5,   5.5,  "MAPK/NFAT",    "complex",   "Calcineurin",
  "NFATC1",   5.2,   0.8,  "MAPK/NFAT",    "tf",        "NFAT",
  # MAPK negative regulators (LoF hits)
  "RASA2",    5.0,   5.0,  "MAPK/NFAT",    "adaptor",   "RASA2",
  "DUSP4",    5.0,   2.7,  "MAPK/NFAT",    "adaptor",   "DUSP4",
  # ── PI3K / AKT / mTOR / Cell cycle (x 7.0–10.3) ────────────────────────────
  "PIK3CD",   7.2,   6.5,  "PI3K/AKT",     "kinase",    "PI3Kδ",
  "PIK3R1",   7.8,   7.2,  "PI3K/AKT",     "adaptor",   "p85α",
  "PTEN",     8.5,   6.5,  "PI3K/AKT",     "adaptor",   "PTEN",
  "PDK1",     7.2,   5.7,  "PI3K/AKT",     "kinase",    "PDK1",
  "AKT1",     7.2,   4.9,  "PI3K/AKT",     "kinase",    "AKT",
  "MTOR",     7.2,   4.1,  "PI3K/AKT",     "complex",   "mTOR",
  "RPS6KB1",  6.7,   3.3,  "PI3K/AKT",     "kinase",    "S6K1",
  "EIF4EBP1", 7.7,   3.3,  "PI3K/AKT",     "adaptor",   "4EBP1",
  "FOXO1",    7.2,   0.8,  "PI3K/AKT",     "tf",        "FOXO1",
  "CDKN1B",   8.5,   0.8,  "PI3K/AKT",     "tf",        "p27",
  # DGK negative regulators (LoF hits)
  "DGKZ",     8.5,   5.7,  "PI3K/AKT",     "kinase",    "DGKζ",
  "DGKA",     9.2,   5.0,  "PI3K/AKT",     "kinase",    "DGKα",
  # ── NF-κB / Ubiquitin / Transcription (x 10.5–14.0) ───────────────────────
  "PRKCQ",   10.7,   6.5,  "NF-κB",        "kinase",    "PKCθ",
  "CARD11",  11.2,   5.7,  "NF-κB",        "complex",   "CARMA1",
  "BCL10",   11.7,   4.9,  "NF-κB",        "adaptor",   "BCL10",
  "MALT1",   12.3,   4.9,  "NF-κB",        "adaptor",   "MALT1",
  "IKBKB",   11.5,   4.1,  "NF-κB",        "kinase",    "IKKβ",
  "CHUK",    12.3,   4.1,  "NF-κB",        "kinase",    "IKKα",
  "NFKBIA",  11.9,   3.3,  "NF-κB",        "adaptor",   "IκBα",
  "RELA",    11.4,   0.8,  "NF-κB",        "tf",        "p65",
  "NFKB1",   12.2,   0.8,  "NF-κB",        "tf",        "p50",
  "IRF4",    10.6,   0.8,  "NF-κB",        "tf",        "IRF4",
  "GATA3",   13.0,   0.8,  "NF-κB",        "tf",        "GATA3",
  "JUNB",    13.5,   0.8,  "NF-κB",        "tf",        "JUNB",
  # NF-κB LoF hits
  "TNFAIP3", 13.0,   4.9,  "NF-κB",        "adaptor",   "A20",
  "TCEB2",   13.5,   4.1,  "NF-κB",        "adaptor",   "Elongin-B",
  "ARIH2",   13.5,   3.3,  "NF-κB",        "adaptor",   "ARIH2",
  # ── JAK-STAT / Cytokine (x 14.2–17.3) ─────────────────────────────────────
  "IL2RA",   14.5,   9.7,  "JAK-STAT",     "receptor",  "IL-2Rα",
  "IL2RB",   15.2,   9.7,  "JAK-STAT",     "receptor",  "IL-2Rβ",
  "IL2RG",   15.8,   9.7,  "JAK-STAT",     "receptor",  "γc",
  "IFNGR1",  16.8,   9.7,  "JAK-STAT",     "receptor",  "IFNγR1",
  "JAK1",    14.8,   8.5,  "JAK-STAT",     "kinase",    "JAK1",
  "JAK3",    15.5,   8.5,  "JAK-STAT",     "kinase",    "JAK3",
  "JAK2",    16.8,   8.5,  "JAK-STAT",     "kinase",    "JAK2",
  "STAT5A",  14.8,   7.5,  "JAK-STAT",     "tf",        "STAT5A",
  "STAT5B",  15.5,   7.5,  "JAK-STAT",     "tf",        "STAT5B",
  "STAT1",   16.8,   7.5,  "JAK-STAT",     "tf",        "STAT1",
  "STAT5n",  15.2,   0.8,  "JAK-STAT",     "tf",        "STAT5(n)",
  "STAT1n",  16.8,   0.8,  "JAK-STAT",     "tf",        "STAT1(n)",
  # JAK-STAT LoF hits
  "SOCS1",   14.2,   6.5,  "JAK-STAT",     "adaptor",   "SOCS1",
  "IRF2BP2", 17.2,   0.8,  "JAK-STAT",     "tf",        "IRF2BP2",
  # ── Epigenetic / Chromatin / Apoptosis (x 17.5–20.8) ──────────────────────
  "EZH2",    18.0,   2.2,  "Epi/Apoptosis","complex",   "EZH2",
  "EED",     18.6,   2.9,  "Epi/Apoptosis","adaptor",   "EED",
  "SUZ12",   17.5,   2.9,  "Epi/Apoptosis","adaptor",   "SUZ12",
  "DNMT3A",  19.2,   2.2,  "Epi/Apoptosis","adaptor",   "DNMT3A",
  "TET2",    19.8,   2.2,  "Epi/Apoptosis","adaptor",   "TET2",
  "KDM6A",   18.0,   1.5,  "Epi/Apoptosis","adaptor",   "UTX",
  "ARID1A",  19.2,   1.5,  "Epi/Apoptosis","adaptor",   "ARID1A",
  "SMARCB1", 19.8,   1.5,  "Epi/Apoptosis","adaptor",   "SMARCB1",
  "MYC",     18.6,   0.8,  "Epi/Apoptosis","tf",        "MYC",
  "MEF2D",   17.5,   0.8,  "Epi/Apoptosis","tf",        "MEF2D",
  "TP53",    20.4,   0.8,  "Epi/Apoptosis","tf",        "p53",
  "BCL2",    18.0,   4.5,  "Epi/Apoptosis","adaptor",   "BCL2",
  "BCL2L1",  18.8,   5.2,  "Epi/Apoptosis","adaptor",   "BCL-XL",
  "BAX",     19.6,   5.9,  "Epi/Apoptosis","adaptor",   "BAX",
  "MCL1",    17.5,   5.2,  "Epi/Apoptosis","adaptor",   "MCL1",
  "CASP3",   19.0,   3.7,  "Epi/Apoptosis","kinase",    "Casp-3",
  # Epi LoF hits
  "TMEM222", 20.4,   2.9,  "Epi/Apoptosis","adaptor",   "TMEM222",
  # ── Transcription factors (placed in NF-κB / Epi zones) ────────────────────
  "YBX1",    13.8,   1.5,  "NF-κB",        "tf",        "YBX1",
  "HIVEP2",  10.6,   1.5,  "NF-κB",        "tf",        "HIVEP2",
  "IKZF2",   11.5,   1.5,  "NF-κB",        "tf",        "IKZF2",
  "NR4A1",   12.3,   1.5,  "NF-κB",        "tf",        "NR4A1",
  "TOB1",    13.0,   1.5,  "NF-κB",        "tf",        "TOB1",
  "EOMES",   17.5,   1.5,  "Epi/Apoptosis","tf",        "EOMES",
  "TBX21",   18.6,   1.5,  "Epi/Apoptosis","tf",        "T-bet",
  "TCF7",    19.2,   0.8,  "Epi/Apoptosis","tf",        "TCF1",
  "PRDM1",   20.0,   1.5,  "Epi/Apoptosis","tf",        "BLIMP1",
  "BATF",    20.6,   1.5,  "Epi/Apoptosis","tf",        "BATF",
  "BACH2",   20.0,   0.8,  "Epi/Apoptosis","tf",        "BACH2",
  "SUPT4H1", 13.8,   3.3,  "NF-κB",        "adaptor",   "SUPT4H1",
  "TAF6",    14.0,   4.1,  "NF-κB",        "adaptor",   "TAF6"
)

# ── 5. Minimal backbone connectors (non-hit grey nodes) ───────
# These appear ONLY to provide biological context between hits.
# They are small, grey, and labelled in smaller text.
connector_nodes <- tribble(
  ~gene,    ~x,    ~y,    ~module,        ~node_type,  ~label,
  "TCR",    1.5,  10.1,  "TCR Proximal", "complex",   "TCRαβ",
  "pMHCII", 1.5,  10.7,  "TCR Proximal", "ligand",    "pMHC-II",
  "CD4",    2.2,   9.7,  "TCR Proximal", "receptor",  "CD4",
  "GRB2",   2.5,   6.2,  "TCR Proximal", "adaptor",   "GRB2",
  "NFKBIA", 11.9,  3.3,  "NF-κB",       "adaptor",   "IκBα",
  "STAT5n", 15.2,  0.8,  "JAK-STAT",    "tf",        "STAT5(n)",
  "STAT1n", 16.8,  0.8,  "JAK-STAT",    "tf",        "STAT1(n)",
  "gene_exp",10.5,-0.6,  "shared",      "output",    "Gene\nExpression"
)

# ── 6. Backbone edges ─────────────────────────────────────────
backbone_edges <- tribble(
  ~from,     ~to,        ~type,            ~is_cross,
  # pMHC → TCR complex
  "pMHCII",  "TCR",      "binding",        FALSE,
  "TCR",     "CD247",    "binding",        FALSE,
  "TCR",     "CD3D",     "binding",        FALSE,
  "TCR",     "CD3E",     "binding",        FALSE,
  "CD4",     "LCK",      "activation",     FALSE,
  "CD28",    "LCK",      "activation",     FALSE,
  # Checkpoint inhibition
  "PDCD1",   "ZAP70",    "inhibition",     FALSE,
  "LAG3",    "LCK",      "inhibition",     FALSE,
  "CTLA4",   "LCK",      "inhibition",     FALSE,
  "CD5",     "LCK",      "inhibition",     FALSE,
  "UBASH3A", "ZAP70",    "inhibition",     FALSE,
  "PTPN6",   "LCK",      "inhibition",     FALSE,
  "CBLB",    "ZAP70",    "inhibition",     FALSE,
  "CBLB",    "LAT",      "inhibition",     FALSE,
  # Proximal kinases
  "LCK",     "ZAP70",    "phosphorylation",FALSE,
  "FYN",     "ZAP70",    "phosphorylation",FALSE,
  "ZAP70",   "LAT",      "phosphorylation",FALSE,
  "ZAP70",   "LCP2",     "phosphorylation",FALSE,
  # Adaptors
  "LAT",     "PLCG1",    "activation",     FALSE,
  "LAT",     "GRB2",     "binding",        FALSE,
  "LAT",     "ITK",      "activation",     FALSE,
  "LAT",     "GRAP2",    "binding",        FALSE,
  "LCP2",    "VAV1",     "activation",     FALSE,
  "LCP2",    "ITK",      "binding",        FALSE,
  "LCP2",    "RASGRP1",  "activation",     FALSE,
  "ITK",     "PLCG1",    "phosphorylation",FALSE,
  "SH2D1A",  "LCK",      "activation",     FALSE,
  # MAPK cascade
  "PLCG1",   "RASGRP1",  "activation",     FALSE,
  "PLCG1",   "PPP3CA",   "indirect_act",   FALSE,
  "RASGRP1", "KRAS",     "activation",     FALSE,
  "RASA2",   "KRAS",     "inhibition",     FALSE,
  "KRAS",    "RAF1",     "activation",     FALSE,
  "RAF1",    "MAP2K1",   "phosphorylation",FALSE,
  "MAP2K1",  "MAPK3",    "phosphorylation",FALSE,
  "MAPK3",   "MAPK1",    "translocation",  FALSE,
  "DUSP4",   "MAPK3",    "inhibition",     FALSE,
  "PPP3CA",  "NFATC1",   "activation",     FALSE,
  "NFATC1",  "gene_exp", "activation",     FALSE,
  "MAPK1",   "gene_exp", "activation",     FALSE,
  # PI3K branch
  "CD28",    "PIK3CD",   "activation",     TRUE,
  "VAV1",    "PIK3CD",   "activation",     TRUE,
  "CBLB",    "PIK3CD",   "inhibition",     TRUE,
  "PIK3R1",  "PIK3CD",   "binding",        FALSE,
  "PTEN",    "PIK3CD",   "inhibition",     FALSE,
  "DGKZ",    "PIK3CD",   "inhibition",     FALSE,
  "DGKA",    "PIK3CD",   "inhibition",     FALSE,
  "PIK3CD",  "PDK1",     "indirect_act",   FALSE,
  "PDK1",    "AKT1",     "phosphorylation",FALSE,
  "AKT1",    "MTOR",     "activation",     FALSE,
  "AKT1",    "FOXO1",    "phosphorylation",FALSE,
  "AKT1",    "CDKN1B",   "phosphorylation",FALSE,
  "AKT1",    "BCL2",     "activation",     TRUE,
  "MTOR",    "RPS6KB1",  "phosphorylation",FALSE,
  "MTOR",    "EIF4EBP1", "phosphorylation",FALSE,
  "MTOR",    "MYC",      "indirect_act",   TRUE,
  "FOXO1",   "gene_exp", "inhibition",     FALSE,
  "CDKN1B",  "gene_exp", "inhibition",     FALSE,
  # NF-κB
  "PLCG1",   "PRKCQ",    "activation",     TRUE,
  "PRKCQ",   "CARD11",   "phosphorylation",FALSE,
  "CARD11",  "BCL10",    "binding",        FALSE,
  "BCL10",   "MALT1",    "binding",        FALSE,
  "MALT1",   "IKBKB",    "activation",     FALSE,
  "TNFAIP3", "MALT1",    "inhibition",     FALSE,
  "TNFAIP3", "IKBKB",    "inhibition",     FALSE,
  "IKBKB",   "NFKBIA",   "phosphorylation",FALSE,
  "CHUK",    "NFKBIA",   "phosphorylation",FALSE,
  "NFKBIA",  "RELA",     "inhibition",     FALSE,
  "NFKBIA",  "NFKB1",    "inhibition",     FALSE,
  "RELA",    "gene_exp", "activation",     FALSE,
  "NFKB1",   "gene_exp", "activation",     FALSE,
  "IRF4",    "gene_exp", "activation",     FALSE,
  "GATA3",   "gene_exp", "activation",     FALSE,
  "JUNB",    "gene_exp", "activation",     FALSE,
  "TCEB2",   "IKBKB",    "inhibition",     FALSE,
  "ARIH2",   "NFKBIA",   "indirect_act",   FALSE,
  # JAK-STAT
  "IL2RB",   "JAK1",     "activation",     FALSE,
  "IL2RG",   "JAK3",     "activation",     FALSE,
  "IFNGR1",  "JAK2",     "activation",     FALSE,
  "JAK1",    "STAT5A",   "phosphorylation",FALSE,
  "JAK3",    "STAT5B",   "phosphorylation",FALSE,
  "JAK2",    "STAT1",    "phosphorylation",FALSE,
  "SOCS1",   "JAK1",     "inhibition",     FALSE,
  "SOCS1",   "JAK3",     "inhibition",     FALSE,
  "SOCS1",   "AKT1",     "inhibition",     TRUE,
  "STAT5A",  "STAT5n",   "translocation",  FALSE,
  "STAT5B",  "STAT5n",   "translocation",  FALSE,
  "STAT5n",  "gene_exp", "activation",     FALSE,
  "STAT1",   "STAT1n",   "translocation",  FALSE,
  "STAT1n",  "gene_exp", "activation",     FALSE,
  "IRF2BP2", "gene_exp", "inhibition",     FALSE,
  # Epigenetic
  "EED",     "EZH2",     "binding",        FALSE,
  "SUZ12",   "EZH2",     "binding",        FALSE,
  "KDM6A",   "EZH2",     "inhibition",     FALSE,
  "ARID1A",  "EZH2",     "inhibition",     FALSE,
  "SMARCB1", "EZH2",     "inhibition",     FALSE,
  "TET2",    "DNMT3A",   "inhibition",     FALSE,
  "EZH2",    "gene_exp", "inhibition",     FALSE,
  "DNMT3A",  "gene_exp", "inhibition",     FALSE,
  "MYC",     "gene_exp", "activation",     FALSE,
  "MEF2D",   "gene_exp", "activation",     FALSE,
  "BCL2",    "BAX",      "inhibition",     FALSE,
  "BCL2L1",  "BAX",      "inhibition",     FALSE,
  "MCL1",    "BAX",      "inhibition",     FALSE,
  "BAX",     "CASP3",    "activation",     FALSE,
  "TP53",    "BAX",      "activation",     FALSE,
  "TP53",    "BCL2",     "inhibition",     FALSE,
  "YBX1",    "gene_exp", "activation",     FALSE,
  "HIVEP2",  "gene_exp", "activation",     FALSE,
  "IKZF2",   "gene_exp", "inhibition",     FALSE,
  "NR4A1",   "gene_exp", "activation",     FALSE,
  "TOB1",    "gene_exp", "inhibition",     FALSE,
  "EOMES",   "gene_exp", "activation",     FALSE,
  "TBX21",   "gene_exp", "activation",     FALSE,
  "TCF7",    "gene_exp", "activation",     FALSE,
  "PRDM1",   "gene_exp", "inhibition",     FALSE,
  "BATF",    "gene_exp", "activation",     FALSE,
  "BACH2",   "gene_exp", "inhibition",     FALSE,
  "SUPT4H1", "gene_exp", "activation",     FALSE,
  "TAF6",    "gene_exp", "activation",     FALSE
)

# ── 7. Build unified node table ───────────────────────────────
# All nodes: hit_positions + connector_nodes (deduped)
# Join with CRISPRi data

# Check which connector nodes are already hits (assign hit aesthetics there)
all_positioned <- bind_rows(
  hit_positions %>% mutate(is_connector = FALSE),
  connector_nodes %>% filter(!gene %in% hit_positions$gene) %>%
    mutate(is_connector = TRUE)
)

nodes_plot <- all_positioned %>%
  left_join(crispr_all %>% select(gene, hit_dir, best_lfc, best_fdr, nlp_fdr),
            by = "gene") %>%
  mutate(
    hit_dir    = replace_na(hit_dir, "none"),
    nlp_fdr    = replace_na(nlp_fdr, 0),
    fill_col   = case_when(
      hit_dir == "depleted" ~ "#2166AC",
      hit_dir == "enriched" ~ "#D6604D",
      is_connector           ~ "#D5D8DC",
      TRUE                   ~ "#BDC3C7"),
    border_col = case_when(
      hit_dir == "depleted" ~ "#053061",
      hit_dir == "enriched" ~ "#67000D",
      TRUE                   ~ "#909497"),
    # Hits: size scales with significance; non-hits: small fixed size
    node_r     = case_when(
      hit_dir == "none" & is_connector ~ 0.18,
      hit_dir == "none"               ~ 0.20,
      TRUE ~ rescale(pmin(nlp_fdr, 6), to=c(0.26, 0.50), from=c(0, 6))),
    label_size = case_when(
      hit_dir != "none"  ~ 2.2,
      is_connector       ~ 1.6,
      TRUE               ~ 1.8),
    label_face = if_else(hit_dir != "none", "bold", "plain"),
    label_col  = if_else(hit_dir != "none", "#111111", "#666666")
  )

# How many hits are now on diagram?
hits_on_diagram <- nodes_plot %>% filter(hit_dir != "none") %>%
  distinct(gene, hit_dir, best_lfc)
cat(sprintf("\nHits on diagram: %d / %d (%.0f%%)\n",
            nrow(hits_on_diagram), nrow(hits),
            100*nrow(hits_on_diagram)/nrow(hits)))
cat(paste(sprintf("  %-12s %s lfc=%.2f",
                  hits_on_diagram$gene, hits_on_diagram$hit_dir,
                  hits_on_diagram$best_lfc), collapse="\n"), "\n")

# Hits NOT on diagram (overflow strip)
hits_overflow <- hits %>%
  filter(!gene %in% nodes_plot$gene)
cat(sprintf("Overflow strip: %d hits\n", nrow(hits_overflow)))

# ── 8. Overflow strip for remaining hits ─────────────────────
# Group by functional category, sorted by |LFC|
overflow_group <- function(gene) {
  case_when(
    grepl("^RP[SL]|^FAU$|MRPL|MRPS", gene)               ~ "Ribosome",
    grepl("^EIF|^EEF|TRMT|GSPT",     gene)               ~ "Translation",
    grepl("^POLR|^TAF|SUPT|NELF|^GTF|^MED|^BRD",gene)   ~ "Transcription machinery",
    grepl("^NOP|^DDX|^DHX|^SNRP|^SF3|^LSM|^U2AF",gene)  ~ "RNA processing",
    gene %in% c("CCND2","CDK4","CDK6","RB1","E2F1")       ~ "Cell cycle",
    gene %in% c("DAD1","DAPK1","BIRC2","BIRC3","XIAP")    ~ "Apoptosis",
    gene %in% c("GNA13","GNAI2","GNAQ","GNB1")             ~ "G-protein",
    gene %in% c("PCBP2","AGO1","AGO2","DICER1")            ~ "RNA regulators",
    gene %in% c("FIBP","GLRX","NDUFB10","C2orf68",
                "RPRD1B","BIRC6","GNA13")                  ~ "Metabolic/Other",
    TRUE ~ "Signalling"
  )
}

if (nrow(hits_overflow) > 0) {
  overflow <- hits_overflow %>%
    mutate(grp = overflow_group(gene)) %>%
    arrange(grp, hit_dir, desc(abs(best_lfc)))

  # Layout: tile each gene in a row, grouped by category
  grp_levels <- unique(overflow$grp)
  overflow <- overflow %>%
    group_by(grp) %>%
    mutate(x_in_grp = row_number()) %>%
    ungroup()

  # Assign x start per group (pack groups sequentially)
  grp_widths <- overflow %>% group_by(grp) %>%
    summarise(n=n(), .groups="drop") %>%
    mutate(x_start = cumsum(c(0, head(n+1,-1))) * 0.88 + 0.5)
  overflow <- overflow %>%
    left_join(grp_widths %>% select(grp, x_start), by="grp") %>%
    mutate(x_tile = (x_in_grp-1)*0.88 + x_start,
           fill_col   = if_else(hit_dir=="depleted","#2166AC","#D6604D"),
           border_col = if_else(hit_dir=="depleted","#053061","#67000D"))
}

# ── 9. Edge coordinates ───────────────────────────────────────
# Build position lookup from all nodes
pos_lookup <- nodes_plot %>% select(gene, x, y)

edges_coord <- backbone_edges %>%
  left_join(pos_lookup, by=c("from"="gene")) %>% rename(x1=x, y1=y) %>%
  left_join(pos_lookup, by=c("to"="gene"))   %>% rename(x2=x, y2=y) %>%
  filter(!is.na(x1), !is.na(x2), !(x1==x2 & y1==y2))

edge_style <- tribble(
  ~type,             ~ecol,      ~elty,
  "activation",      "#1F618D",  "solid",
  "phosphorylation", "#935116",  "solid",
  "inhibition",      "#922B21",  "solid",
  "indirect_act",    "#1F618D",  "dashed",
  "binding",         "#5D6D7E",  "dashed",
  "translocation",   "#6C3483",  "dotted"
)
edges_coord <- edges_coord %>%
  left_join(edge_style, by="type") %>%
  mutate(ecol = replace_na(ecol,"#888888"),
         elty = replace_na(elty,"solid"))

# ── 10. Drawing helpers ───────────────────────────────────────
draw_seg <- function(df, ecol_v, elty_v, arr_v, lwd=0.38, al=0.80) {
  sub <- df %>% filter(ecol==ecol_v, elty==elty_v, !is_cross)
  if (nrow(sub)==0) return(NULL)
  arr <- switch(arr_v,
    "open"  = arrow(type="open", length=unit(3.8,"pt"), angle=22),
    "blunt" = arrow(type="open", length=unit(3.0,"pt"), angle=90, ends="last"),
    NULL)
  geom_segment(data=sub, aes(x=x1,y=y1,xend=x2,yend=y2),
    colour=ecol_v, linetype=elty_v, linewidth=lwd,
    alpha=al, arrow=arr, lineend="round", show.legend=FALSE)
}

draw_arc <- function(df, ecol_v, elty_v, arr_v, lwd=0.32) {
  sub <- df %>% filter(ecol==ecol_v, elty==elty_v, is_cross)
  if (nrow(sub)==0) return(NULL)
  arr <- switch(arr_v,
    "open"  = arrow(type="open", length=unit(3.5,"pt"), angle=22),
    "blunt" = arrow(type="open", length=unit(3.0,"pt"), angle=90, ends="last"),
    NULL)
  geom_curve(data=sub, aes(x=x1,y=y1,xend=x2,yend=y2),
    colour=ecol_v, linetype=elty_v, linewidth=lwd,
    curvature=0.22, alpha=0.50, arrow=arr,
    lineend="round", show.legend=FALSE)
}

# ── 11. Module header positions ───────────────────────────────
module_info <- tribble(
  ~name,            ~xmid,  ~col,
  "TCR Proximal",    1.7,   "#2471A3",
  "MAPK/NFAT",       5.2,   "#1E8449",
  "PI3K/AKT",        8.5,   "#7D3C98",
  "NF-κB",          12.2,   "#D35400",
  "JAK-STAT",       15.7,   "#B03A2E",
  "Epi/Apoptosis",  19.2,   "#148F77"
)

# ── 12. Compartment boxes ─────────────────────────────────────
XMIN <- -0.5; XMAX <- 21.5
comp_boxes <- tribble(
  ~clabel,           ~ymin, ~ymax,  ~cfill,
  "Extracellular",   10.2,  11.1,  "#D6EAF8",
  "Plasma Membrane",  9.0,  10.15, "#D5F5E3",
  "Cytoplasm",        1.1,   8.95, "#FDFEFE",
  "Nucleus",         -0.75,  1.05, "#FEF9E7",
  "Mitochondria",     4.2,   6.3,  "#FDEDEC"
) %>% mutate(xmin=XMIN, xmax=XMAX)

mod_dividers <- tibble(
  x    = c(3.4, 6.9, 10.4, 14.1, 17.4),
  ymin = -0.7, ymax = 10.1)

# ── 13. Build ggplot ─────────────────────────────────────────
p_main <- ggplot() +

  # Compartment fills
  geom_rect(data=comp_boxes,
    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=clabel),
    colour="#CCCCCC", linewidth=0.20, alpha=0.40) +
  scale_fill_manual(name="Compartment",
    values=setNames(comp_boxes$cfill, comp_boxes$clabel)) +
  new_scale_fill() +

  # Module dividers
  geom_segment(data=mod_dividers,
    aes(x=x, xend=x, y=ymin, yend=ymax),
    colour="#C0C0C0", linewidth=0.25, linetype="dashed") +

  # Module headers
  annotate("rect",
    xmin=module_info$xmid-1.4, xmax=module_info$xmid+1.4,
    ymin=10.55, ymax=10.95,
    fill=alpha(module_info$col, 0.12), colour=module_info$col,
    linewidth=0.45) +
  annotate("text", x=module_info$xmid, y=10.75,
    label=module_info$name, colour=module_info$col,
    size=2.5, fontface="bold", hjust=0.5) +

  # ── Edges (segments first, then arcs on top) ──
  draw_seg(edges_coord,"#1F618D","solid", "open")  +
  draw_seg(edges_coord,"#935116","solid", "open")  +
  draw_seg(edges_coord,"#922B21","solid", "blunt") +
  draw_seg(edges_coord,"#1F618D","dashed","open")  +
  draw_seg(edges_coord,"#5D6D7E","dashed","none")  +
  draw_seg(edges_coord,"#6C3483","dotted","open")  +
  draw_arc(edges_coord,"#1F618D","solid", "open")  +
  draw_arc(edges_coord,"#935116","solid", "open")  +
  draw_arc(edges_coord,"#922B21","solid", "blunt") +
  draw_arc(edges_coord,"#1F618D","dashed","open")  +

  # ── Connector (non-hit) nodes — drawn first (beneath hits) ──
  {
    cn <- nodes_plot %>% filter(is_connector | hit_dir=="none")
    rn <- cn %>% filter(node_type %in% c("receptor","complex","output","ligand"))
    kn <- cn %>% filter(node_type=="kinase")
    an <- cn %>% filter(node_type=="adaptor")
    tn <- cn %>% filter(node_type=="tf")
    list(
      if(nrow(rn)>0) geom_tile(data=rn,
        aes(x=x,y=y,width=node_r*2.2,height=node_r,fill=fill_col,colour=border_col),
        linewidth=0.45, show.legend=FALSE) else NULL,
      if(nrow(kn)>0) geom_regon(data=kn,
        aes(x0=x,y0=y,r=node_r,angle=pi/6,sides=6,fill=fill_col,colour=border_col),
        linewidth=0.45, show.legend=FALSE) else NULL,
      if(nrow(an)>0) geom_ellipse(data=an,
        aes(x0=x,y0=y,a=node_r*1.5,b=node_r*0.7,angle=0,fill=fill_col,colour=border_col),
        linewidth=0.45, show.legend=FALSE) else NULL,
      if(nrow(tn)>0) geom_regon(data=tn,
        aes(x0=x,y0=y,r=node_r*1.2,angle=pi/4,sides=4,fill=fill_col,colour=border_col),
        linewidth=0.45, show.legend=FALSE) else NULL
    )
  } +
  scale_fill_identity() + scale_colour_identity() +
  new_scale_fill() + new_scale_colour() +

  # ── Hit nodes — drawn ON TOP (larger, coloured) ──
  {
    hn <- nodes_plot %>% filter(hit_dir != "none")
    rn <- hn %>% filter(node_type %in% c("receptor","complex","output","ligand"))
    kn <- hn %>% filter(node_type=="kinase")
    an <- hn %>% filter(node_type=="adaptor")
    tn <- hn %>% filter(node_type=="tf")
    list(
      if(nrow(rn)>0) geom_tile(data=rn,
        aes(x=x,y=y,width=node_r*2.2,height=node_r,fill=fill_col,colour=border_col),
        linewidth=0.80, show.legend=FALSE) else NULL,
      if(nrow(kn)>0) geom_regon(data=kn,
        aes(x0=x,y0=y,r=node_r,angle=pi/6,sides=6,fill=fill_col,colour=border_col),
        linewidth=0.80, show.legend=FALSE) else NULL,
      if(nrow(an)>0) geom_ellipse(data=an,
        aes(x0=x,y0=y,a=node_r*1.5,b=node_r*0.7,angle=0,fill=fill_col,colour=border_col),
        linewidth=0.80, show.legend=FALSE) else NULL,
      if(nrow(tn)>0) geom_regon(data=tn,
        aes(x0=x,y0=y,r=node_r*1.2,angle=pi/4,sides=4,fill=fill_col,colour=border_col),
        linewidth=0.80, show.legend=FALSE) else NULL
    )
  } +
  scale_fill_identity() + scale_colour_identity() +
  new_scale_fill() + new_scale_colour() +

  # ── Labels: hits (bold, repel) ──
  geom_label_repel(
    data = nodes_plot %>% filter(hit_dir != "none"),
    aes(x=x, y=y, label=label, fontface="bold"),
    size=2.0, label.padding=unit(0.06,"lines"), label.size=0.12,
    fill=alpha("white",0.85), colour="#111111",
    box.padding=0.10, point.padding=0.10,
    segment.linewidth=0.20, segment.colour="#AAAAAA",
    max.overlaps=80, seed=42, show.legend=FALSE) +

  # ── Labels: connectors / non-hits (smaller, plain) ──
  geom_text(
    data = nodes_plot %>% filter(hit_dir == "none"),
    aes(x=x, y=y, label=label),
    size=1.55, colour="#777777", fontface="plain",
    vjust=-1.0, lineheight=0.85) +

  # ── Overflow strip: remaining hits ──
  { if (exists("overflow") && nrow(hits_overflow)>0) {
    STRIP_Y <- -1.8
    GRP_Y   <- -2.25
    grp_pos <- overflow %>%
      group_by(grp) %>%
      summarise(gx = mean(x_tile), .groups="drop")
    list(
      # Group labels
      annotate("text", x=grp_pos$gx, y=GRP_Y,
        label=grp_pos$grp, size=1.75, colour="#555555",
        fontface="bold.italic", hjust=0.5),
      # Tiles
      geom_tile(data=overflow,
        aes(x=x_tile, y=STRIP_Y, fill=fill_col, colour=border_col),
        width=0.82, height=0.38, linewidth=0.45, show.legend=FALSE),
      scale_fill_identity(),
      scale_colour_identity(),
      new_scale_fill(), new_scale_colour(),
      # Gene names inside tiles
      geom_text(data=overflow,
        aes(x=x_tile, y=STRIP_Y, label=gene),
        size=1.35, colour="white", fontface="bold", hjust=0.5),
      # Strip header
      annotate("text", x=-0.3, y=STRIP_Y+0.05,
        label="Additional hits (not in backbone pathway):",
        size=1.9, colour="#333333", fontface="bold", hjust=0)
    )
  } else NULL } +

  # ── Edge legend ──
  annotate("segment",x=0.0,xend=0.6,y=-1.15,yend=-1.15,colour="#1F618D",linewidth=0.6,
    arrow=arrow(type="open",length=unit(3.5,"pt"),angle=22)) +
  annotate("text",x=0.7,y=-1.15,label="Activation",    size=1.85,hjust=0) +
  annotate("segment",x=0.0,xend=0.6,y=-1.40,yend=-1.40,colour="#935116",linewidth=0.6,
    arrow=arrow(type="open",length=unit(3.5,"pt"),angle=22)) +
  annotate("text",x=0.7,y=-1.40,label="Phosphorylation",size=1.85,hjust=0) +
  annotate("segment",x=3.0,xend=3.6,y=-1.15,yend=-1.15,colour="#922B21",linewidth=0.6,
    arrow=arrow(type="open",length=unit(3.0,"pt"),angle=90,ends="last")) +
  annotate("text",x=3.7,y=-1.15,label="Inhibition",    size=1.85,hjust=0) +
  annotate("segment",x=3.0,xend=3.6,y=-1.40,yend=-1.40,colour="#1F618D",linewidth=0.6,linetype="dashed") +
  annotate("text",x=3.7,y=-1.40,label="Indirect",      size=1.85,hjust=0) +
  annotate("segment",x=6.2,xend=6.8,y=-1.15,yend=-1.15,colour="#5D6D7E",linewidth=0.6,linetype="dashed") +
  annotate("text",x=6.9,y=-1.15,label="Binding",       size=1.85,hjust=0) +
  annotate("segment",x=6.2,xend=6.8,y=-1.40,yend=-1.40,colour="#6C3483",linewidth=0.6,linetype="dotted",
    arrow=arrow(type="open",length=unit(3.5,"pt"),angle=22)) +
  annotate("text",x=6.9,y=-1.40,label="Translocation", size=1.85,hjust=0) +
  # CRISPRi legend
  annotate("point",x=9.8, y=-1.15,colour="#053061",fill="#2166AC",size=3.0,shape=22,stroke=0.9) +
  annotate("text", x=10.1,y=-1.15,label=sprintf("Depleted/GoF (n=%d)",nrow(hits_dep)),size=1.85,hjust=0) +
  annotate("point",x=9.8, y=-1.40,colour="#67000D",fill="#D6604D",size=3.0,shape=22,stroke=0.9) +
  annotate("text", x=10.1,y=-1.40,label=sprintf("Enriched/LoF (n=%d)",nrow(hits_enr)),size=1.85,hjust=0) +
  annotate("point",x=14.0,y=-1.15,colour="#909497",fill="#D5D8DC",size=2.2,shape=22,stroke=0.6) +
  annotate("text", x=14.3,y=-1.15,label="Backbone (not hit)",size=1.85,hjust=0) +
  annotate("text", x=14.0,y=-1.40,size=1.7,hjust=0,colour="#666666",fontface="italic",
    label="⬡ Kinase  ◇ TF  ○ Adaptor  ▭ Receptor/Complex  | size ∝ −log₁₀FDR") +

  # Axes / theme
  coord_cartesian(xlim=c(XMIN, XMAX),
                  ylim=c(if(nrow(hits_overflow)>0) -2.6 else -1.7, 11.2),
                  clip="off") +
  labs(
    title    = "Integrated T Cell Signalling — CRISPRi Screen (Carnevale et al.)",
    subtitle = sprintf(
      "%d GoF hits (blue) | %d LoF hits (red) | %d on pathway | %d in strip\nTop ORA: %s",
      nrow(hits_dep), nrow(hits_enr),
      nrow(hits_on_diagram),
      nrow(hits_overflow),
      str_trunc(ora_subtitle, 90)),
    caption  = "CRISPRi: depleted=GoF (knockdown impairs T cells); enriched=LoF (knockdown enhances T cells)"
  ) +
  theme_void(base_size=10) +
  theme(
    plot.title    = element_text(size=11, face="bold",  hjust=0.5, margin=margin(b=2)),
    plot.subtitle = element_text(size=6.8, colour="#444444", hjust=0.5, margin=margin(b=1)),
    plot.caption  = element_text(size=6.2, colour="#888888", hjust=0.5, margin=margin(t=2)),
    plot.margin   = margin(8,6,14,6),
    legend.position = "right",
    legend.title  = element_text(size=8, face="bold"),
    legend.text   = element_text(size=7)
  )

# ── 14. Export ────────────────────────────────────────────────
cat("\nSaving output...\n")
ggsave(paste0(OUT_PREFIX,".pdf"), p_main,
       width=28, height=16, units="in", device=cairo_pdf)
ggsave(paste0(OUT_PREFIX,".png"), p_main,
       width=28, height=16, units="in", dpi=170)
cat(sprintf("  %s.pdf / .png\n", OUT_PREFIX))

# ORA dotplot
if (nrow(ora_df) >= 3) {
  ora_top <- ora_df %>% slice_head(n=20) %>%
    mutate(Description = str_trunc(Description, 60),
           Description = fct_reorder(Description, -p.adjust))
  p_ora <- ggplot(ora_top,
    aes(x=GeneRatio, y=Description, size=Count, colour=p.adjust)) +
    geom_point() +
    scale_colour_gradient(low="#C0392B", high="#AEB6BF", name="FDR") +
    scale_size_continuous(range=c(2,8), name="Gene count") +
    labs(title="Reactome ORA — all CRISPRi hits", x="Gene Ratio", y=NULL) +
    theme_bw(base_size=10) +
    theme(axis.text.y = element_text(size=8))
  ggsave(paste0(OUT_PREFIX,"_ORA.pdf"), p_ora,
         width=10, height=7, units="in", device=cairo_pdf)
  cat(sprintf("  %s_ORA.pdf\n", OUT_PREFIX))
}

# Full hit summary table
hit_summary <- bind_rows(
  hits_on_diagram %>%
    left_join(nodes_plot %>% distinct(gene,module,node_type), by="gene") %>%
    mutate(location="pathway"),
  if(nrow(hits_overflow)>0)
    hits_overflow %>% mutate(module="strip", node_type=NA, location="strip")
  else NULL
) %>% arrange(location, hit_dir, best_fdr)
write.csv(hit_summary, paste0(OUT_PREFIX,"_hits_table.csv"), row.names=FALSE)
cat(sprintf("  %s_hits_table.csv  (%d rows)\n", OUT_PREFIX, nrow(hit_summary)))

cat("\n── Done ──\n")
cat(sprintf("On-pathway hits : %d / %d (%.0f%%)\n",
            nrow(hits_on_diagram), nrow(hits),
            100*nrow(hits_on_diagram)/nrow(hits)))
cat(sprintf("Strip hits      : %d\n", nrow(hits_overflow)))
cat(sprintf("Total shown     : %d / %d (%.0f%%)\n",
            nrow(hits_on_diagram)+nrow(hits_overflow), nrow(hits),
            100*(nrow(hits_on_diagram)+nrow(hits_overflow))/nrow(hits)))
