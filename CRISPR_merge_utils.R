# ============================================================
# CRISPR_merge_utils.R
#
# Utilities for importing and merging multiple MAGeCK
# CRISPRi/CRISPRko/CRISPRa gene-level screens, then
# producing a merged node table compatible with
# CRISPR_to_Tcell_pathway.R.
#
# Key functions:
#   load_mageck()    — read one MAGeCK gene_summary.txt
#   merge_assays()   — union/intersect multiple assay objects
#   assay_summary()  — concordance report
#
# Usage (add after "── 1. Settings" in main script):
#   source("C:/Users/difen/Rcode/CRISPR_merge_utils.R")
#   a1 <- load_mageck(INPUT,  name = "Carnevale_CRISPRi", convention = "CRISPRi")
#   a2 <- load_mageck("path/to/screen2.txt", name = "Screen2_CRISPRko",
#                     convention = "CRISPRko")
#   mg <- merge_assays(list(a1, a2), method = "union")
#   dat       <- mg$dat          # drop-in replacement for dat in main script
#   hits      <- mg$hits
#   hits_dep  <- mg$hits_dep
#   hits_enr  <- mg$hits_enr
# ============================================================

# ── Dependencies (assume already loaded by main script) ─────
# dplyr, tibble, tidyr, scales, stringr, purrr

# ── Convention guide ─────────────────────────────────────────
# "CRISPRi"  : knockdown. depleted gene → LOSS of gene product
#              neg.lfc < 0 (depleted) means gene was NEEDED → GoF hit (blue)
#              pos.lfc > 0 (enriched) means gene suppressed T cell → LoF hit (red)
#              (Carnevale convention)
#
# "CRISPRko" : knockout. same direction mapping as CRISPRi for most screens
#              but some labs flip terminology. Use "CRISPRko" here if the
#              screen is knockout-based with the same MAGeCK sign convention.
#              If your lab calls depleted=LoF, set convention="CRISPRko_flip".
#
# "CRISPRa"  : activation. depleted = gene overexpression IMPAIRED cells → LoF
#              signs are opposite: neg.lfc < 0 → enriched hit direction.
#              Set convention="CRISPRa" to auto-flip.

# ─────────────────────────────────────────────────────────────
# load_mageck()
# ─────────────────────────────────────────────────────────────
#' Read a MAGeCK gene_summary.txt and return a standardised assay object.
#'
#' @param path       Path to MAGeCK gene_summary file.
#' @param name       Short label for this assay (used in merged columns).
#' @param fdr_cut    FDR threshold (default 0.10).
#' @param lfc_cut    |LFC| threshold (default 0.20).
#' @param convention One of "CRISPRi", "CRISPRko", "CRISPRko_flip", "CRISPRa".
#' @param cell_type  Optional free-text cell type label for plots.
#'
#' @return A named list:
#'   $name, $convention, $cell_type,
#'   $dat      (all genes, standardised columns),
#'   $hits     (significant genes only),
#'   $hits_dep (depleted / GoF hits),
#'   $hits_enr (enriched / LoF hits)
load_mageck <- function(path,
                        name          = basename(path),
                        fdr_cut       = 0.10,
                        lfc_cut       = 0.20,
                        convention    = c("CRISPRi", "CRISPRko",
                                          "CRISPRko_flip", "CRISPRa"),
                        cell_type     = NA_character_,
                        toupper_genes = FALSE) {

  convention <- match.arg(convention)
  stopifnot(file.exists(path))

  raw <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("[%s] Loaded: %d genes from %s\n", name, nrow(raw), path))

  # ── Column auto-detection ────────────────────────────────────
  # MAGeCK ≥0.5.9 uses 'id'; older may use 'Gene' or first column
  gene_col <- intersect(c("id", "Gene", "gene"), colnames(raw))[1]
  if (is.na(gene_col)) gene_col <- colnames(raw)[1]

  # neg/pos column detection (some files prefix with group names)
  neg_lfc_col <- grep("neg.*lfc|neg\\.lfc",  colnames(raw), value = TRUE, ignore.case = TRUE)[1]
  neg_fdr_col <- grep("neg.*fdr|neg\\.fdr",  colnames(raw), value = TRUE, ignore.case = TRUE)[1]
  pos_lfc_col <- grep("pos.*lfc|pos\\.lfc",  colnames(raw), value = TRUE, ignore.case = TRUE)[1]
  pos_fdr_col <- grep("pos.*fdr|pos\\.fdr",  colnames(raw), value = TRUE, ignore.case = TRUE)[1]

  missing_cols <- is.na(c(neg_lfc_col, neg_fdr_col, pos_lfc_col, pos_fdr_col))
  if (any(missing_cols)) {
    stop(sprintf(
      "[%s] Cannot detect MAGeCK columns. Found: %s",
      name, paste(colnames(raw), collapse = ", ")))
  }

  # ── Standardise ─────────────────────────────────────────────
  dat <- raw %>%
    transmute(
      gene    = .data[[gene_col]],
      neg_lfc = .data[[neg_lfc_col]],
      neg_fdr = .data[[neg_fdr_col]],
      pos_lfc = .data[[pos_lfc_col]],
      pos_fdr = .data[[pos_fdr_col]]
    ) %>%
    filter(!is.na(gene), gene != "")

  # ── Normalise gene symbols to uppercase (for mouse→human mapping) ───
  if (toupper_genes) {
    dat <- dat %>% mutate(gene = toupper(gene))
    cat(sprintf("[%s] Gene names converted to uppercase (toupper_genes=TRUE)\n", name))
  }

  # ── Apply convention (flip if needed) ───────────────────────
  # CRISPRa: activation → flip which direction = "depleted" hit
  # CRISPRko_flip: lab uses reversed terminology
  if (convention %in% c("CRISPRa", "CRISPRko_flip")) {
    dat <- dat %>%
      rename(neg_lfc_raw = neg_lfc, neg_fdr_raw = neg_fdr,
             pos_lfc_raw = pos_lfc, pos_fdr_raw = pos_fdr) %>%
      mutate(neg_lfc = pos_lfc_raw, neg_fdr = pos_fdr_raw,
             pos_lfc = neg_lfc_raw, pos_fdr = neg_fdr_raw) %>%
      select(gene, neg_lfc, neg_fdr, pos_lfc, pos_fdr)
    cat(sprintf("[%s] Columns flipped for convention='%s'\n", name, convention))
  }

  # ── Classify hits ────────────────────────────────────────────
  dat <- dat %>%
    mutate(
      hit_dir  = case_when(
        neg_fdr < fdr_cut & neg_lfc < -lfc_cut ~ "depleted",
        pos_fdr < fdr_cut & pos_lfc >  lfc_cut ~ "enriched",
        TRUE ~ "none"),
      best_lfc = case_when(
        hit_dir == "depleted" ~ neg_lfc,
        hit_dir == "enriched" ~ pos_lfc,
        TRUE ~ 0),
      best_fdr = pmin(neg_fdr, pos_fdr),
      nlp      = -log10(pmax(pmin(neg_fdr, pos_fdr), 1e-6)),
      assay    = name
    )

  hits     <- dat %>% filter(hit_dir != "none")
  hits_dep <- hits %>% filter(hit_dir == "depleted") %>% arrange(neg_fdr, neg_lfc)
  hits_enr <- hits %>% filter(hit_dir == "enriched") %>% arrange(pos_fdr, desc(pos_lfc))

  cat(sprintf("[%s] Hits: %d depleted (GoF/blue), %d enriched (LoF/red)\n",
              name, nrow(hits_dep), nrow(hits_enr)))

  list(
    name       = name,
    convention = convention,
    cell_type  = cell_type,
    fdr_cut    = fdr_cut,
    lfc_cut    = lfc_cut,
    dat        = dat,
    hits       = hits,
    hits_dep   = hits_dep,
    hits_enr   = hits_enr
  )
}


# ─────────────────────────────────────────────────────────────
# merge_assays()
# ─────────────────────────────────────────────────────────────
#' Merge two or more assay objects into a single gene-level table.
#'
#' For each gene the function computes:
#'   - Per-assay LFC, FDR, direction
#'   - Concordance across assays (all agree / mixed / single)
#'   - Meta-significance: max(−log10 FDR) across assays
#'   - Merged hit_dir for node colouring (see `consensus_rule`)
#'   - Visual encoding columns: fill_col, border_col, node_style
#'
#' @param assay_list  Named list of load_mageck() outputs.
#' @param method      "union" = gene is a hit in ANY assay (default);
#'                    "intersect" = gene must be a hit in ALL assays.
#' @param consensus_rule
#'   "majority" = direction with most assay votes wins (default);
#'   "strict"   = only concordant genes carry a direction (discordant → "conflict");
#'   "first"    = direction from the first assay in assay_list.
#' @param min_assays  Minimum number of assays in which a gene must be a hit
#'                    (only applies when method = "union"; default = 1).
#'
#' @return A named list:
#'   $dat       : gene × (assay columns) merged tibble – drop-in for main `dat`
#'   $hits      : significant genes
#'   $hits_dep  : depleted consensus hits
#'   $hits_enr  : enriched consensus hits
#'   $concordance: gene-level concordance table
#'   $assay_names: character vector of assay names in order
merge_assays <- function(assay_list,
                         method         = c("union", "intersect"),
                         consensus_rule = c("majority", "strict", "first"),
                         min_assays     = 1L) {

  method         <- match.arg(method)
  consensus_rule <- match.arg(consensus_rule)
  stopifnot(length(assay_list) >= 2)

  assay_names <- vapply(assay_list, `[[`, character(1), "name")
  n_assays    <- length(assay_list)
  cat(sprintf("\nMerging %d assays: %s\n", n_assays, paste(assay_names, collapse = " + ")))

  # ── 1. Wide join: one row per gene, one block of cols per assay ──
  wide <- assay_list[[1]]$dat %>%
    select(gene, neg_lfc, neg_fdr, pos_lfc, pos_fdr, hit_dir, best_lfc, best_fdr, nlp) %>%
    rename_with(~ paste0(., "__", assay_names[1]), -gene)

  for (i in seq_along(assay_list)[-1]) {
    aname <- assay_names[i]
    chunk <- assay_list[[i]]$dat %>%
      select(gene, neg_lfc, neg_fdr, pos_lfc, pos_fdr, hit_dir, best_lfc, best_fdr, nlp) %>%
      rename_with(~ paste0(., "__", aname), -gene)
    wide <- full_join(wide, chunk, by = "gene")
  }

  # ── 2. Concordance columns ───────────────────────────────────
  dir_cols <- paste0("hit_dir__", assay_names)   # one per assay

  # Helper: Fisher's combined p-value for a vector of independent p-values
  fisher_combine <- function(p_vec) {
    p_vec <- p_vec[!is.na(p_vec) & p_vec > 0 & p_vec <= 1]
    if (length(p_vec) == 0) return(NA_real_)
    if (length(p_vec) == 1) return(p_vec[1])
    chi2 <- -2 * sum(log(p_vec))
    pchisq(chi2, df = 2L * length(p_vec), lower.tail = FALSE)
  }

  concordance <- wide %>%
    mutate(across(all_of(dir_cols), ~ replace_na(.x, "none"))) %>%
    rowwise() %>%
    mutate(
      # Which assays call this gene a hit
      n_hit       = sum(c_across(all_of(dir_cols)) != "none"),
      n_depleted  = sum(c_across(all_of(dir_cols)) == "depleted"),
      n_enriched  = sum(c_across(all_of(dir_cols)) == "enriched"),
      # Concordance label
      concordance = case_when(
        n_hit == 0          ~ "not_hit",
        n_hit == 1          ~ "single_assay",
        n_depleted == n_hit ~ "concordant_depleted",
        n_enriched == n_hit ~ "concordant_enriched",
        TRUE                ~ "discordant"
      ),
      # max(nlp) — used as fallback and for per-assay node sizing
      meta_nlp = {
        nlp_vals <- c_across(paste0("nlp__", assay_names))
        dir_vals <- c_across(all_of(dir_cols))
        nlp_vals[dir_vals == "none"] <- 0
        max(nlp_vals, na.rm = TRUE)
      },
      # Fisher's combined p across assays where gene is a hit
      # Rewards genes replicated in multiple screens; node size reflects this
      fisher_nlp = {
        fdr_vals <- c_across(paste0("best_fdr__", assay_names))
        dir_vals <- c_across(all_of(dir_cols))
        p_hit    <- fdr_vals[dir_vals != "none"]
        fp       <- fisher_combine(p_hit)
        if (is.na(fp)) meta_nlp else -log10(max(fp, 1e-10))
      },
      # Meta-LFC: from the highest-significance assay
      meta_lfc = {
        nlp_vals <- c_across(paste0("nlp__", assay_names))
        lfc_vals <- c_across(paste0("best_lfc__", assay_names))
        dir_vals <- c_across(all_of(dir_cols))
        nlp_vals[dir_vals == "none"] <- -Inf
        lfc_vals[which.max(nlp_vals)]
      }
    ) %>%
    ungroup()

  # ── 3. Consensus direction ───────────────────────────────────
  # Pre-compute best_fdr outside mutate() to avoid .data[[cn]] scoping
  # issues when cn is a lapply argument inside dplyr::mutate().
  best_fdr_vec <- {
    all_fdr_cols <- c(paste0("neg_fdr__", assay_names),
                      paste0("pos_fdr__", assay_names))
    Reduce(pmin,
           lapply(all_fdr_cols,
                  function(cn) replace_na(concordance[[cn]], 1)))
  }

  concordance <- concordance %>%
    mutate(
      hit_dir = case_when(
        concordance == "not_hit"  ~ "none",
        consensus_rule == "first" ~ {
          first_dir <- .data[[paste0("hit_dir__", assay_names[1])]]
          if_else(first_dir != "none", first_dir, "none")
        },
        consensus_rule == "strict" ~ case_when(
          concordance == "concordant_depleted" ~ "depleted",
          concordance == "concordant_enriched" ~ "enriched",
          TRUE ~ "conflict"
        ),
        # majority (default): direction with most votes wins
        n_depleted > n_enriched ~ "depleted",
        n_enriched > n_depleted ~ "enriched",
        TRUE                    ~ "conflict"   # exact tie
      ),
      best_fdr = best_fdr_vec,          # pre-computed above
      nlp      = fisher_nlp,
      best_lfc = meta_lfc,
      assay    = "merged"
    )

  # ── 4. Apply hit selection method ───────────────────────────
  if (method == "union") {
    hit_genes <- concordance %>%
      filter(n_hit >= min_assays) %>%
      pull(gene)
  } else {  # intersect
    hit_genes <- concordance %>%
      filter(n_hit == n_assays) %>%
      pull(gene)
  }

  dat_merged <- concordance %>%
    select(gene, hit_dir, best_lfc, best_fdr, nlp, assay,
           concordance, n_hit, n_depleted, n_enriched,
           all_of(paste0("hit_dir__",  assay_names)),
           all_of(paste0("best_lfc__", assay_names)),
           all_of(paste0("nlp__",      assay_names)))

  hits_merged <- dat_merged %>%
    filter(gene %in% hit_genes, hit_dir %in% c("depleted","enriched","conflict"))

  hits_dep <- hits_merged %>% filter(hit_dir == "depleted") %>%
    arrange(desc(n_depleted), desc(nlp))
  hits_enr <- hits_merged %>% filter(hit_dir == "enriched") %>%
    arrange(desc(n_enriched), desc(nlp))

  cat(sprintf(
    "Merged hits [%s, min_assays=%d]: %d depleted · %d enriched · %d conflict\n",
    method, min_assays, nrow(hits_dep), nrow(hits_enr),
    sum(hits_merged$hit_dir == "conflict")))
  cat(sprintf("Concordance: %d concordant · %d single-assay · %d discordant\n",
    sum(grepl("^concordant", hits_merged$concordance)),
    sum(hits_merged$concordance == "single_assay"),
    sum(hits_merged$concordance == "discordant")))

  list(
    dat            = dat_merged,
    hits           = hits_merged,
    hits_dep       = hits_dep,
    hits_enr       = hits_enr,
    concordance    = concordance,
    assay_names    = assay_names,
    n_assays       = n_assays,
    method         = method,
    consensus_rule = consensus_rule
  )
}


# ─────────────────────────────────────────────────────────────
# ASSAY_PALETTE  — per-assay colour families (up to 4 screens)
# ─────────────────────────────────────────────────────────────
# GoF = cool hues (blue → green → purple → teal)
# LoF = warm hues (red  → orange → pink  → amber)
# Concordant = darkest shade of GoF/LoF (highest confidence)
#
# Row order matches assay_names order in ASSAY_CONFIG.
ASSAY_PALETTE <- tibble::tribble(
  ~gof_fill,   ~gof_border,  ~lof_fill,   ~lof_border,
  "#4292C6",   "#2171B5",    "#FB6A4A",   "#CB181D",   # assay 1 – blue / red
  "#238B45",   "#005824",    "#E08214",   "#8C4A00",   # assay 2 – green / orange
  "#6A3D9A",   "#3F007D",    "#F768A1",   "#AE017E",   # assay 3 – purple / pink
  "#1D9E9E",   "#006D6D",    "#D95F02",   "#993D00"    # assay 4 – teal / burnt-orange
)

# ─────────────────────────────────────────────────────────────
# node_colours_merged()
# ─────────────────────────────────────────────────────────────
#' Map merged concordance → fill/border colours, with per-assay
#' distinction for single-screen hits.
#'
#' Encoding
#'   concordant GoF  (all screens agree) → deep blue   #08519C
#'   concordant LoF  (all screens agree) → deep red    #CB181D
#'   assay-N-only GoF                    → assay colour from ASSAY_PALETTE row N
#'   assay-N-only LoF                    → assay colour from ASSAY_PALETTE row N
#'   discordant / conflict               → gold        #D4AC0D
#'
#' @param hits_merged  $hits from merge_assays().
#' @param assay_names  Character vector of assay names (mg$assay_names).
node_colours_merged <- function(hits_merged, assay_names) {

  dir_cols <- paste0("hit_dir__", assay_names)
  n        <- length(assay_names)

  # Palette rows for each assay (recycle if >4 assays)
  pal <- ASSAY_PALETTE[((seq_len(n) - 1L) %% nrow(ASSAY_PALETTE)) + 1L, ]
  pal$assay <- assay_names

  # Identify which single assay is the source for single-assay hits
  source_tbl <- hits_merged %>%
    dplyr::select(gene, concordance, hit_dir,
                  dplyr::all_of(dir_cols)) %>%
    rowwise() %>%
    dplyr::mutate(
      source_assay = {
        dirs    <- c_across(dplyr::all_of(dir_cols))
        hit_idx <- which(dirs != "none")
        if (length(hit_idx) == 1L) assay_names[hit_idx] else NA_character_
      }
    ) %>%
    ungroup() %>%
    dplyr::left_join(pal, by = c("source_assay" = "assay"))

  source_tbl %>%
    dplyr::transmute(
      gene,
      fill_col = dplyr::case_when(
        concordance == "concordant_depleted"                   ~ "#08519C",
        concordance == "concordant_enriched"                   ~ "#CB181D",
        concordance == "single_assay" & hit_dir == "depleted"  ~ gof_fill,
        concordance == "single_assay" & hit_dir == "enriched"  ~ lof_fill,
        TRUE                                                   ~ "#D4AC0D"
      ),
      border_col = dplyr::case_when(
        concordance == "concordant_depleted"                   ~ "#08306B",
        concordance == "concordant_enriched"                   ~ "#67000D",
        concordance == "single_assay" & hit_dir == "depleted"  ~ gof_border,
        concordance == "single_assay" & hit_dir == "enriched"  ~ lof_border,
        TRUE                                                   ~ "#9A7D0A"
      )
    )
}


# ─────────────────────────────────────────────────────────────
# assay_summary()
# ─────────────────────────────────────────────────────────────
#' Print and return a concordance summary table.
#'
#' @param mg  Output of merge_assays().
#' @return A tibble with one row per gene (hits only), sorted by concordance.
assay_summary <- function(mg) {
  dir_cols <- paste0("hit_dir__",  mg$assay_names)
  lfc_cols <- paste0("best_lfc__", mg$assay_names)
  nlp_cols <- paste0("nlp__",      mg$assay_names)

  out <- mg$hits %>%
    select(gene, concordance, hit_dir, nlp, best_lfc,
           all_of(dir_cols), all_of(lfc_cols), all_of(nlp_cols)) %>%
    arrange(concordance, desc(nlp))

  cat("\n── Merged assay concordance summary ─────────────────────\n")
  tbl <- out %>%
    group_by(concordance) %>%
    summarise(n = n(), genes = paste(head(gene, 6), collapse = ", "), .groups = "drop")
  print(tbl, n = Inf)

  # Discordant genes (worth inspecting)
  disc <- out %>% filter(concordance == "discordant")
  if (nrow(disc) > 0) {
    cat(sprintf("\n%d discordant genes (direction differs between assays):\n", nrow(disc)))
    print(disc %>% select(gene, all_of(dir_cols), nlp), n = Inf)
  }

  invisible(out)
}


# ─────────────────────────────────────────────────────────────
# plot_overlap()
# ─────────────────────────────────────────────────────────────
#' Quick UpSet / Venn-style bar chart of hit overlap between assays.
#'
#' @param mg  Output of merge_assays().
#' @return A ggplot object (bar chart of concordance × direction).
plot_overlap <- function(mg) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))

  count_tbl <- mg$hits %>%
    count(concordance, hit_dir) %>%
    mutate(
      label     = sprintf("%d", n),
      fill_col  = case_when(
        hit_dir == "depleted" & grepl("concordant|single.*dep", concordance) ~ "#2166AC",
        hit_dir == "enriched" & grepl("concordant|single.*enr", concordance) ~ "#D6604D",
        hit_dir == "depleted" ~ "#6BAED6",
        hit_dir == "enriched" ~ "#FC8D59",
        TRUE ~ "#D4AC0D"),
      concordance = factor(concordance,
        levels = c("concordant_depleted","concordant_enriched",
                   "single_assay","discordant","conflict"))
    )

  ggplot(count_tbl, aes(x = concordance, y = n, fill = fill_col)) +
    geom_col(width = 0.6, colour = "#333333", linewidth = 0.3) +
    geom_text(aes(label = label), vjust = -0.4, size = 3.2, fontface = "bold") +
    scale_fill_identity() +
    scale_x_discrete(labels = c(
      concordant_depleted = "Concordant\nGoF (depleted)",
      concordant_enriched = "Concordant\nLoF (enriched)",
      single_assay        = "Single-assay\nonly",
      discordant          = "Discordant",
      conflict            = "Tied / conflict")) +
    labs(
      title = sprintf("CRISPRi/CRISPRko hit overlap: %s",
                      paste(mg$assay_names, collapse = " vs ")),
      x = NULL, y = "Number of genes") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(size = 9))
}


# ─────────────────────────────────────────────────────────────
# build_dotmatrix()
# ─────────────────────────────────────────────────────────────
#' Build a long-format per-(gene × assay) table for the dot-matrix overlay.
#'
#' Returns one row per placed-node × assay. Each row carries the x/y
#' position of one small indicator dot drawn just below the node.
#' Dot fill colour encodes direction in that assay; position (left→right)
#' encodes assay identity.
#'
#' @param mg     Output of merge_assays().
#' @param nodes  The `nodes` tibble from the main script (must have x, y, gene).
#' @param y_offset Vertical offset below node centre (default −0.42).
#' @param x_step  Horizontal spacing between dots for N assays (default 0.22).
#'
#' @return A tibble with columns:
#'   gene, assay, assay_i, x_dot, y_dot,
#'   dot_fill, dot_border, dot_dir, dot_nlp
build_dotmatrix <- function(mg,
                            nodes,
                            y_offset = -0.42,
                            x_step   =  0.22) {
  assay_names <- mg$assay_names
  n_assays    <- length(assay_names)

  # Direction → colour lookup
  dir_pal <- tibble::tribble(
    ~dot_dir,    ~dot_fill,  ~dot_border,
    "depleted",  "#2166AC",  "#08306B",    # blue  – GoF
    "enriched",  "#D6604D",  "#67000D",    # red   – LoF
    "none",      "#AAAAAA",  "#777777"     # grey  – not a hit in this assay
  )

  purrr::imap_dfr(assay_names, function(aname, idx) {
    dir_col <- paste0("hit_dir__", aname)
    nlp_col <- paste0("nlp__",     aname)

    mg$dat %>%
      dplyr::select(
        gene,
        dot_dir = dplyr::all_of(dir_col),
        dot_nlp = dplyr::all_of(nlp_col)
      ) %>%
      dplyr::mutate(
        dot_dir = tidyr::replace_na(dot_dir, "none"),
        dot_nlp = tidyr::replace_na(dot_nlp, 0),
        assay   = aname,
        assay_i = idx
      ) %>%
      # Keep only nodes that are physically placed on the canvas
      dplyr::inner_join(
        nodes %>% dplyr::select(gene, x, y),
        by = "gene"
      ) %>%
      dplyr::mutate(
        # Centre the dot group horizontally under the node
        x_dot = x + (idx - (n_assays + 1L) / 2) * x_step,
        y_dot = y + y_offset
      )
  }) %>%
  dplyr::left_join(dir_pal, by = "dot_dir")
}


# ─────────────────────────────────────────────────────────────
# dotmatrix_legend_items()
# ─────────────────────────────────────────────────────────────
#' Return annotation parameters for a dot-matrix legend strip.
#'
#' Use the returned list inside an `annotate()` call in the main ggplot
#' to label the per-assay dot columns.
#'
#' @param mg     Output of merge_assays().
#' @param xl     x-coordinate of the legend left edge.
#' @param yl     y-coordinate of the legend row.
#' @param x_step Dot spacing (should match build_dotmatrix x_step).
#'
#' @return A list of ggplot annotation layers (annotate point + text).
dotmatrix_legend_items <- function(mg, xl = 22.0, yl = -1.15, x_step = 0.22) {
  assay_names <- mg$assay_names
  n_assays    <- length(assay_names)

  items <- list(
    ggplot2::annotate("text", x = xl - 0.15, y = yl,
                      label = "Screens:", size = 1.85,
                      colour = "#444444", hjust = 1, fontface = "bold")
  )

  for (i in seq_along(assay_names)) {
    xd      <- xl + (i - 1L) * (x_step + 0.9)
    pal_row <- ASSAY_PALETTE[((i - 1L) %% nrow(ASSAY_PALETTE)) + 1L, ]
    items <- c(items, list(
      ggplot2::annotate("point", x = xd, y = yl,
                        colour = pal_row$gof_border, fill = pal_row$gof_fill,
                        size = 2.2, shape = 21, stroke = 0.5),
      ggplot2::annotate("text", x = xd + 0.12, y = yl,
                        label = assay_names[i], size = 1.75,
                        colour = "#333333", hjust = 0)
    ))
  }

  # Colour key for dot fill
  items <- c(items, list(
    ggplot2::annotate("point", x = xl,            y = yl - 0.28,
                      colour = "#08306B", fill = "#2166AC",
                      size = 2.2, shape = 21, stroke = 0.5),
    ggplot2::annotate("text",  x = xl + 0.12,     y = yl - 0.28,
                      label = "GoF in screen", size = 1.75,
                      colour = "#333333", hjust = 0),
    ggplot2::annotate("point", x = xl + 1.8,      y = yl - 0.28,
                      colour = "#67000D", fill = "#D6604D",
                      size = 2.2, shape = 21, stroke = 0.5),
    ggplot2::annotate("text",  x = xl + 1.92,     y = yl - 0.28,
                      label = "LoF in screen", size = 1.75,
                      colour = "#333333", hjust = 0),
    ggplot2::annotate("point", x = xl + 3.6,      y = yl - 0.28,
                      colour = "#777777", fill = "#AAAAAA",
                      size = 2.2, shape = 21, stroke = 0.5),
    ggplot2::annotate("text",  x = xl + 3.72,     y = yl - 0.28,
                      label = "Not a hit",    size = 1.75,
                      colour = "#333333", hjust = 0)
  ))

  items
}
