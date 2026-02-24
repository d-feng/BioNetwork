############################################################
# TCGA-PAAD Survival Analysis v2 (Consolidated)
#
# Data sources (switchable via DATA_SOURCE setting):
#   TCGAbiolinks:     GDC STAR-Counts → tpm_unstrand → log2(TPM+1)
#   curatedTCGAData:  Bioc MAE → RNASeq2GeneNorm (RSEM UQ) → log2(RSEM+1)
#
# Validated against:
#   UCSC Xena star_tpm:  exact match (TCGAbiolinks tpm_unstrand)
#   TIMER2:              close match (uses counts-based, HR~1.289 vs TPM HR~1.330)
#   KMplot Q1 vs Q4:     different model (categorical, HR~1.83)
#
# Pipeline:
#   1-3) Load expression + clinical (dispatched by DATA_SOURCE)
#   4) Merge ESTIMATE TME scores
#   5) Helper: extract HR, CI, LRT/Wald/Log-rank p-values
#   6) Build gene-TME comparison table (unadj vs adjusted)
#   7) Gene survival + TME-adjusted Cox models
#   8) KM plots (Nature-style, log-rank p on plot)
#   9) Forest plot (LRT p-values)
#  10) Stratified Cox with strata() for TME
#  11) Nomogram system (C-index, calibration, time-ROC)
#  12) Run all analyses
############################################################

# ==========================================================
# 0) Packages
# ==========================================================
suppressPackageStartupMessages({
  pkgs <- c(
    "TCGAbiolinks", "SummarizedExperiment", "dplyr", "tibble", "stringr",
    "survival", "survminer", "maxstat", "ggplot2", "readr",
    "rms", "Hmisc", "timeROC"
  )
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) BiocManager::install("TCGAbiolinks")
  for (pkg in c("curatedTCGAData", "MultiAssayExperiment", "TCGAutils")) {
    if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  }
  if (!requireNamespace("estimate", quietly = TRUE)) {
    install.packages("estimate", repos = "http://R-Forge.R-project.org")
  }

  lapply(pkgs, library, character.only = TRUE)
})

# ==========================================================
# USER SETTINGS
# ==========================================================
PROJECT       <- "TCGA-PAAD"
DATA_SOURCE   <- "curatedTCGAData"        # "TCGAbiolinks" or "curatedTCGAData"
TPM_ASSAY     <- "tpm_unstrand"           # TCGAbiolinks only; validated: matches UCSC Xena star_tpm
ESTIMATE_CSV  <- "PAAD_Estimate.csv"     # auto-computed if not found
TME_COVARS    <- c("estimate_score", "stromal_score")
GENES_TO_TEST <- c("KRAS", "EGFR", "MAPK1", "MAPK3", "PIK3CA", "AKT1")
STRAT_METHOD  <- "quantile"
Q_CUTS        <- c(0.25, 0.75)
OUT_CSV       <- "PAAD_gene_unadj_vs_TMEadj.csv"

# ==========================================================
# 1-3) Data loaders: TCGAbiolinks or curatedTCGAData
#      Both return list(expr, surv_df, source_label)
#      expr:    numeric matrix, rownames=gene symbols, colnames=sample barcodes, log2-transformed
#      surv_df: tibble with sample, patient, vital_status, OS_time, OS_event, age_years, gender, stage_simple
# ==========================================================

# --- Loader A: TCGAbiolinks (GDC STAR-Counts) ---
load_from_tcgabiolinks <- function(project = PROJECT, tpm_assay = TPM_ASSAY) {

  cat("=== Loading data from TCGAbiolinks (GDC) ===\n")

  query <- GDCquery(
    project       = project,
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type  = "STAR - Counts"
  )
  GDCdownload(query)
  se <- GDCprepare(query)
  cat("Available assays:", paste(SummarizedExperiment::assayNames(se), collapse = ", "), "\n")

  if (!tpm_assay %in% SummarizedExperiment::assayNames(se)) {
    stop(sprintf("TPM_ASSAY='%s' not found. Available: %s",
                 tpm_assay, paste(SummarizedExperiment::assayNames(se), collapse = ", ")))
  }

  expr <- SummarizedExperiment::assay(se, tpm_assay)

  # Primary tumor only
  barcodes_tp <- TCGAquery_SampleTypes(colnames(expr), typesample = "TP")
  expr <- expr[, barcodes_tp, drop = FALSE]
  cat("Primary tumor samples:", ncol(expr), "\n")

  # Map gene symbols
  rowdat <- as.data.frame(SummarizedExperiment::rowData(se))
  if (!"gene_name" %in% colnames(rowdat)) stop("rowData(se)$gene_name not found.")
  rownames(expr) <- rowdat$gene_name

  # Collapse duplicate gene symbols by mean
  expr_df <- as.data.frame(expr) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(gene) & gene != "") %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean), .groups = "drop") %>%
    tibble::column_to_rownames("gene")
  expr <- as.matrix(expr_df)

  # log2(TPM + 1)
  expr <- log2(expr + 1)
  cat("Expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
  cat("Transform: log2(tpm_unstrand + 1)\n")

  # Clinical
  clin <- GDCquery_clinic(project = project, type = "clinical")
  clin2 <- clin %>%
    dplyr::mutate(
      patient   = dplyr::coalesce(submitter_id, bcr_patient_barcode),
      age_years = age_at_diagnosis / 365.25
    ) %>%
    dplyr::select(patient, vital_status, days_to_death, days_to_last_follow_up,
                  age_years, gender, ajcc_pathologic_stage) %>%
    dplyr::distinct(patient, .keep_all = TRUE)

  patient_id <- substr(colnames(expr), 1, 12)
  surv_df <- tibble::tibble(sample = colnames(expr), patient = patient_id) %>%
    dplyr::left_join(clin2, by = "patient") %>%
    dplyr::mutate(
      OS_time  = dplyr::if_else(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
      OS_event = dplyr::if_else(tolower(vital_status) == "dead", 1L, 0L),
      gender   = factor(gender),
      stage_simple = dplyr::case_when(
        stringr::str_detect(ajcc_pathologic_stage, stringr::regex("^stage i\\b",   ignore_case = TRUE)) ~ "I",
        stringr::str_detect(ajcc_pathologic_stage, stringr::regex("^stage ii\\b",  ignore_case = TRUE)) ~ "II",
        stringr::str_detect(ajcc_pathologic_stage, stringr::regex("^stage iii\\b", ignore_case = TRUE)) ~ "III",
        stringr::str_detect(ajcc_pathologic_stage, stringr::regex("^stage iv\\b",  ignore_case = TRUE)) ~ "IV",
        TRUE ~ NA_character_
      ),
      stage_simple = factor(stage_simple, levels = c("I", "II", "III", "IV"))
    ) %>%
    dplyr::filter(!is.na(OS_time) & OS_time > 0)

  expr <- expr[, surv_df$sample, drop = FALSE]
  list(expr = expr, surv_df = surv_df, source_label = "GDC tpm_unstrand + log2(TPM+1)")
}

# --- Loader B: curatedTCGAData (MultiAssayExperiment) ---
load_from_curatedtcga <- function(disease_code = sub("TCGA-", "", PROJECT)) {

  cat("=== Loading data from curatedTCGAData (MAE) ===\n")

  requireNamespace("curatedTCGAData", quietly = TRUE)
  requireNamespace("MultiAssayExperiment", quietly = TRUE)
  requireNamespace("TCGAutils", quietly = TRUE)

  # Download MAE
  mae <- curatedTCGAData::curatedTCGAData(
    diseaseCode = disease_code,
    assays      = "RNASeq2GeneNorm",
    version     = "2.0.1",
    dry.run     = FALSE
  )
  cat("MAE experiments:", paste(names(MultiAssayExperiment::experiments(mae)), collapse = ", "), "\n")

  # Extract expression matrix
  assay_name <- grep("RNASeq2GeneNorm", names(MultiAssayExperiment::experiments(mae)), value = TRUE)[1]
  if (is.na(assay_name)) stop("RNASeq2GeneNorm assay not found in MAE.")
  se <- MultiAssayExperiment::experiments(mae)[[assay_name]]
  expr <- as.matrix(SummarizedExperiment::assay(se))
  cat("Raw expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")

  # Filter primary tumors (sample type code "01")
  is_primary <- TCGAutils::TCGAsampleSelect(colnames(expr), sampleCodes = "01")
  expr <- expr[, is_primary, drop = FALSE]
  cat("Primary tumor samples:", ncol(expr), "\n")

  # Gene symbols already in rownames; collapse duplicates by mean
  expr_df <- as.data.frame(expr) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(gene) & gene != "") %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean), .groups = "drop") %>%
    tibble::column_to_rownames("gene")
  expr <- as.matrix(expr_df)

  # log2(RSEM + 1)
  expr <- log2(expr + 1)
  cat("Expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
  cat("Transform: log2(RSEM_UQ + 1)\n")
  cat("NOTE: RSEM upper-quartile normalization differs from STAR-Counts TPM.\n")
  cat("      Survival directions should agree; exact HR values may differ.\n")

  # --- Clinical from colData(mae) ---
  clin_raw <- as.data.frame(MultiAssayExperiment::colData(mae))
  cat("colData columns:", ncol(clin_raw), "\n")

  # Defensive column detection
  col_names <- colnames(clin_raw)

  # Vital status: numeric (0=alive, 1=dead) in curatedTCGAData
  vs_col <- intersect(c("vital_status"), col_names)[1]
  if (is.na(vs_col)) stop("vital_status column not found in colData. Available: ",
                           paste(head(col_names, 30), collapse = ", "))

  # Days to death / follow-up
  dtd_col <- intersect(c("days_to_death"), col_names)[1]
  dtf_col <- intersect(c("days_to_last_followup", "days_to_last_follow_up"), col_names)[1]

  # Age (in years for curatedTCGAData, in days for GDC)
  age_col <- intersect(c("age_at_initial_pathologic_diagnosis",
                          "age_at_diagnosis"), col_names)[1]

  # Gender
  gender_col <- intersect(c("gender"), col_names)[1]

  # Stage
  stage_col <- intersect(c("pathologic_stage", "ajcc_pathologic_stage",
                            "clinical_stage"), col_names)[1]

  cat("Mapped clinical columns:\n")
  cat("  vital_status:", vs_col, "\n")
  cat("  days_to_death:", dtd_col, "\n")
  cat("  days_to_followup:", dtf_col, "\n")
  cat("  age:", age_col, "\n")
  cat("  gender:", gender_col, "\n")
  cat("  stage:", stage_col, "\n")

  # Safe column extraction (returns NA vector if column missing)
  safe_col <- function(df, col_name) {
    if (is.na(col_name) || !col_name %in% colnames(df)) return(rep(NA_real_, nrow(df)))
    as.numeric(df[[col_name]])
  }
  safe_col_chr <- function(df, col_name) {
    if (is.na(col_name) || !col_name %in% colnames(df)) return(rep(NA_character_, nrow(df)))
    as.character(df[[col_name]])
  }

  clin2 <- clin_raw %>%
    tibble::rownames_to_column("patient") %>%
    dplyr::mutate(
      .vs     = as.numeric(.data[[vs_col]]),
      .dtd    = safe_col(clin_raw, dtd_col),
      .dtf    = safe_col(clin_raw, dtf_col),
      .age    = safe_col(clin_raw, age_col),
      .gender = safe_col_chr(clin_raw, gender_col),
      .stage  = safe_col_chr(clin_raw, stage_col)
    ) %>%
    dplyr::select(patient, .vs, .dtd, .dtf, .age, .gender, .stage) %>%
    dplyr::distinct(patient, .keep_all = TRUE)

  # Detect age unit: curatedTCGAData stores age in years (typically < 120)
  # GDC stores in days (typically > 1000)
  age_values <- na.omit(clin2$.age)
  age_in_days <- length(age_values) > 0 && stats::median(age_values) > 200
  if (age_in_days) {
    cat("  Age appears to be in days, converting to years\n")
    clin2$.age <- clin2$.age / 365.25
  }

  # Build surv_df
  patient_id <- substr(colnames(expr), 1, 12)
  surv_df <- tibble::tibble(sample = colnames(expr), patient = patient_id) %>%
    dplyr::left_join(clin2, by = "patient") %>%
    dplyr::mutate(
      OS_time  = dplyr::if_else(!is.na(.dtd), .dtd, .dtf),
      OS_event = dplyr::if_else(.vs == 1, 1L, 0L),
      vital_status = dplyr::if_else(.vs == 1, "Dead", "Alive"),
      age_years = .age,
      gender   = factor(.gender),
      stage_simple = dplyr::case_when(
        stringr::str_detect(.stage, stringr::regex("^stage i[ab]?$",   ignore_case = TRUE)) ~ "I",
        stringr::str_detect(.stage, stringr::regex("^stage ii[abc]?$", ignore_case = TRUE)) ~ "II",
        stringr::str_detect(.stage, stringr::regex("^stage iii[abc]?$",ignore_case = TRUE)) ~ "III",
        stringr::str_detect(.stage, stringr::regex("^stage iv[abc]?$", ignore_case = TRUE)) ~ "IV",
        TRUE ~ NA_character_
      ),
      stage_simple = factor(stage_simple, levels = c("I", "II", "III", "IV"))
    ) %>%
    dplyr::select(-dplyr::starts_with(".")) %>%
    dplyr::filter(!is.na(OS_time) & OS_time > 0)

  expr <- expr[, surv_df$sample, drop = FALSE]
  list(expr = expr, surv_df = surv_df,
       source_label = "curatedTCGAData RSEM_UQ + log2(RSEM+1)")
}

# --- Dispatcher ---
load_tcga_data <- function(data_source = DATA_SOURCE, ...) {

  result <- switch(data_source,
    "TCGAbiolinks"    = load_from_tcgabiolinks(...),
    "curatedTCGAData" = load_from_curatedtcga(...),
    stop("Unknown DATA_SOURCE: '", data_source, "'. Use 'TCGAbiolinks' or 'curatedTCGAData'.")
  )

  expr    <- result$expr
  surv_df <- result$surv_df

  # Validate downstream contract
  stopifnot(is.matrix(expr), is.numeric(expr))
  stopifnot(all(surv_df$sample %in% colnames(expr)))
  required_cols <- c("sample", "patient", "OS_time", "OS_event",
                     "age_years", "gender", "stage_simple")
  missing_cols <- setdiff(required_cols, names(surv_df))
  if (length(missing_cols) > 0) stop("surv_df missing columns: ", paste(missing_cols, collapse = ", "))

  cat("\n", strrep("=", 50), "\n")
  cat("  DATA LOADING SUMMARY\n")
  cat(strrep("=", 50), "\n")
  cat("  Source:         ", data_source, "\n")
  cat("  Label:          ", result$source_label, "\n")
  cat("  Genes:          ", nrow(expr), "\n")
  cat("  Samples:        ", ncol(expr), "\n")
  cat("  Samples w/ OS:  ", nrow(surv_df), "\n")
  cat("  OS events:      ", sum(surv_df$OS_event == 1), "\n")
  cat("  Barcode example:", surv_df$sample[1], "\n")
  cat("  Expr range:      [", sprintf("%.2f", min(expr)), ",",
      sprintf("%.2f", max(expr)), "]\n")
  cat("  Stage dist:     ", paste(names(table(surv_df$stage_simple)),
                                   table(surv_df$stage_simple), sep = "=", collapse = ", "), "\n")
  cat(strrep("=", 50), "\n\n")

  result
}

# ==========================================================
# 1-3) Execute data loading
# ==========================================================
.data  <- load_tcga_data()
expr    <- .data$expr
surv_df <- .data$surv_df
.source_label <- .data$source_label
rm(.data)

cat("PAAD tumor samples with OS:", nrow(surv_df), "\n")
cat("Events:", sum(surv_df$OS_event == 1), "\n")

# ==========================================================
# 4) Compute ESTIMATE TME scores from expression matrix
#    Uses the 'estimate' R package (Yoshihara et al. 2013)
#    Input:  expr (log2-transformed, gene symbols as rownames)
#    Output: stromal_score, immune_score, estimate_score added to surv_df
#
#    If ESTIMATE_CSV exists, load from file instead (backwards-compatible)
# ==========================================================
compute_estimate_scores <- function(expr_mat, samples) {
  requireNamespace("estimate", quietly = TRUE)

  # Write expression matrix to GCT format (required by estimate package)
  gct_in  <- tempfile(fileext = ".gct")
  gct_flt <- tempfile(fileext = ".gct")
  gct_out <- tempfile(fileext = ".gct")
  on.exit(unlink(c(gct_in, gct_flt, gct_out)), add = TRUE)

  n_genes   <- nrow(expr_mat)
  n_samples <- ncol(expr_mat)

  cat("  Writing GCT file:", n_genes, "genes x", n_samples, "samples\n")

  con <- file(gct_in, "w")
  writeLines("#1.2", con)
  writeLines(paste(n_genes, n_samples, sep = "\t"), con)
  header <- paste(c("Name", "Description", colnames(expr_mat)), collapse = "\t")
  writeLines(header, con)
  for (i in seq_len(n_genes)) {
    row_data <- paste(c(rownames(expr_mat)[i], "na",
                        as.character(expr_mat[i, ])), collapse = "\t")
    writeLines(row_data, con)
  }
  close(con)

  # Filter to ESTIMATE common genes (~10,412)
  cat("  Filtering common genes...\n")
  estimate::filterCommonGenes(
    input.f  = gct_in,
    output.f = gct_flt,
    id       = "GeneSymbol"
  )

  # Compute ESTIMATE scores
  cat("  Computing ESTIMATE scores...\n")
  estimate::estimateScore(
    input.ds  = gct_flt,
    output.ds = gct_out,
    platform  = "illumina"    # RNA-seq: use "illumina" (no purity calc)
  )

  # Read scores from GCT output
  score_lines <- readLines(gct_out)
  # GCT: line 1 = #1.2, line 2 = dims, line 3+ = header + data
  score_raw <- utils::read.delim(text = paste(score_lines[-(1:2)], collapse = "\n"),
                                  stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  score_raw <- score_raw[, -1, drop = FALSE]   # drop Description column

  # Transpose: rows = score names, cols = samples → we want samples as rows
  score_t <- as.data.frame(t(score_raw), stringsAsFactors = FALSE)
  score_t$sample <- rownames(score_t)
  rownames(score_t) <- NULL

  # Normalize column names
  names(score_t) <- names(score_t) %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("_+$", "") %>%
    tolower()

  # Ensure numeric
  for (col in c("stromalscore", "immunescore", "estimatescore", "tumorpurity")) {
    if (col %in% names(score_t)) score_t[[col]] <- as.numeric(score_t[[col]])
  }

  # Standardize column names
  rename_map <- c(
    "stromal_score"  = "stromalscore",
    "immune_score"   = "immunescore",
    "estimate_score" = "estimatescore",
    "tumor_purity"   = "tumorpurity"
  )
  for (nm in names(rename_map)) {
    if (rename_map[[nm]] %in% names(score_t) && !(nm %in% names(score_t))) {
      score_t <- score_t %>% dplyr::rename(!!nm := !!rename_map[[nm]])
    }
  }

  cat("  ESTIMATE scores computed for", nrow(score_t), "samples\n")
  score_t
}

# --- Compute or load ESTIMATE scores ---
if (file.exists(ESTIMATE_CSV)) {
  cat("Loading ESTIMATE scores from file:", ESTIMATE_CSV, "\n")
  score_raw <- readr::read_csv(ESTIMATE_CSV, show_col_types = FALSE)
  score <- score_raw
  names(score) <- names(score) %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("_+$", "") %>%
    tolower()
  id_col <- intersect(c("id", "sample", "barcode", "tcga_id"), names(score))[1]
  if (is.na(id_col)) stop("Could not find ID column in ESTIMATE CSV.")
  score <- score %>% dplyr::rename(sample = !!id_col)
  rename_map <- c("stromal_score" = "stromalscore", "immune_score" = "immunescore",
                   "estimate_score" = "estimatescore", "tumor_purity" = "tumorpurity")
  for (nm in names(rename_map)) {
    if (rename_map[[nm]] %in% names(score) && !(nm %in% names(score))) {
      score <- score %>% dplyr::rename(!!nm := !!rename_map[[nm]])
    }
  }
} else {
  cat("No", ESTIMATE_CSV, "found. Computing ESTIMATE scores from expression matrix...\n")
  score <- compute_estimate_scores(expr, surv_df$sample)
}

# Merge on 15-char barcode (handles both 16-char curatedTCGAData and 28-char GDC barcodes)
surv_df <- surv_df %>%
  dplyr::mutate(sample_key = substr(sample, 1, 15))

score2 <- score %>%
  dplyr::mutate(sample_key = substr(sample, 1, 15))

surv_df <- surv_df %>%
  dplyr::left_join(
    score2 %>% dplyr::select(sample_key,
      dplyr::any_of(c("stromal_score", "immune_score", "estimate_score"))),
    by = "sample_key"
  ) %>%
  dplyr::select(-sample_key)

cat("NA estimate_score:", sum(is.na(surv_df$estimate_score)), "\n")
cat("NA stromal_score :", sum(is.na(surv_df$stromal_score)), "\n")

# Save computed scores for future use
if (!file.exists(ESTIMATE_CSV)) {
  score_save <- score %>%
    dplyr::select(sample, dplyr::any_of(c("stromal_score", "immune_score",
                                            "estimate_score", "tumor_purity")))
  readr::write_csv(score_save, ESTIMATE_CSV)
  cat("Saved computed scores to:", ESTIMATE_CSV, "\n")
}

# ==========================================================
# 5) Helper: extract HR, CI, three p-values (LRT, Wald, Log-rank)
#    LRT   = robust for small samples (N~88)
#    Wald  = standard Cox output
#    Log-rank = non-parametric, unadjusted only
# ==========================================================
extract_groupHigh_stats <- function(fit, dat = NULL) {
  na_out <- c(HR = NA_real_, L = NA_real_, U = NA_real_,
              p_lrt = NA_real_, p_wald = NA_real_, p_logrank = NA_real_)
  if (is.null(fit)) return(na_out)

  s  <- summary(fit)
  rn <- rownames(s$coefficients)

  idx <- which(rn == "groupHigh")
  if (length(idx) == 0) idx <- grep("^group", rn)
  if (length(idx) == 0) return(na_out)
  idx <- idx[1]

  beta <- s$coefficients[idx, "coef"]
  se   <- s$coefficients[idx, "se(coef)"]
  HR   <- exp(beta)
  L    <- exp(beta - 1.96 * se)
  U    <- exp(beta + 1.96 * se)

  # Wald p
  p_wald <- as.numeric(s$coefficients[idx, "Pr(>|z|)"])

  # LRT p (robust for small samples)
  p_lrt <- tryCatch({
    d1 <- drop1(fit, test = "Chisq")
    group_row <- grep("^group$", rownames(d1))
    if (length(group_row) > 0) {
      as.numeric(d1[group_row, "Pr(>Chi)"])
    } else {
      as.numeric(s$logtest["pvalue"])
    }
  }, error = function(e) {
    as.numeric(s$logtest["pvalue"])
  })

  # Log-rank p (unadjusted only)
  p_logrank <- NA_real_
  if (!is.null(dat) && nrow(s$coefficients) == 1) {
    p_logrank <- tryCatch({
      sdiff <- survival::survdiff(survival::Surv(OS_time, OS_event) ~ group, data = dat)
      as.numeric(1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1))
    }, error = function(e) NA_real_)
  }

  c(HR = as.numeric(HR), L = as.numeric(L), U = as.numeric(U),
    p_lrt = as.numeric(p_lrt), p_wald = as.numeric(p_wald),
    p_logrank = as.numeric(p_logrank))
}

# ==========================================================
# 6) Build gene-TME comparison table (unadj vs adjusted)
#    No clinical metadata — only TME scores as covariates
# ==========================================================
build_gene_TME_compare_table <- function(genes,
                                         q       = c(0.25, 0.75),
                                         out_csv = "PAAD_gene_unadj_vs_TMEadj.csv") {

  genes <- unique(genes)
  genes_present <- genes[genes %in% rownames(expr)]
  if (length(genes_present) == 0) stop("None of the genes found in expr.")

  results <- vector("list", length(genes_present))
  names(results) <- genes_present

  for (g in genes_present) {

    dat <- surv_df %>%
      dplyr::mutate(gene_expr = as.numeric(expr[g, sample])) %>%
      dplyr::filter(!is.na(gene_expr))

    qs       <- stats::quantile(dat$gene_expr, probs = q, na.rm = TRUE)
    cut_low  <- as.numeric(qs[1])
    cut_high <- as.numeric(qs[2])

    dat <- dat %>%
      dplyr::mutate(group = dplyr::case_when(
        gene_expr <= cut_low  ~ "Low",
        gene_expr >= cut_high ~ "High",
        TRUE ~ NA_character_
      )) %>%
      dplyr::filter(!is.na(group))

    dat$group <- factor(dat$group, levels = c("Low", "High"))

    n_total  <- nrow(dat)
    n_events <- sum(dat$OS_event == 1, na.rm = TRUE)

    # Cox models: unadjusted, + estimate_score, + stromal_score
    cox_unadj <- tryCatch(
      survival::coxph(survival::Surv(OS_time, OS_event) ~ group, data = dat),
      error = function(e) NULL)
    cox_est <- tryCatch(
      survival::coxph(survival::Surv(OS_time, OS_event) ~ group + estimate_score, data = dat),
      error = function(e) NULL)
    cox_str <- tryCatch(
      survival::coxph(survival::Surv(OS_time, OS_event) ~ group + stromal_score, data = dat),
      error = function(e) NULL)

    unadj <- extract_groupHigh_stats(cox_unadj, dat = dat)
    est   <- extract_groupHigh_stats(cox_est)
    str   <- extract_groupHigh_stats(cox_str)

    results[[g]] <- data.frame(
      Gene = g, Method = "quantile", N_total = n_total, Events = n_events,
      Cutoff_low = cut_low, Cutoff_high = cut_high,
      HR_unadj = unadj["HR"], CI_low_unadj = unadj["L"], CI_high_unadj = unadj["U"],
      p_lrt_unadj = unadj["p_lrt"], p_wald_unadj = unadj["p_wald"],
      p_logrank_unadj = unadj["p_logrank"],
      HR_adj_estimate = est["HR"], CI_low_adj_estimate = est["L"],
      CI_high_adj_estimate = est["U"],
      p_lrt_adj_estimate = est["p_lrt"], p_wald_adj_estimate = est["p_wald"],
      HR_adj_stromal = str["HR"], CI_low_adj_stromal = str["L"],
      CI_high_adj_stromal = str["U"],
      p_lrt_adj_stromal = str["p_lrt"], p_wald_adj_stromal = str["p_wald"],
      stringsAsFactors = FALSE, row.names = NULL
    )
  }

  out <- dplyr::bind_rows(results)
  readr::write_csv(out, out_csv)
  cat("Saved:", out_csv, "\n")
  out
}

# ==========================================================
# 7) Run gene survival + TME-adjusted Cox
# ==========================================================
run_gene_survival_TME <- function(gene,
                                  q              = c(0.25, 0.75),
                                  covar_estimate = "estimate_score",
                                  covar_stromal  = "stromal_score") {

  if (!gene %in% rownames(expr)) stop("Gene not found: ", gene)

  dat <- surv_df %>%
    dplyr::mutate(gene_expr = as.numeric(expr[gene, sample])) %>%
    dplyr::filter(!is.na(gene_expr), !is.na(OS_time), OS_time > 0)

  if (nrow(dat) < 10) stop("Too few samples.")

  qs          <- stats::quantile(dat$gene_expr, probs = q, na.rm = TRUE)
  cutoff_low  <- as.numeric(qs[1])
  cutoff_high <- as.numeric(qs[2])

  dat <- dat %>%
    dplyr::mutate(group = dplyr::case_when(
      gene_expr <= cutoff_low  ~ "Low",
      gene_expr >= cutoff_high ~ "High",
      TRUE ~ NA_character_
    )) %>%
    dplyr::filter(!is.na(group))

  dat$group <- factor(dat$group, levels = c("Low", "High"))

  if (nrow(dat) < 10 || length(unique(dat$group)) < 2 ||
      sum(dat$OS_event == 1) < 1) stop("Insufficient data after stratification.")

  cox_unadj <- survival::coxph(survival::Surv(OS_time, OS_event) ~ group, data = dat)

  ok_covar <- function(x) {
    v <- dat[[x]]; if (is.null(v) || all(is.na(v))) return(FALSE)
    length(unique(na.omit(v))) >= 2
  }

  cox_adj_estimate <- if (ok_covar(covar_estimate)) {
    survival::coxph(
      stats::as.formula(paste0("survival::Surv(OS_time, OS_event) ~ group + ", covar_estimate)),
      data = dat)
  } else NULL

  cox_adj_stromal <- if (ok_covar(covar_stromal)) {
    survival::coxph(
      stats::as.formula(paste0("survival::Surv(OS_time, OS_event) ~ group + ", covar_stromal)),
      data = dat)
  } else NULL

  km_fit <- survival::survfit(survival::Surv(OS_time, OS_event) ~ group, data = dat)

  list(gene = gene, method = "quantile",
       cutoff_low = cutoff_low, cutoff_high = cutoff_high,
       data = dat, km_fit = km_fit,
       cox_unadj = cox_unadj,
       cox_adj_estimate = cox_adj_estimate,
       cox_adj_stromal = cox_adj_stromal)
}

# ==========================================================
# 8) KM plot (Nature-style, log-rank p on plot)
# ==========================================================
plot_gene_km <- function(gene, q = c(0.25, 0.75),
                         out_png = NULL, dpi = 600, show_hr = TRUE) {

  res <- run_gene_survival_TME(gene, q = q)
  dat <- res$data

  title_txt <- paste0(gene, " (Q", q[1]*100, " vs Q", q[2]*100, ")")

  # HR + p-values for subtitle
  extract_group_hr <- function(fit, fit_dat = NULL) {
    if (is.null(fit)) return(NULL)
    s <- summary(fit); rn <- rownames(s$coefficients)
    idx <- which(grepl("^groupHigh$", rn))
    if (!length(idx)) return(NULL)
    beta <- s$coefficients[idx, "coef"]; se <- s$coefficients[idx, "se(coef)"]
    hr <- exp(beta); ci_l <- exp(beta - 1.96 * se); ci_u <- exp(beta + 1.96 * se)
    p_wald <- as.numeric(s$coefficients[idx, "Pr(>|z|)"])
    p_lrt <- tryCatch({
      d1 <- drop1(fit, test = "Chisq")
      gr <- grep("^group$", rownames(d1))
      if (length(gr) > 0) as.numeric(d1[gr, "Pr(>Chi)"]) else as.numeric(s$logtest["pvalue"])
    }, error = function(e) as.numeric(s$logtest["pvalue"]))
    p_logrank <- NA_real_
    if (!is.null(fit_dat) && nrow(s$coefficients) == 1) {
      p_logrank <- tryCatch({
        sdiff <- survival::survdiff(survival::Surv(OS_time, OS_event) ~ group, data = fit_dat)
        as.numeric(1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1))
      }, error = function(e) NA_real_)
    }
    list(HR = hr, L = ci_l, U = ci_u, p_lrt = p_lrt, p_wald = p_wald, p_logrank = p_logrank)
  }

  sub_txt <- NULL
  if (show_hr) {
    u <- extract_group_hr(res$cox_unadj, fit_dat = dat)
    a <- extract_group_hr(res$cox_adj_estimate)
    b <- extract_group_hr(res$cox_adj_stromal)
    fmt <- function(x, label) {
      if (is.null(x)) return(paste0(label, ": NA"))
      p_txt <- paste0("LRT p=", signif(x$p_lrt, 3), ", Wald p=", signif(x$p_wald, 3))
      if (!is.na(x$p_logrank)) p_txt <- paste0(p_txt, ", Log-rank p=", signif(x$p_logrank, 3))
      paste0(label, ": HR=", sprintf("%.2f", x$HR),
             " (", sprintf("%.2f", x$L), "\u2013", sprintf("%.2f", x$U), "), ", p_txt)
    }
    sub_txt <- paste(fmt(u, "Unadj"), fmt(a, "Adj+Estimate"), fmt(b, "Adj+Stroma"), sep = "\n")
  }

  # Log-rank p for annotation
  logrank_p <- tryCatch({
    sdiff <- survival::survdiff(survival::Surv(OS_time, OS_event) ~ group, data = dat)
    as.numeric(1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1))
  }, error = function(e) NA_real_)

  fit <- survival::survfit(survival::Surv(OS_time, OS_event) ~ group, data = dat)

  g <- survminer::ggsurvplot(
    fit, data = dat,
    palette = c("black", "grey50"), linetype = c("solid", "dashed"),
    conf.int = FALSE, pval = FALSE,
    risk.table = TRUE, risk.table.col = "black",
    ggtheme = ggplot2::theme_classic(base_size = 16),
    title = title_txt, subtitle = sub_txt,
    legend.title = "", legend.labs = c("Low", "High"),
    xlab = "Time (days)", ylab = "Overall survival probability"
  )

  g$plot <- g$plot +
    ggplot2::annotate("text",
      x = max(dat$OS_time, na.rm = TRUE) * 0.55, y = 0.95,
      label = paste0("Log-rank p = ", signif(logrank_p, 3)),
      size = 5, hjust = 0) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.0),
      axis.line = ggplot2::element_line(linewidth = 0.9),
      axis.ticks = ggplot2::element_line(linewidth = 0.9),
      legend.position = "right")

  if (is.null(out_png)) out_png <- paste0("KM_", gene, "_PAAD.png")
  grDevices::png(out_png, width = 8, height = 7, units = "in", res = dpi)
  print(g)
  grDevices::dev.off()
  cat("Saved:", out_png, "\n")
}

# ==========================================================
# 9) Forest plot (LRT p-values)
# ==========================================================
plot_forest <- function(tme_table,
                        out_png = "Forest_PAAD_unadj_vs_TMEadj.png", dpi = 600) {

  forest_df <- dplyr::bind_rows(
    tme_table %>% dplyr::transmute(Gene, Model = "Unadjusted",
      HR = HR_unadj, CI_low = CI_low_unadj, CI_high = CI_high_unadj, p = p_lrt_unadj),
    tme_table %>% dplyr::transmute(Gene, Model = "Adj + ESTIMATE",
      HR = HR_adj_estimate, CI_low = CI_low_adj_estimate,
      CI_high = CI_high_adj_estimate, p = p_lrt_adj_estimate),
    tme_table %>% dplyr::transmute(Gene, Model = "Adj + Stromal",
      HR = HR_adj_stromal, CI_low = CI_low_adj_stromal,
      CI_high = CI_high_adj_stromal, p = p_lrt_adj_stromal)
  ) %>% dplyr::filter(!is.na(HR))

  forest_df$Gene  <- factor(forest_df$Gene, levels = rev(unique(tme_table$Gene)))
  forest_df$Model <- factor(forest_df$Model,
    levels = c("Unadjusted", "Adj + ESTIMATE", "Adj + Stromal"))

  forest_df <- forest_df %>%
    dplyr::mutate(
      sig_label = dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""),
      hr_label = paste0(sprintf("%.2f", HR),
        " (", sprintf("%.2f", CI_low), "\u2013", sprintf("%.2f", CI_high), ") ", sig_label))

  pd <- ggplot2::position_dodge(width = 0.6)

  p <- ggplot2::ggplot(forest_df, ggplot2::aes(x = HR, y = Gene, color = Model, shape = Model)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    ggplot2::geom_point(position = pd, size = 3.5) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low, xmax = CI_high),
      position = pd, height = 0.2, linewidth = 0.7) +
    ggplot2::geom_text(ggplot2::aes(x = CI_high + 0.05, label = hr_label),
      position = pd, hjust = 0, size = 3, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c("Unadjusted" = "black",
      "Adj + ESTIMATE" = "#E7298A", "Adj + Stromal" = "#2E9FDF")) +
    ggplot2::scale_shape_manual(values = c("Unadjusted" = 16,
      "Adj + ESTIMATE" = 17, "Adj + Stromal" = 15)) +
    ggplot2::scale_x_continuous(trans = "log2",
      breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4")) +
    ggplot2::labs(title = "TCGA-PAAD: Cox HR (High vs Low, Q25 vs Q75)",
      subtitle = "LRT p-value | * p<0.05, ** p<0.01, *** p<0.001",
      x = "Hazard Ratio (log2 scale)", y = "", color = "Model", shape = "Model") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey30"),
      axis.text.y = ggplot2::element_text(face = "bold", size = 12)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 80, 10, 10))

  grDevices::png(out_png, width = 10, height = 6, units = "in", res = dpi)
  print(p)
  grDevices::dev.off()
  cat("Saved:", out_png, "\n")
}

# ==========================================================
# 10) Stratified Cox with strata() for TME
# ==========================================================
run_stratified_cox <- function(genes, q = c(0.25, 0.75),
                               tme_var = "estimate_score",
                               tme_split = "median",
                               out_csv = "PAAD_stratified_cox.csv") {

  genes <- unique(genes)
  genes_present <- genes[genes %in% rownames(expr)]
  if (length(genes_present) == 0) stop("No genes found.")

  results <- vector("list", length(genes_present))
  names(results) <- genes_present

  for (g in genes_present) {
    dat <- surv_df %>%
      dplyr::mutate(gene_expr = as.numeric(expr[g, sample])) %>%
      dplyr::filter(!is.na(gene_expr), !is.na(.data[[tme_var]]))

    qs <- stats::quantile(dat$gene_expr, probs = q, na.rm = TRUE)
    cut_low <- as.numeric(qs[1]); cut_high <- as.numeric(qs[2])

    dat <- dat %>%
      dplyr::mutate(group = dplyr::case_when(
        gene_expr <= cut_low ~ "Low", gene_expr >= cut_high ~ "High",
        TRUE ~ NA_character_)) %>%
      dplyr::filter(!is.na(group))
    dat$group <- factor(dat$group, levels = c("Low", "High"))

    tme_cut <- if (tme_split == "mean") mean(dat[[tme_var]], na.rm = TRUE) else
      stats::median(dat[[tme_var]], na.rm = TRUE)
    dat$tme_group <- factor(ifelse(dat[[tme_var]] >= tme_cut, "High", "Low"),
                            levels = c("Low", "High"))

    n_total <- nrow(dat); n_events <- sum(dat$OS_event == 1, na.rm = TRUE)

    # 4 models
    cox_unadj <- tryCatch(
      survival::coxph(survival::Surv(OS_time, OS_event) ~ group, data = dat),
      error = function(e) NULL)

    fml_adj <- stats::as.formula(paste0("survival::Surv(OS_time, OS_event) ~ group + ", tme_var))
    cox_adj <- tryCatch(survival::coxph(fml_adj, data = dat), error = function(e) NULL)

    cox_strata <- tryCatch(
      survival::coxph(survival::Surv(OS_time, OS_event) ~ group + survival::strata(tme_group),
                      data = dat),
      error = function(e) NULL)

    fml_int <- stats::as.formula(paste0("survival::Surv(OS_time, OS_event) ~ group * ", tme_var))
    cox_interact <- tryCatch(survival::coxph(fml_int, data = dat), error = function(e) NULL)

    extract_stats <- function(fit, model_label) {
      if (is.null(fit)) {
        return(data.frame(Gene = g, Model = model_label, N = n_total, Events = n_events,
          HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
          p_lrt = NA_real_, p_wald = NA_real_, stringsAsFactors = FALSE))
      }
      s <- summary(fit); rn <- rownames(s$coefficients)
      idx <- which(rn == "groupHigh")
      if (length(idx) == 0) idx <- grep("^group", rn)
      if (length(idx) == 0) {
        return(data.frame(Gene = g, Model = model_label, N = n_total, Events = n_events,
          HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
          p_lrt = NA_real_, p_wald = NA_real_, stringsAsFactors = FALSE))
      }
      idx <- idx[1]
      beta <- s$coefficients[idx, "coef"]; se <- s$coefficients[idx, "se(coef)"]
      HR <- exp(beta); L <- exp(beta - 1.96 * se); U <- exp(beta + 1.96 * se)
      p_wald <- as.numeric(s$coefficients[idx, "Pr(>|z|)"])
      p_lrt <- tryCatch({
        d1 <- drop1(fit, test = "Chisq"); gr <- grep("^group$", rownames(d1))
        if (length(gr) > 0) as.numeric(d1[gr, "Pr(>Chi)"]) else as.numeric(s$logtest["pvalue"])
      }, error = function(e) as.numeric(s$logtest["pvalue"]))
      data.frame(Gene = g, Model = model_label, N = n_total, Events = n_events,
        HR = HR, CI_low = L, CI_high = U, p_lrt = p_lrt, p_wald = p_wald,
        stringsAsFactors = FALSE)
    }

    p_interact <- NA_real_
    if (!is.null(cox_interact)) {
      s_int <- summary(cox_interact); rn_int <- rownames(s_int$coefficients)
      int_idx <- grep(":", rn_int)
      if (length(int_idx) > 0) p_interact <- as.numeric(s_int$coefficients[int_idx[1], "Pr(>|z|)"])
    }

    rows <- dplyr::bind_rows(
      extract_stats(cox_unadj,   "Unadjusted"),
      extract_stats(cox_adj,     paste0("Adjusted + ", tme_var)),
      extract_stats(cox_strata,  paste0("Stratified by ", tme_var)),
      extract_stats(cox_interact, paste0("Interaction (group*", tme_var, ")"))
    )
    rows$p_interaction <- c(NA_real_, NA_real_, NA_real_, p_interact)
    rows$TME_cutoff    <- c(NA_real_, NA_real_, tme_cut, NA_real_)
    rows$HR_fmt <- ifelse(is.na(rows$HR), "NA",
      paste0(sprintf("%.2f", rows$HR), " (", sprintf("%.2f", rows$CI_low),
             "\u2013", sprintf("%.2f", rows$CI_high), ")"))

    results[[g]] <- rows
  }

  out <- dplyr::bind_rows(results)
  readr::write_csv(out, out_csv)
  cat("\n=== Stratified Cox Summary ===\n")
  print(out %>% dplyr::select(Gene, Model, N, Events, HR_fmt, p_lrt, p_wald, p_interaction))
  cat("Saved:", out_csv, "\n")
  out
}

# ==========================================================
# 11) Nomogram system (C-index, calibration, time-ROC)
# ==========================================================

# 11a) Build nomogram data
build_nomogram_data <- function(genes, time_days_max = NULL) {
  genes <- genes[genes %in% rownames(expr)]
  if (length(genes) == 0) stop("No genes found.")
  nom_df <- surv_df
  for (g in genes) nom_df[[g]] <- as.numeric(expr[g, nom_df$sample])
  nom_df <- nom_df %>%
    dplyr::mutate(OS_months = OS_time / 30.44,
                  age_years = ifelse(is.na(age_years), NA_real_, age_years))
  required_cols <- c("OS_months", "OS_event", "age_years", "stage_simple",
                     "estimate_score", "stromal_score", genes)
  nom_df <- nom_df %>% dplyr::filter(dplyr::if_all(dplyr::all_of(required_cols), ~ !is.na(.)))
  if (!is.null(time_days_max)) {
    time_months_max <- time_days_max / 30.44
    nom_df <- nom_df %>% dplyr::mutate(
      OS_event = ifelse(OS_months > time_months_max, 0L, OS_event),
      OS_months = pmin(OS_months, time_months_max))
  }
  cat("Nomogram data: N =", nrow(nom_df), "samples,", sum(nom_df$OS_event == 1), "events\n")
  nom_df
}

# 11b) Fit rms::cph model
fit_nomogram_model <- function(nom_df, genes,
                               include_clinical = TRUE, include_tme = TRUE) {
  dd <- rms::datadist(nom_df); options(datadist = "dd"); assign("dd", dd, envir = .GlobalEnv)
  gene_terms <- paste(genes, collapse = " + ")
  clin_terms <- if (include_clinical) "age_years + stage_simple" else NULL
  tme_terms  <- if (include_tme) "estimate_score + stromal_score" else NULL
  all_terms <- paste(c(gene_terms, clin_terms, tme_terms), collapse = " + ")
  fml <- stats::as.formula(paste0("Surv(OS_months, OS_event) ~ ", all_terms))
  cat("Nomogram formula:", deparse(fml), "\n")
  fit <- rms::cph(fml, data = nom_df, x = TRUE, y = TRUE, surv = TRUE, time.inc = 12)
  cat("Model fitted. N =", fit$stats["Obs"], ", Events =", fit$stats["Events"], "\n")
  fit
}

# 11c) Plot nomogram
plot_nomogram <- function(fit, time_points = c(12, 24, 36),
                          out_png = "Nomogram_PAAD.png", dpi = 600) {
  surv_obj <- rms::Survival(fit)
  nom <- rms::nomogram(fit,
    fun = lapply(time_points, function(t) function(x) surv_obj(t, x)),
    fun.at = seq(0.1, 0.9, by = 0.1),
    funlabel = paste0(time_points / 12, "-Year OS"),
    lp = TRUE, maxscale = 100)
  grDevices::png(out_png, width = 12, height = 10, units = "in", res = dpi)
  plot(nom, xfrac = 0.3, cex.axis = 0.8, cex.var = 1.0,
       col.grid = gray(c(0.8, 0.95)), main = "TCGA-PAAD Prognostic Nomogram")
  grDevices::dev.off()
  cat("Saved:", out_png, "\n"); nom
}

# 11d) C-index with bootstrap
compute_cindex <- function(fit, nom_df, n_boot = 1000,
                           out_csv = "Nomogram_Cindex_PAAD.csv") {
  c_apparent <- as.numeric(fit$stats["C"])
  set.seed(42); v <- rms::validate(fit, method = "boot", B = n_boot)
  dxy_row <- v["Dxy", ]
  c_corrected <- 0.5 + dxy_row["index.corrected"] / 2
  boot_c <- numeric(n_boot); set.seed(42)
  for (i in seq_len(n_boot)) {
    idx <- sample(nrow(nom_df), replace = TRUE); boot_dat <- nom_df[idx, ]
    boot_c[i] <- tryCatch({
      dd_b <- rms::datadist(boot_dat); options(datadist = "dd_b")
      fit_b <- rms::cph(fit$sformula, data = boot_dat, x = TRUE, y = TRUE, surv = TRUE)
      as.numeric(fit_b$stats["C"])
    }, error = function(e) NA_real_)
  }
  options(datadist = "dd")
  boot_c <- boot_c[!is.na(boot_c)]
  ci_low <- as.numeric(quantile(boot_c, 0.025)); ci_high <- as.numeric(quantile(boot_c, 0.975))
  cindex_df <- data.frame(
    Metric = c("C-index (apparent)", "C-index (optimism-corrected)",
               "Bootstrap 95% CI lower", "Bootstrap 95% CI upper",
               "Dxy (original)", "Dxy (corrected)", "Optimism"),
    Value = c(c_apparent, c_corrected, ci_low, ci_high,
              dxy_row["index.orig"], dxy_row["index.corrected"], dxy_row["optimism"]),
    stringsAsFactors = FALSE)
  readr::write_csv(cindex_df, out_csv)
  cat("\n=== C-index Results ===\n"); print(cindex_df); cat("Saved:", out_csv, "\n")
  cindex_df
}

# 11e) Calibration curves
plot_calibration <- function(fit, nom_df, time_points = c(12, 24, 36),
                             n_boot = 200, out_png = "Calibration_PAAD.png", dpi = 600) {
  n_tp <- length(time_points)
  grDevices::png(out_png, width = 5 * n_tp, height = 5.5, units = "in", res = dpi)
  par(mfrow = c(1, n_tp), mar = c(5, 5, 3, 1))
  for (tp in time_points) {
    dd_cal <- rms::datadist(nom_df); options(datadist = "dd_cal")
    assign("dd_cal", dd_cal, envir = .GlobalEnv)
    fit_cal <- rms::cph(fit$sformula, data = nom_df, x = TRUE, y = TRUE, surv = TRUE, time.inc = tp)
    set.seed(42)
    cal <- tryCatch(rms::calibrate(fit_cal, u = tp, B = n_boot, cmethod = "KM"),
      error = function(e) { cat("Calibration failed for", tp, "months\n"); NULL })
    if (!is.null(cal)) {
      plot(cal, xlab = "Predicted Probability", ylab = "Observed Probability (KM)",
           main = paste0(round(tp/12, 1), "-Year OS Calibration"), subtitles = TRUE, cex.subtitle = 0.7)
      abline(0, 1, lty = 2, col = "red", lwd = 1.5)
    } else { plot.new(); title(main = paste0(round(tp/12, 1), "-Year (failed)")) }
  }
  grDevices::dev.off(); options(datadist = "dd"); cat("Saved:", out_png, "\n")
}

# 11f) Time-dependent ROC
plot_time_roc <- function(fit, nom_df, time_points = c(12, 24, 36),
                          out_png = "TimeROC_PAAD.png",
                          out_csv = "TimeROC_AUC_PAAD.csv", dpi = 600) {
  lp <- predict(fit, type = "lp")
  troc <- timeROC::timeROC(T = nom_df$OS_months, delta = nom_df$OS_event,
    marker = lp, cause = 1, weighting = "marginal", times = time_points, iid = TRUE)
  auc_vals <- troc$AUC; ci_list <- timeROC::confint(troc, level = 0.95)
  auc_df <- data.frame(
    Time_months = time_points, Time_years = time_points / 12,
    AUC = as.numeric(auc_vals),
    CI_low = as.numeric(ci_list$CI_AUC[, 1]),
    CI_high = as.numeric(ci_list$CI_AUC[, 2]), stringsAsFactors = FALSE)
  readr::write_csv(auc_df, out_csv)
  cat("\n=== Time-dependent AUC ===\n"); print(auc_df); cat("Saved:", out_csv, "\n")
  colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
  grDevices::png(out_png, width = 7, height = 7, units = "in", res = dpi)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "1 - Specificity", ylab = "Sensitivity",
       main = "TCGA-PAAD: Time-Dependent ROC Curves")
  abline(0, 1, lty = 2, col = "grey50")
  for (i in seq_along(time_points)) {
    tp_idx <- which(troc$times == time_points[i])
    lines(troc$FP[, tp_idx], troc$TP[, tp_idx], col = colors[i], lwd = 2.5)
  }
  legend("bottomright",
    legend = paste0(time_points/12, "-yr AUC = ", sprintf("%.3f", auc_df$AUC),
      " (", sprintf("%.3f", auc_df$CI_low), "\u2013", sprintf("%.3f", auc_df$CI_high), ")"),
    col = colors[seq_along(time_points)], lwd = 2.5, cex = 0.9, bty = "n")
  grDevices::dev.off(); cat("Saved:", out_png, "\n"); auc_df
}

# 11g) Master nomogram runner
run_nomogram_system <- function(genes = GENES_TO_TEST, time_points = c(12, 24, 36),
                                n_boot_cindex = 1000, n_boot_calib = 200,
                                include_clinical = TRUE, include_tme = TRUE,
                                prefix = "PAAD") {
  cat("\n", strrep("=", 60), "\n  NOMOGRAM SYSTEM: TCGA-PAAD\n", strrep("=", 60), "\n\n")
  cat("--- Step 1: Building nomogram data ---\n")
  nom_df <- build_nomogram_data(genes)
  cat("\n--- Step 2: Fitting Cox model ---\n")
  fit <- fit_nomogram_model(nom_df, genes, include_clinical, include_tme)
  cat("\n--- Model Summary ---\n"); print(fit)
  cat("\n--- Step 3: Nomogram ---\n")
  nom <- plot_nomogram(fit, time_points, out_png = paste0("Nomogram_", prefix, ".png"))
  cat("\n--- Step 4: C-index ---\n")
  cindex <- compute_cindex(fit, nom_df, n_boot_cindex, paste0("Nomogram_Cindex_", prefix, ".csv"))
  cat("\n--- Step 5: Calibration ---\n")
  plot_calibration(fit, nom_df, time_points, n_boot_calib, paste0("Calibration_", prefix, ".png"))
  cat("\n--- Step 6: Time-dependent ROC ---\n")
  auc <- plot_time_roc(fit, nom_df, time_points,
    paste0("TimeROC_", prefix, ".png"), paste0("TimeROC_AUC_", prefix, ".csv"))
  cat("\n", strrep("=", 60), "\n  NOMOGRAM SYSTEM COMPLETE\n", strrep("=", 60), "\n")
  invisible(list(nom_df = nom_df, fit = fit, nomogram = nom, cindex = cindex, auc = auc))
}

# ==========================================================
# 12) Run all analyses
# ==========================================================
gene_list <- GENES_TO_TEST

# Comparison table
tme_table <- build_gene_TME_compare_table(genes = gene_list, q = Q_CUTS, out_csv = OUT_CSV)
print(tme_table)

# Forest plot
plot_forest(tme_table, out_png = "Forest_PAAD_unadj_vs_TMEadj.png")

# KM plots
for (g in gene_list) plot_gene_km(gene = g, q = Q_CUTS)

# Stratified Cox: estimate_score
strata_est <- run_stratified_cox(genes = gene_list, q = Q_CUTS,
  tme_var = "estimate_score", tme_split = "median",
  out_csv = "PAAD_stratified_cox_estimate.csv")

# Stratified Cox: stromal_score
strata_str <- run_stratified_cox(genes = gene_list, q = Q_CUTS,
  tme_var = "stromal_score", tme_split = "median",
  out_csv = "PAAD_stratified_cox_stromal.csv")

# Nomogram system
nom_results <- run_nomogram_system(genes = gene_list,
  time_points = c(12, 24, 36), n_boot_cindex = 1000, n_boot_calib = 200,
  include_clinical = TRUE, include_tme = TRUE, prefix = "PAAD")

cat("\n=== All done! ===\n")
cat("Data source:", DATA_SOURCE, "-", .source_label, "\n")
if (DATA_SOURCE == "TCGAbiolinks") {
  cat("Validated: matches UCSC Xena star_tpm exactly\n")
} else {
  cat("Note: curatedTCGAData uses RSEM upper-quartile normalization (differs from STAR-Counts TPM)\n")
}
cat("Table:          ", OUT_CSV, "\n")
cat("Forest plot:     Forest_PAAD_unadj_vs_TMEadj.png\n")
cat("KM plots:        KM_<gene>_PAAD.png\n")
cat("Strata Cox:      PAAD_stratified_cox_estimate.csv\n")
cat("Strata Cox:      PAAD_stratified_cox_stromal.csv\n")
cat("Nomogram:        Nomogram_PAAD.png\n")
cat("C-index:         Nomogram_Cindex_PAAD.csv\n")
cat("Calibration:     Calibration_PAAD.png\n")
cat("Time-dep ROC:    TimeROC_PAAD.png\n")
cat("Time-dep AUC:    TimeROC_AUC_PAAD.csv\n")
