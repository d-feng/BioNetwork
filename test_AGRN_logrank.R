############################################################
# Self-contained AGRN Validation — Q1 vs Q4 (TCGA-PAAD)
#
# Reference (KMplot online tool):
#   RNAseq ID:       AGRN
#   Survival:        OS
#   Trichotomization: Q1 vs Q4
#   Cutoff value:    10083.5
#   Expression range: 700 – 48534
#   HR = 1.83 (0.98–3.4), logrank P = 0.054
#   N = 44 / 44
#   Follow-up threshold: all
#   Censored at threshold: checked
#   Compute median over entire database: false
#
# Strategy:
#   1. Download TCGA-PAAD STAR-Counts fresh
#   2. Extract AGRN from EVERY available assay
#      (unstranded, tpm_unstrand, fpkm_unstrand, fpkm_uq_unstrand)
#   3. Test each assay with RAW (no log) + log2(TPM+1) scales
#   4. Print side-by-side comparison vs reference
#   5. Generate KM plots for every assay
############################################################

# ==========================================================
# 0) Packages
# ==========================================================
suppressPackageStartupMessages({
  pkgs <- c("TCGAbiolinks", "SummarizedExperiment", "dplyr",
            "tibble", "stringr", "survival", "survminer", "ggplot2",
            "data.table")
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) BiocManager::install("TCGAbiolinks")
  lapply(pkgs, library, character.only = TRUE)
})

# ==========================================================
# 1) Download TCGA-PAAD STAR-Counts
# ==========================================================
PROJECT <- "TCGA-PAAD"

query <- GDCquery(
  project       = PROJECT,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type  = "STAR - Counts"
)
GDCdownload(query)
se <- GDCprepare(query)

# ==========================================================
# 2) Clinical / OS endpoint
# ==========================================================
clin <- GDCquery_clinic(project = PROJECT, type = "clinical")

clin2 <- clin %>%
  dplyr::mutate(
    patient   = dplyr::coalesce(submitter_id, bcr_patient_barcode),
    age_years = age_at_diagnosis / 365.25
  ) %>%
  dplyr::select(
    patient, vital_status, days_to_death, days_to_last_follow_up,
    age_years, gender, ajcc_pathologic_stage
  ) %>%
  dplyr::distinct(patient, .keep_all = TRUE)

# ==========================================================
# 3) Shared helpers
# ==========================================================

# Primary tumor barcodes
barcodes_tp <- TCGAquery_SampleTypes(colnames(se), typesample = "TP")
cat("Primary tumor samples:", length(barcodes_tp), "\n")

# Gene symbols from rowData
rowdat <- as.data.frame(SummarizedExperiment::rowData(se))
gene_symbol <- rowdat$gene_name

# List all available assays
all_assays <- SummarizedExperiment::assayNames(se)
cat("Available assays:", paste(all_assays, collapse = ", "), "\n\n")

# Extract one assay → matrix (genes × samples), collapse dup symbols
extract_assay_mat <- function(se, assay_name, barcodes, symbols) {
  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    cat("  Assay '", assay_name, "' not found — skipping\n")
    return(NULL)
  }
  mat <- SummarizedExperiment::assay(se, assay_name)
  mat <- mat[, barcodes, drop = FALSE]
  rownames(mat) <- symbols
  # Collapse duplicates by mean
  df <- as.data.frame(mat) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(gene) & gene != "") %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean), .groups = "drop") %>%
    tibble::column_to_rownames("gene")
  as.matrix(df)
}

# Build surv_df for a given expression matrix
build_surv <- function(mat, clin2) {
  patient_id <- substr(colnames(mat), 1, 12)
  sdf <- tibble::tibble(sample = colnames(mat), patient = patient_id) %>%
    dplyr::left_join(clin2, by = "patient") %>%
    dplyr::mutate(
      OS_time   = dplyr::if_else(!is.na(days_to_death),
                                 days_to_death, days_to_last_follow_up),
      OS_event  = dplyr::if_else(tolower(vital_status) == "dead", 1L, 0L),
      OS_months = OS_time / 30.44
    ) %>%
    dplyr::filter(!is.na(OS_time) & OS_time > 0)
  sdf
}

# ==========================================================
# 4) AGRN quick-look across all assays (raw scale)
# ==========================================================
gene <- "AGRN"
cat("==== AGRN expression range across assays (raw, no log) ====\n")

for (a in all_assays) {
  mat <- extract_assay_mat(se, a, barcodes_tp, gene_symbol)
  if (!is.null(mat) && gene %in% rownames(mat)) {
    vals <- as.numeric(mat[gene, ])
    cat(sprintf("  %-22s  min=%10.1f  max=%10.1f  median=%10.1f  Q25=%10.1f  Q75=%10.1f\n",
                a,
                min(vals, na.rm = TRUE),
                max(vals, na.rm = TRUE),
                median(vals, na.rm = TRUE),
                quantile(vals, 0.25, na.rm = TRUE),
                quantile(vals, 0.75, na.rm = TRUE)))
  }
}
cat("  REFERENCE              min=     700.0  max=   48534.0  cutoff=  10083.5\n\n")

# ==========================================================
# 5) Run full survival test per assay × scale
#    Scales: "raw" = no transform, "log2p1" = log2(x+1)
# ==========================================================

run_agrn_survival <- function(mat, sdf, assay_name, scale_name, gene) {

  # Gene expression
  sdf$gene_expr <- as.numeric(mat[gene, sdf$sample])
  sdf <- sdf %>% dplyr::filter(!is.na(gene_expr))

  expr_min <- min(sdf$gene_expr, na.rm = TRUE)
  expr_max <- max(sdf$gene_expr, na.rm = TRUE)

  # Q1 vs Q4
  qs     <- stats::quantile(sdf$gene_expr, probs = c(0.25, 0.75), na.rm = TRUE)
  cut_Q1 <- as.numeric(qs[1])
  cut_Q4 <- as.numeric(qs[2])

  sdf <- sdf %>%
    dplyr::mutate(group = dplyr::case_when(
      gene_expr <= cut_Q1 ~ "Low",
      gene_expr >= cut_Q4 ~ "High",
      TRUE ~ NA_character_
    )) %>%
    dplyr::filter(!is.na(group))

  sdf$group <- factor(sdf$group, levels = c("Low", "High"))

  n_lo <- sum(sdf$group == "Low")
  n_hi <- sum(sdf$group == "High")

  # Cox HR
  cox_fit <- survival::coxph(survival::Surv(OS_months, OS_event) ~ group, data = sdf)
  cox_sum <- summary(cox_fit)
  hr   <- cox_sum$conf.int["groupHigh", "exp(coef)"]
  ci_l <- cox_sum$conf.int["groupHigh", "lower .95"]
  ci_u <- cox_sum$conf.int["groupHigh", "upper .95"]

  # Log-rank
  sdiff <- survival::survdiff(survival::Surv(OS_months, OS_event) ~ group, data = sdf)
  logrank_p <- 1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1)

  # KM plot
  fit <- survival::survfit(survival::Surv(OS_months, OS_event) ~ group, data = sdf)
  hr_label <- sprintf("HR = %.2f (%.2f\u2013%.2f)\nLog-rank p = %s",
                       hr, ci_l, ci_u, signif(logrank_p, 3))

  tag <- paste0(assay_name, "_", scale_name)

  g <- survminer::ggsurvplot(
    fit, data = sdf,
    palette    = c("black", "grey50"),
    linetype   = c("solid", "dashed"),
    conf.int   = FALSE,
    pval       = FALSE,
    risk.table = TRUE,
    risk.table.col = "black",
    ggtheme    = ggplot2::theme_classic(base_size = 16),
    title      = paste0(gene, " \u2014 Q1 vs Q4 [", tag, "]  (TCGA-PAAD)"),
    legend.title = "",
    legend.labs  = c("Low (Q1)", "High (Q4)"),
    xlab = "Time (months)",
    ylab = "Overall survival probability"
  )

  g$plot <- g$plot +
    ggplot2::annotate("text",
      x = max(sdf$OS_months, na.rm = TRUE) * 0.45, y = 0.95,
      label = hr_label, size = 4.5, hjust = 0
    ) +
    ggplot2::theme(
      panel.border    = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.0),
      axis.line       = ggplot2::element_line(linewidth = 0.9),
      axis.ticks      = ggplot2::element_line(linewidth = 0.9),
      legend.position = "right"
    )

  out_file <- paste0("KM_AGRN_Q1Q4_", tag, ".png")
  grDevices::png(out_file, width = 8, height = 7, units = "in", res = 600)
  print(g)
  grDevices::dev.off()

  list(
    assay     = assay_name,
    scale     = scale_name,
    tag       = tag,
    hr        = hr,
    ci_l      = ci_l,
    ci_u      = ci_u,
    p         = logrank_p,
    n_lo      = n_lo,
    n_hi      = n_hi,
    range_min = expr_min,
    range_max = expr_max,
    cut_Q1    = cut_Q1,
    cut_Q4    = cut_Q4,
    file      = out_file
  )
}

# ---------- Loop over assays and scales ----------
results <- list()

for (a in all_assays) {
  mat_raw <- extract_assay_mat(se, a, barcodes_tp, gene_symbol)
  if (is.null(mat_raw) || !gene %in% rownames(mat_raw)) next

  sdf <- build_surv(mat_raw, clin2)

  # A) RAW scale
  cat(sprintf("\n--- Running: %s / raw ---\n", a))
  res_raw <- tryCatch(
    run_agrn_survival(mat_raw, sdf, a, "raw", gene),
    error = function(e) { message("  Error: ", e$message); NULL }
  )
  if (!is.null(res_raw)) results[[length(results) + 1]] <- res_raw

  # B) log2(x + 1) scale
  mat_log <- log2(mat_raw + 1)
  cat(sprintf("--- Running: %s / log2p1 ---\n", a))
  res_log <- tryCatch(
    run_agrn_survival(mat_log, sdf, a, "log2p1", gene),
    error = function(e) { message("  Error: ", e$message); NULL }
  )
  if (!is.null(res_log)) results[[length(results) + 1]] <- res_log
}

# ==========================================================
# 6) Summary comparison table
# ==========================================================
cat("\n\n")
cat("==========================================================================\n")
cat("  AGRN Q1 vs Q4  —  Which assay/scale matches the reference?\n")
cat("==========================================================================\n")
cat(sprintf("%-30s %6s %14s %8s %5s/%5s %12s %12s %10s\n",
            "Assay_Scale", "HR", "95% CI", "p", "N_lo", "N_hi",
            "Range_min", "Range_max", "Q75_cut"))
cat(paste(rep("-", 115), collapse = ""), "\n")

for (r in results) {
  cat(sprintf("%-30s %6.2f (%5.2f\u2013%5.2f) %8.4f %5d/%5d %12.1f %12.1f %10.1f\n",
              r$tag, r$hr, r$ci_l, r$ci_u, r$p,
              r$n_lo, r$n_hi, r$range_min, r$range_max, r$cut_Q4))
}
cat(paste(rep("-", 115), collapse = ""), "\n")
cat(sprintf("%-30s %6.2f (%5.2f\u2013%5.1f) %8.4f %5d/%5d %12.1f %12.1f %10.1f\n",
            "REFERENCE", 1.83, 0.98, 3.4, 0.054,
            44, 44, 700.0, 48534.0, 10083.5))

cat("\n")
cat("Look for the row whose Range_min/max and Q75_cut best match the REFERENCE.\n")
cat("That assay+scale is what the KMplot tool uses.\n")
cat("KM plots saved for each combination.\n\n")

# ==========================================================
# 7) Highlight best match
# ==========================================================
if (length(results) > 0) {
  # Score each result: distance from reference range + cutoff
  ref_min <- 700; ref_max <- 48534; ref_cut <- 10083.5
  scores <- sapply(results, function(r) {
    abs(log(r$range_min / ref_min)) +
      abs(log(r$range_max / ref_max)) +
      abs(log(r$cut_Q4 / ref_cut))
  })
  best_idx <- which.min(scores)
  best <- results[[best_idx]]
  cat(sprintf(">>> BEST MATCH:  %s\n", best$tag))
  cat(sprintf("    HR = %.2f (%.2f\u2013%.2f),  p = %.4f,  N = %d/%d\n",
              best$hr, best$ci_l, best$ci_u, best$p, best$n_lo, best$n_hi))
  cat(sprintf("    Range: %.1f – %.1f,  Q75 cutoff: %.1f\n",
              best$range_min, best$range_max, best$cut_Q4))
  cat(sprintf("    KM plot: %s\n", best$file))
}

# ==========================================================
# 8) CONTINUOUS Cox model (TIMER2 validation)
#    Model: Surv(OS, EVENT) ~ AGRN  (no dichotomization)
#
#    TIMER2 reference (n=179, 93 events):
#      coef=0.254  HR=1.289  se=0.1  95%CI 1.06–1.567
#      Wald p=0.011   LRT p=0.0082   Score p=0.0104
# ==========================================================

cat("\n\n")
cat("##########################################################\n")
cat("  CONTINUOUS Cox model — Surv(OS, EVENT) ~ AGRN\n")
cat("  (Matches TIMER2 approach: no dichotomization)\n")
cat("##########################################################\n\n")

cont_results <- list()

for (a in all_assays) {
  mat_raw <- extract_assay_mat(se, a, barcodes_tp, gene_symbol)
  if (is.null(mat_raw) || !gene %in% rownames(mat_raw)) next

  sdf <- build_surv(mat_raw, clin2)

  # --- RAW scale ---
  sdf_raw <- sdf
  sdf_raw$gene_expr <- as.numeric(mat_raw[gene, sdf_raw$sample])
  sdf_raw <- sdf_raw %>% dplyr::filter(!is.na(gene_expr))

  cox_raw <- survival::coxph(
    survival::Surv(OS_time, OS_event) ~ gene_expr, data = sdf_raw
  )
  s_raw <- summary(cox_raw)

  cont_results[[paste0(a, "_raw")]] <- list(
    tag      = paste0(a, "_raw"),
    n        = s_raw$n,
    events   = s_raw$nevent,
    coef     = s_raw$coefficients["gene_expr", "coef"],
    hr       = s_raw$coefficients["gene_expr", "exp(coef)"],
    se       = s_raw$coefficients["gene_expr", "se(coef)"],
    ci_l     = s_raw$conf.int["gene_expr", "lower .95"],
    ci_u     = s_raw$conf.int["gene_expr", "upper .95"],
    p_wald   = s_raw$coefficients["gene_expr", "Pr(>|z|)"],
    p_lrt    = as.numeric(s_raw$logtest["pvalue"]),
    p_score  = as.numeric(s_raw$sctest["pvalue"]),
    rsq      = s_raw$rsq[1],
    rsq_max  = s_raw$rsq[2]
  )

  # --- log2(x+1) scale ---
  mat_log <- log2(mat_raw + 1)
  sdf_log <- sdf
  sdf_log$gene_expr <- as.numeric(mat_log[gene, sdf_log$sample])
  sdf_log <- sdf_log %>% dplyr::filter(!is.na(gene_expr))

  cox_log <- survival::coxph(
    survival::Surv(OS_time, OS_event) ~ gene_expr, data = sdf_log
  )
  s_log <- summary(cox_log)

  cont_results[[paste0(a, "_log2p1")]] <- list(
    tag      = paste0(a, "_log2p1"),
    n        = s_log$n,
    events   = s_log$nevent,
    coef     = s_log$coefficients["gene_expr", "coef"],
    hr       = s_log$coefficients["gene_expr", "exp(coef)"],
    se       = s_log$coefficients["gene_expr", "se(coef)"],
    ci_l     = s_log$conf.int["gene_expr", "lower .95"],
    ci_u     = s_log$conf.int["gene_expr", "upper .95"],
    p_wald   = s_log$coefficients["gene_expr", "Pr(>|z|)"],
    p_lrt    = as.numeric(s_log$logtest["pvalue"]),
    p_score  = as.numeric(s_log$sctest["pvalue"]),
    rsq      = s_log$rsq[1],
    rsq_max  = s_log$rsq[2]
  )
}

# --- Print continuous model summary table ---
cat(sprintf("%-26s %5s %4s %7s %6s %6s %10s %10s %8s %8s %8s\n",
            "Assay_Scale", "N", "Evt", "coef", "HR", "se", "95%CI_l",
            "95%CI_u", "p_wald", "p_lrt", "p_score"))
cat(paste(rep("-", 125), collapse = ""), "\n")

for (r in cont_results) {
  cat(sprintf("%-26s %5d %4d %7.4f %6.3f %6.4f %10.4f %10.4f %8.1e %8.1e %8.1e\n",
              r$tag, r$n, r$events, r$coef, r$hr, r$se,
              r$ci_l, r$ci_u, r$p_wald, r$p_lrt, r$p_score))
}
cat(paste(rep("-", 125), collapse = ""), "\n")
cat(sprintf("%-26s %5d %4d %7.4f %6.3f %6.4f %10.4f %10.4f %8.1e %8.1e %8.1e\n",
            "TIMER2 REFERENCE", 178, 93, 0.254, 1.289, 0.100,
            1.060, 1.567, 1.09e-02, 8.2e-03, 1.04e-02))

cat("\nMatch criteria: coef~0.254, HR~1.289, se~0.1, N~178, events~93\n")

# --- Highlight best match for continuous model ---
if (length(cont_results) > 0) {
  ref_coef <- 0.254; ref_se <- 0.1; ref_hr <- 1.289
  cscores <- sapply(cont_results, function(r) {
    abs(r$coef - ref_coef) + abs(r$se - ref_se) + abs(r$hr - ref_hr)
  })
  best_c <- cont_results[[which.min(cscores)]]

  cat(sprintf("\n>>> BEST CONTINUOUS MATCH:  %s\n", best_c$tag))
  cat(sprintf("    N = %d,  Events = %d\n", best_c$n, best_c$events))
  cat(sprintf("    coef = %.4f,  HR = %.3f,  se = %.4f\n", best_c$coef, best_c$hr, best_c$se))
  cat(sprintf("    95%%CI: %.4f \u2013 %.4f\n", best_c$ci_l, best_c$ci_u))
  cat(sprintf("    Wald p = %.2e,  LRT p = %.2e,  Score p = %.2e\n",
              best_c$p_wald, best_c$p_lrt, best_c$p_score))
  cat(sprintf("    R-squared = %.4f (max = %.3f)\n", best_c$rsq, best_c$rsq_max))
}

# ==========================================================
# 9) UCSC Xena log2(TPM+1) vs GDC tpm_unstrand log2(TPM+1)
#    Source: TCGA-PAAD.star_tpm.tsv (Ensembl IDs, 16-char barcodes)
#    AGRN Ensembl ID: ENSG00000188157
# ==========================================================

cat("\n\n")
cat("##########################################################\n")
cat("  UCSC Xena log2(TPM+1) vs GDC tpm_unstrand log2(TPM+1)\n")
cat("##########################################################\n\n")

UCSC_FILE <- "C:/Users/difen/Rcode/TCGA-PAAD.star_tpm.tsv"
if (!file.exists(UCSC_FILE)) {
  UCSC_FILE <- "TCGA-PAAD.star_tpm.tsv"  # try working directory
}

if (!file.exists(UCSC_FILE)) {
  cat("WARNING: UCSC file not found. Skipping Section 9.\n")
} else {

  cat("Reading UCSC file:", UCSC_FILE, "\n")
  ucsc_raw <- data.table::fread(UCSC_FILE, header = TRUE, sep = "\t",
                                 check.names = FALSE, data.table = FALSE)

  # Row names = Ensembl IDs (strip version for matching)
  ucsc_ensembl <- ucsc_raw[[1]]
  ucsc_ensembl_base <- sub("\\.\\d+$", "", ucsc_ensembl)
  ucsc_mat <- as.matrix(ucsc_raw[, -1])
  rownames(ucsc_mat) <- ucsc_ensembl

  cat("UCSC matrix:", nrow(ucsc_mat), "genes x", ncol(ucsc_mat), "samples\n")

  # --- Map AGRN Ensembl ID ---
  agrn_ens <- "ENSG00000188157"
  agrn_rows <- grep(agrn_ens, rownames(ucsc_mat))
  if (length(agrn_rows) == 0) {
    cat("ERROR: AGRN (", agrn_ens, ") not found in UCSC file\n")
  } else {
    agrn_ens_full <- rownames(ucsc_mat)[agrn_rows[1]]
    cat("AGRN Ensembl ID in UCSC:", agrn_ens_full, "\n")

    ucsc_agrn <- as.numeric(ucsc_mat[agrn_rows[1], ])
    names(ucsc_agrn) <- colnames(ucsc_mat)

    cat(sprintf("UCSC AGRN log2(TPM+1):  min=%.3f  max=%.3f  median=%.3f\n",
                min(ucsc_agrn, na.rm=TRUE), max(ucsc_agrn, na.rm=TRUE),
                median(ucsc_agrn, na.rm=TRUE)))

    # --- GDC tpm_unstrand log2(TPM+1) for AGRN ---
    mat_gdc_tpm <- extract_assay_mat(se, "tpm_unstrand", barcodes_tp, gene_symbol)
    mat_gdc_log2 <- log2(mat_gdc_tpm + 1)
    gdc_agrn <- as.numeric(mat_gdc_log2["AGRN", ])
    names(gdc_agrn) <- colnames(mat_gdc_log2)

    cat(sprintf("GDC  AGRN log2(TPM+1):  min=%.3f  max=%.3f  median=%.3f\n\n",
                min(gdc_agrn, na.rm=TRUE), max(gdc_agrn, na.rm=TRUE),
                median(gdc_agrn, na.rm=TRUE)))

    # --- Match samples by patient ID (first 12 chars) ---
    # UCSC barcodes: TCGA-XX-XXXX-01A (16 chars)
    # GDC  barcodes: TCGA-XX-XXXX-01A-11R-A24H-07 (28 chars)
    ucsc_patients <- substr(names(ucsc_agrn), 1, 15)
    gdc_patients  <- substr(names(gdc_agrn), 1, 15)

    # Find shared samples (match on 15-char barcode)
    shared_keys <- intersect(ucsc_patients, gdc_patients)
    cat("UCSC samples:", length(ucsc_agrn), "\n")
    cat("GDC  samples:", length(gdc_agrn), "\n")
    cat("Shared (15-char match):", length(shared_keys), "\n\n")

    if (length(shared_keys) > 0) {
      # Build matched vectors
      ucsc_idx <- match(shared_keys, ucsc_patients)
      gdc_idx  <- match(shared_keys, gdc_patients)

      ucsc_vals <- ucsc_agrn[ucsc_idx]
      gdc_vals  <- gdc_agrn[gdc_idx]

      # Correlation
      cor_val <- cor(ucsc_vals, gdc_vals, use = "complete.obs")
      diff_vals <- ucsc_vals - gdc_vals
      cat(sprintf("Pearson correlation:    r = %.6f\n", cor_val))
      cat(sprintf("Mean difference:        %.6f\n", mean(diff_vals, na.rm=TRUE)))
      cat(sprintf("Max  abs difference:    %.6f\n", max(abs(diff_vals), na.rm=TRUE)))
      cat(sprintf("SD   difference:        %.6f\n", sd(diff_vals, na.rm=TRUE)))

      # Show a few matched values
      cat("\nSample matched values (first 10):\n")
      cat(sprintf("%-18s  %12s  %12s  %10s\n", "Barcode_15", "UCSC", "GDC", "Diff"))
      cat(paste(rep("-", 58), collapse = ""), "\n")
      show_n <- min(10, length(shared_keys))
      for (i in seq_len(show_n)) {
        cat(sprintf("%-18s  %12.4f  %12.4f  %10.4f\n",
                    shared_keys[i], ucsc_vals[i], gdc_vals[i], diff_vals[i]))
      }

      # --- Scatter plot: UCSC vs GDC ---
      scatter_df <- data.frame(UCSC = ucsc_vals, GDC = gdc_vals)
      p_scatter <- ggplot2::ggplot(scatter_df, ggplot2::aes(x = GDC, y = UCSC)) +
        ggplot2::geom_point(alpha = 0.6, size = 2) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        ggplot2::labs(
          title = paste0("AGRN log2(TPM+1): UCSC Xena vs GDC  (r=", round(cor_val, 4), ")"),
          x = "GDC tpm_unstrand log2(TPM+1)",
          y = "UCSC Xena STAR-TPM log2(TPM+1)"
        ) +
        ggplot2::annotate("text", x = min(gdc_vals)*1.05, y = max(ucsc_vals)*0.95,
                          label = sprintf("r = %.4f\nN = %d\nMean diff = %.4f",
                                          cor_val, length(shared_keys),
                                          mean(diff_vals, na.rm=TRUE)),
                          hjust = 0, size = 4) +
        ggplot2::theme_classic(base_size = 14)

      ggplot2::ggsave("AGRN_UCSC_vs_GDC_scatter.png", p_scatter,
                      width = 7, height = 6, dpi = 600)
      cat("\nSaved: AGRN_UCSC_vs_GDC_scatter.png\n")

      # --- Run survival with UCSC data (continuous + Q1 vs Q4) ---
      cat("\n--- UCSC AGRN Survival (continuous Cox) ---\n")

      # Keep only primary tumor (01A) from UCSC
      ucsc_tp_idx <- grep("-01A$", names(ucsc_agrn))
      ucsc_agrn_tp <- ucsc_agrn[ucsc_tp_idx]

      # Build surv_df for UCSC
      ucsc_patient <- substr(names(ucsc_agrn_tp), 1, 12)
      sdf_ucsc <- tibble::tibble(
        sample  = names(ucsc_agrn_tp),
        patient = ucsc_patient
      ) %>%
        dplyr::left_join(clin2, by = "patient") %>%
        dplyr::mutate(
          OS_time  = dplyr::if_else(!is.na(days_to_death),
                                     days_to_death, days_to_last_follow_up),
          OS_event = dplyr::if_else(tolower(vital_status) == "dead", 1L, 0L),
          OS_months = OS_time / 30.44,
          gene_expr = ucsc_agrn_tp[sample]
        ) %>%
        dplyr::filter(!is.na(OS_time) & OS_time > 0 & !is.na(gene_expr))

      cat("UCSC primary tumor with OS:", nrow(sdf_ucsc), "\n")
      cat("Events:", sum(sdf_ucsc$OS_event == 1), "\n")

      # Continuous Cox
      cox_ucsc <- survival::coxph(
        survival::Surv(OS_time, OS_event) ~ gene_expr, data = sdf_ucsc
      )
      s_ucsc <- summary(cox_ucsc)

      cat(sprintf("\nUCSC Continuous Cox:\n"))
      cat(sprintf("  N = %d,  Events = %d\n", s_ucsc$n, s_ucsc$nevent))
      cat(sprintf("  coef = %.4f,  HR = %.3f,  se = %.4f\n",
                  s_ucsc$coefficients["gene_expr","coef"],
                  s_ucsc$coefficients["gene_expr","exp(coef)"],
                  s_ucsc$coefficients["gene_expr","se(coef)"]))
      cat(sprintf("  95%%CI: %.4f \u2013 %.4f\n",
                  s_ucsc$conf.int["gene_expr","lower .95"],
                  s_ucsc$conf.int["gene_expr","upper .95"]))
      cat(sprintf("  Wald p = %.2e,  LRT p = %.2e,  Score p = %.2e\n",
                  s_ucsc$coefficients["gene_expr","Pr(>|z|)"],
                  as.numeric(s_ucsc$logtest["pvalue"]),
                  as.numeric(s_ucsc$sctest["pvalue"])))

      # Q1 vs Q4 with UCSC
      cat("\n--- UCSC AGRN Survival (Q1 vs Q4) ---\n")
      qs_ucsc <- stats::quantile(sdf_ucsc$gene_expr, probs = c(0.25, 0.75), na.rm = TRUE)
      cat(sprintf("  Q25 cutoff: %.4f,  Q75 cutoff: %.4f\n",
                  as.numeric(qs_ucsc[1]), as.numeric(qs_ucsc[2])))

      sdf_ucsc_q <- sdf_ucsc %>%
        dplyr::mutate(group = dplyr::case_when(
          gene_expr <= as.numeric(qs_ucsc[1]) ~ "Low",
          gene_expr >= as.numeric(qs_ucsc[2]) ~ "High",
          TRUE ~ NA_character_
        )) %>%
        dplyr::filter(!is.na(group))
      sdf_ucsc_q$group <- factor(sdf_ucsc_q$group, levels = c("Low", "High"))

      cat(sprintf("  Low (Q1): %d,  High (Q4): %d\n",
                  sum(sdf_ucsc_q$group=="Low"), sum(sdf_ucsc_q$group=="High")))

      cox_ucsc_q <- survival::coxph(
        survival::Surv(OS_months, OS_event) ~ group, data = sdf_ucsc_q
      )
      s_ucsc_q <- summary(cox_ucsc_q)

      sdiff_ucsc <- survival::survdiff(
        survival::Surv(OS_months, OS_event) ~ group, data = sdf_ucsc_q
      )
      lr_ucsc <- 1 - pchisq(sdiff_ucsc$chisq, df = length(sdiff_ucsc$n) - 1)

      cat(sprintf("  HR = %.2f (%.2f\u2013%.2f)\n",
                  s_ucsc_q$conf.int["groupHigh","exp(coef)"],
                  s_ucsc_q$conf.int["groupHigh","lower .95"],
                  s_ucsc_q$conf.int["groupHigh","upper .95"]))
      cat(sprintf("  Log-rank p = %s\n", signif(lr_ucsc, 4)))

      # KM plot for UCSC Q1 vs Q4
      fit_ucsc <- survival::survfit(
        survival::Surv(OS_months, OS_event) ~ group, data = sdf_ucsc_q
      )
      hr_ucsc_label <- sprintf(
        "HR = %.2f (%.2f\u2013%.2f)\nLog-rank p = %s",
        s_ucsc_q$conf.int["groupHigh","exp(coef)"],
        s_ucsc_q$conf.int["groupHigh","lower .95"],
        s_ucsc_q$conf.int["groupHigh","upper .95"],
        signif(lr_ucsc, 3)
      )

      g_ucsc <- survminer::ggsurvplot(
        fit_ucsc, data = sdf_ucsc_q,
        palette    = c("black", "grey50"),
        linetype   = c("solid", "dashed"),
        conf.int   = FALSE,
        pval       = FALSE,
        risk.table = TRUE,
        risk.table.col = "black",
        ggtheme    = ggplot2::theme_classic(base_size = 16),
        title      = "AGRN \u2014 Q1 vs Q4 [UCSC Xena log2(TPM+1)]  (TCGA-PAAD)",
        legend.title = "",
        legend.labs  = c("Low (Q1)", "High (Q4)"),
        xlab = "Time (months)",
        ylab = "Overall survival probability"
      )

      g_ucsc$plot <- g_ucsc$plot +
        ggplot2::annotate("text",
          x = max(sdf_ucsc_q$OS_months, na.rm = TRUE) * 0.45, y = 0.95,
          label = hr_ucsc_label, size = 4.5, hjust = 0
        ) +
        ggplot2::theme(
          panel.border    = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.0),
          axis.line       = ggplot2::element_line(linewidth = 0.9),
          axis.ticks      = ggplot2::element_line(linewidth = 0.9),
          legend.position = "right"
        )

      grDevices::png("KM_AGRN_Q1Q4_UCSC_log2p1.png", width = 8, height = 7,
                     units = "in", res = 600)
      print(g_ucsc)
      grDevices::dev.off()
      cat("Saved: KM_AGRN_Q1Q4_UCSC_log2p1.png\n")
    }
  }

  # ===========================================================
  # 9b) Final cross-source comparison table
  # ===========================================================
  cat("\n\n")
  cat("##########################################################\n")
  cat("  FINAL COMPARISON: GDC vs UCSC vs TIMER2 vs KMplot\n")
  cat("##########################################################\n\n")

  cat(sprintf("%-28s %6s %7s %6s %10s %8s %8s\n",
              "Source", "N", "Events", "HR", "95%CI", "p_wald", "p_lrt"))
  cat(paste(rep("-", 85), collapse = ""), "\n")

  # GDC continuous (best match from Section 8)
  if (exists("best_c")) {
    cat(sprintf("%-28s %6d %7d %6.3f %10s %8.1e %8.1e\n",
                paste0("GDC ", best_c$tag, " (cont)"),
                best_c$n, best_c$events, best_c$hr,
                sprintf("%.2f\u2013%.2f", best_c$ci_l, best_c$ci_u),
                best_c$p_wald, best_c$p_lrt))
  }

  # UCSC continuous
  if (exists("s_ucsc")) {
    cat(sprintf("%-28s %6d %7d %6.3f %10s %8.1e %8.1e\n",
                "UCSC Xena log2p1 (cont)",
                s_ucsc$n, s_ucsc$nevent,
                s_ucsc$coefficients["gene_expr","exp(coef)"],
                sprintf("%.2f\u2013%.2f",
                        s_ucsc$conf.int["gene_expr","lower .95"],
                        s_ucsc$conf.int["gene_expr","upper .95"]),
                s_ucsc$coefficients["gene_expr","Pr(>|z|)"],
                as.numeric(s_ucsc$logtest["pvalue"])))
  }

  # TIMER2 reference
  cat(sprintf("%-28s %6d %7d %6.3f %10s %8.1e %8.1e\n",
              "TIMER2 (cont) REFERENCE",
              178, 93, 1.289, "1.06\u20131.57", 1.09e-02, 8.2e-03))

  # KMplot Q1 vs Q4 reference
  cat(sprintf("%-28s %6s %7s %6.2f %10s %8s %8s\n",
              "KMplot Q1Q4 REFERENCE",
              "44/44", "—", 1.83, "0.98\u20133.4", "—", "p_lr=0.054"))

  cat(paste(rep("-", 85), collapse = ""), "\n")
}

cat("\nDone.\n")
