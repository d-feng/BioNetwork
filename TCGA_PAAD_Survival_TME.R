############################################################
# TCGA-PAAD Survival (Unadjusted vs TME-adjusted)
# - Downloads TCGA-PAAD STAR-Counts, extracts TPM assay
# - Builds OS endpoint
# - Merges PAAD_Estimate.csv (estimate_score, stromal_score)
# - Runs gene-by-gene survival (Quantile stratification)
# - Uses LRT p-values (robust for small samples, N~88)
# - Outputs side-by-side table: Unadjusted vs TME-adjusted HR
# - Kaplan-Meier survival curves saved as PNG
############################################################

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

  lapply(pkgs, library, character.only = TRUE)
})

# ---------------------------
# USER SETTINGS
# ---------------------------
PROJECT     <- "TCGA-PAAD"
TPM_ASSAY   <- "tpm_unstrand"
ESTIMATE_CSV <- "PAAD_Estimate.csv"
TME_COVARS  <- c("estimate_score", "stromal_score")
GENES_TO_TEST <- c("KRAS","EGFR","MAPK1","MAPK3","PIK3CA","AKT1")
STRAT_METHOD  <- "quantile"
Q_CUTS        <- c(0.25, 0.75)
OUT_CSV       <- "PAAD_gene_unadj_vs_TMEadj.csv"

# ---------------------------
# 1) Download expression (STAR - Counts) + clinical
# ---------------------------
query <- GDCquery(
  project       = PROJECT,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type  = "STAR - Counts"
)
GDCdownload(query)
se <- GDCprepare(query)

# ---------------------------
# 2) Extract TPM, keep tumor, map symbols, collapse duplicates, log2(TPM+1)
# ---------------------------
if (!TPM_ASSAY %in% SummarizedExperiment::assayNames(se)) {
  stop(
    sprintf("TPM_ASSAY='%s' not found. Available assays: %s",
            TPM_ASSAY, paste(SummarizedExperiment::assayNames(se), collapse = ", "))
  )
}

expr <- SummarizedExperiment::assay(se, TPM_ASSAY)

# Keep primary tumor only
barcodes_tp <- TCGAquery_SampleTypes(colnames(expr), typesample = "TP")
expr <- expr[, barcodes_tp, drop = FALSE]

# Map gene symbols
rowdat <- as.data.frame(SummarizedExperiment::rowData(se))
if (!"gene_name" %in% colnames(rowdat)) {
  stop("rowData(se)$gene_name not found. Check rowData(se) columns.")
}
gene_symbol <- rowdat$gene_name
rownames(expr) <- gene_symbol

# Collapse duplicate symbols by mean
expr_df <- as.data.frame(expr) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(!is.na(gene) & gene != "") %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean), .groups = "drop") %>%
  tibble::column_to_rownames("gene")
expr <- as.matrix(expr_df)

# Log2(TPM+1)
expr <- log2(expr + 1)

# ---------------------------
# 3) Clinical + OS endpoint
# ---------------------------
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
    stage_simple = factor(stage_simple, levels = c("I","II","III","IV"))
  ) %>%
  dplyr::filter(!is.na(OS_time) & OS_time > 0)

# Align expression to survival samples
expr <- expr[, surv_df$sample, drop = FALSE]
cat("PAAD tumor samples with OS:", nrow(surv_df), "\n")

# ---------------------------
# 4) Read ESTIMATE scores and merge
# ---------------------------
score_raw <- readr::read_csv(ESTIMATE_CSV, show_col_types = FALSE)
score <- score_raw
names(score) <- names(score) %>%
  stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
  stringr::str_replace_all("_+$", "") %>%
  tolower()

# Identify ID column
id_col <- intersect(c("id", "sample", "barcode", "tcga_id"), names(score))[1]
if (is.na(id_col)) {
  stop("Could not find an ID column in score file. Expected one of: id, sample, barcode, tcga_id")
}
score <- score %>% dplyr::rename(sample = !!id_col)

# Normalize ESTIMATE column names
rename_map <- c(
  "stromal_score"  = "stromalscore",
  "immune_score"   = "immunescore",
  "estimate_score" = "estimatescore",
  "tumor_purity"   = "tumorpurity"
)
for (nm in names(rename_map)) {
  if (rename_map[[nm]] %in% names(score) && !(nm %in% names(score))) {
    score <- score %>% dplyr::rename(!!nm := !!rename_map[[nm]])
  }
}

# Merge on first 15 characters of barcode
surv_df <- surv_df %>%
  dplyr::mutate(sample_key = substr(sample, 1, 15))

score2 <- score %>%
  dplyr::mutate(sample_key = substr(sample, 1, 15))

surv_df <- surv_df %>%
  dplyr::left_join(
    score2 %>% dplyr::select(sample_key, stromal_score, immune_score, estimate_score),
    by = "sample_key"
  ) %>%
  dplyr::select(-sample_key)

cat("NA estimate_score:", sum(is.na(surv_df$estimate_score)), "\n")
cat("NA stromal_score :", sum(is.na(surv_df$stromal_score)), "\n")

# ---------------------------
# 5) Helper: extract HR, CI, and three p-values (LRT, Wald, Log-rank)
#    LRT   = recommended for small samples (N~88), works with covariates
#    Wald  = standard Cox output, included for reference
#    Log-rank = non-parametric, unadjusted only (NA for adjusted models)
# ---------------------------
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

  HR <- exp(beta)
  L  <- exp(beta - 1.96 * se)
  U  <- exp(beta + 1.96 * se)

  # Wald p-value (from Cox summary)
  p_wald <- as.numeric(s$coefficients[idx, "Pr(>|z|)"])

  # LRT p-value for the group term (robust for small samples)
  p_lrt <- tryCatch({
    d1 <- drop1(fit, test = "Chisq")
    group_row <- grep("^group$", rownames(d1))
    if (length(group_row) > 0) {
      as.numeric(d1[group_row, "Pr(>Chi)"])
    } else {
      # Single-variable model: use overall LRT
      as.numeric(s$logtest["pvalue"])
    }
  }, error = function(e) {
    as.numeric(s$logtest["pvalue"])
  })

  # Log-rank p-value (non-parametric, unadjusted only)
  # Only valid when dat is provided and model has group as only covariate
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

# ---------------------------
# 6) Build comparison table: Unadjusted vs TME-adjusted
# ---------------------------
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

    # Cox models
    cox_unadj <- tryCatch(
      survival::coxph(survival::Surv(OS_time, OS_event) ~ group, data = dat),
      error = function(e) NULL
    )
    cox_est <- tryCatch(
      survival::coxph(survival::Surv(OS_time, OS_event) ~ group + estimate_score, data = dat),
      error = function(e) NULL
    )
    cox_str <- tryCatch(
      survival::coxph(survival::Surv(OS_time, OS_event) ~ group + stromal_score, data = dat),
      error = function(e) NULL
    )

    unadj <- extract_groupHigh_stats(cox_unadj, dat = dat)  # pass dat for log-rank
    est   <- extract_groupHigh_stats(cox_est)               # adjusted: log-rank = NA
    str   <- extract_groupHigh_stats(cox_str)               # adjusted: log-rank = NA

    results[[g]] <- data.frame(
      Gene             = g,
      Method           = "quantile",
      N_total          = n_total,
      Events           = n_events,
      Cutoff_low       = cut_low,
      Cutoff_high      = cut_high,

      HR_unadj              = unadj["HR"],
      CI_low_unadj          = unadj["L"],
      CI_high_unadj         = unadj["U"],
      p_lrt_unadj           = unadj["p_lrt"],
      p_wald_unadj          = unadj["p_wald"],
      p_logrank_unadj       = unadj["p_logrank"],

      HR_adj_estimate       = est["HR"],
      CI_low_adj_estimate   = est["L"],
      CI_high_adj_estimate  = est["U"],
      p_lrt_adj_estimate    = est["p_lrt"],
      p_wald_adj_estimate   = est["p_wald"],

      HR_adj_stromal        = str["HR"],
      CI_low_adj_stromal    = str["L"],
      CI_high_adj_stromal   = str["U"],
      p_lrt_adj_stromal     = str["p_lrt"],
      p_wald_adj_stromal    = str["p_wald"],

      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  }

  out <- dplyr::bind_rows(results)
  readr::write_csv(out, out_csv)
  cat("Saved:", out_csv, "\n")
  out
}

# ---------------------------
# 7) Run gene survival + TME-adjusted analysis
# ---------------------------
run_gene_survival_TME <- function(gene,
                                  q              = c(0.25, 0.75),
                                  covar_estimate = "estimate_score",
                                  covar_stromal  = "stromal_score") {

  if (!gene %in% rownames(expr)) stop("Gene not found in expr: ", gene)

  dat <- surv_df %>%
    dplyr::mutate(gene_expr = as.numeric(expr[gene, sample])) %>%
    dplyr::filter(!is.na(gene_expr), !is.na(OS_time), OS_time > 0)

  if (nrow(dat) < 10) stop("Too few samples after filtering for gene_expr/OS_time.")

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

  if (nrow(dat) < 10) stop("Too few samples after stratification.")
  if (length(unique(dat$group)) < 2) stop("Only one group after stratification.")
  if (sum(dat$OS_event == 1, na.rm = TRUE) < 1) stop("No events after stratification.")

  # Cox models
  cox_unadj <- survival::coxph(survival::Surv(OS_time, OS_event) ~ group, data = dat)

  ok_covar <- function(x) {
    v <- dat[[x]]
    if (is.null(v)) return(FALSE)
    if (all(is.na(v))) return(FALSE)
    length(unique(na.omit(v))) >= 2
  }

  cox_adj_estimate <- NULL
  if (ok_covar(covar_estimate)) {
    cox_adj_estimate <- survival::coxph(
      stats::as.formula(paste0("survival::Surv(OS_time, OS_event) ~ group + ", covar_estimate)),
      data = dat
    )
  }

  cox_adj_stromal <- NULL
  if (ok_covar(covar_stromal)) {
    cox_adj_stromal <- survival::coxph(
      stats::as.formula(paste0("survival::Surv(OS_time, OS_event) ~ group + ", covar_stromal)),
      data = dat
    )
  }

  km_fit <- survival::survfit(survival::Surv(OS_time, OS_event) ~ group, data = dat)

  list(
    gene             = gene,
    method           = "quantile",
    cutoff_low       = cutoff_low,
    cutoff_high      = cutoff_high,
    data             = dat,
    km_fit           = km_fit,
    cox_unadj        = cox_unadj,
    cox_adj_estimate = cox_adj_estimate,
    cox_adj_stromal  = cox_adj_stromal
  )
}

# ---------------------------
# 8) KM plot function (LRT p-value, Nature-style)
# ---------------------------
plot_gene_km <- function(gene,
                         q       = c(0.25, 0.75),
                         out_png = NULL,
                         dpi     = 600,
                         show_hr = TRUE) {

  res <- run_gene_survival_TME(gene, q = q)
  dat <- res$data

  title_txt <- paste0(gene, " (Q", q[1]*100, " vs Q", q[2]*100, ")")

  # Extract HR + all three p-values for subtitle
  extract_group_hr <- function(fit, fit_dat = NULL) {
    if (is.null(fit)) return(NULL)
    s   <- summary(fit)
    rn  <- rownames(s$coefficients)
    idx <- which(grepl("^groupHigh$", rn))
    if (!length(idx)) return(NULL)

    beta <- s$coefficients[idx, "coef"]
    se   <- s$coefficients[idx, "se(coef)"]
    hr   <- exp(beta)
    ci_l <- exp(beta - 1.96 * se)
    ci_u <- exp(beta + 1.96 * se)

    # Wald p-value
    p_wald <- as.numeric(s$coefficients[idx, "Pr(>|z|)"])

    # LRT p-value
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

    # Log-rank p-value (only for unadjusted / single-covariate models)
    p_logrank <- NA_real_
    if (!is.null(fit_dat) && nrow(s$coefficients) == 1) {
      p_logrank <- tryCatch({
        sdiff <- survival::survdiff(survival::Surv(OS_time, OS_event) ~ group, data = fit_dat)
        as.numeric(1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1))
      }, error = function(e) NA_real_)
    }

    list(HR = as.numeric(hr), L = as.numeric(ci_l), U = as.numeric(ci_u),
         p_lrt = as.numeric(p_lrt), p_wald = as.numeric(p_wald),
         p_logrank = as.numeric(p_logrank))
  }

  sub_txt <- NULL
  if (show_hr) {
    u <- extract_group_hr(res$cox_unadj, fit_dat = dat)  # pass dat for log-rank
    a <- extract_group_hr(res$cox_adj_estimate)
    b <- extract_group_hr(res$cox_adj_stromal)

    fmt <- function(x, label) {
      if (is.null(x)) return(paste0(label, ": NA"))
      p_txt <- paste0("LRT p=", signif(x$p_lrt, 3),
                       ", Wald p=", signif(x$p_wald, 3))
      if (!is.na(x$p_logrank)) {
        p_txt <- paste0(p_txt, ", Log-rank p=", signif(x$p_logrank, 3))
      }
      paste0(label, ": HR=", sprintf("%.2f", x$HR),
             " (", sprintf("%.2f", x$L), "\u2013", sprintf("%.2f", x$U),
             "), ", p_txt)
    }

    sub_txt <- paste(
      fmt(u, "Unadj"),
      fmt(a, "Adj+Estimate"),
      fmt(b, "Adj+Stroma"),
      sep = "\n"
    )
  }

  # All three p-values for unadjusted model (displayed on plot)
  lrt_p <- tryCatch({
    as.numeric(summary(res$cox_unadj)$logtest["pvalue"])
  }, error = function(e) NA_real_)
  wald_p <- tryCatch({
    s <- summary(res$cox_unadj)
    idx <- which(grepl("^groupHigh$", rownames(s$coefficients)))
    as.numeric(s$coefficients[idx, "Pr(>|z|)"])
  }, error = function(e) NA_real_)
  logrank_p <- tryCatch({
    sdiff <- survival::survdiff(survival::Surv(OS_time, OS_event) ~ group, data = dat)
    as.numeric(1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1))
  }, error = function(e) NA_real_)

  fit <- survival::survfit(survival::Surv(OS_time, OS_event) ~ group, data = dat)

  g <- survminer::ggsurvplot(
    fit, data = dat,
    palette      = c("black", "grey50"),
    linetype     = c("solid", "dashed"),
    conf.int     = FALSE,
    pval         = FALSE,
    risk.table   = TRUE,
    risk.table.col = "black",
    ggtheme      = ggplot2::theme_classic(base_size = 16),
    title        = title_txt,
    subtitle     = sub_txt,
    legend.title = "",
    legend.labs  = c("Low", "High"),
    xlab         = "Time (days)",
    ylab         = "Overall survival probability"
  )

  # Add Log-rank p-value annotation only
  g$plot <- g$plot +
    ggplot2::annotate(
      "text",
      x     = max(dat$OS_time, na.rm = TRUE) * 0.55,
      y     = 0.95,
      label = paste0("Log-rank p = ", signif(logrank_p, 3)),
      size  = 5, hjust = 0
    ) +
    ggplot2::theme(
      panel.border  = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.0),
      axis.line     = ggplot2::element_line(linewidth = 0.9),
      axis.ticks    = ggplot2::element_line(linewidth = 0.9),
      legend.position = "right"
    )

  if (is.null(out_png)) out_png <- paste0("KM_", gene, "_PAAD.png")

  grDevices::png(out_png, width = 8, height = 7, units = "in", res = dpi)
  print(g)
  grDevices::dev.off()

  cat("Saved:", out_png, "\n")
}

# ---------------------------
# 9) Forest plot: Unadjusted vs TME-adjusted HRs (using LRT p-values)
# ---------------------------
plot_forest <- function(tme_table,
                        out_png = "Forest_PAAD_unadj_vs_TMEadj.png",
                        dpi     = 600) {

  # Reshape table into long format for plotting
  # Each gene has 3 rows: Unadjusted, Adj+ESTIMATE, Adj+Stromal
  forest_df <- dplyr::bind_rows(
    tme_table %>%
      dplyr::transmute(
        Gene,
        Model   = "Unadjusted",
        HR      = HR_unadj,
        CI_low  = CI_low_unadj,
        CI_high = CI_high_unadj,
        p       = p_lrt_unadj
      ),
    tme_table %>%
      dplyr::transmute(
        Gene,
        Model   = "Adj + ESTIMATE",
        HR      = HR_adj_estimate,
        CI_low  = CI_low_adj_estimate,
        CI_high = CI_high_adj_estimate,
        p       = p_lrt_adj_estimate
      ),
    tme_table %>%
      dplyr::transmute(
        Gene,
        Model   = "Adj + Stromal",
        HR      = HR_adj_stromal,
        CI_low  = CI_low_adj_stromal,
        CI_high = CI_high_adj_stromal,
        p       = p_lrt_adj_stromal
      )
  ) %>%
    dplyr::filter(!is.na(HR))

  # Order: genes on y-axis, models as color/shape
  forest_df$Gene  <- factor(forest_df$Gene, levels = rev(unique(tme_table$Gene)))
  forest_df$Model <- factor(forest_df$Model,
                            levels = c("Unadjusted", "Adj + ESTIMATE", "Adj + Stromal"))

  # Significance label
  forest_df <- forest_df %>%
    dplyr::mutate(
      sig_label = dplyr::case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE      ~ ""
      ),
      hr_label = paste0(
        sprintf("%.2f", HR),
        " (", sprintf("%.2f", CI_low), "\u2013", sprintf("%.2f", CI_high), ")",
        " ", sig_label
      )
    )

  # Dodge position for multiple models per gene
  pd <- ggplot2::position_dodge(width = 0.6)

  p <- ggplot2::ggplot(forest_df,
                       ggplot2::aes(x = HR, y = Gene, color = Model, shape = Model)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    ggplot2::geom_point(position = pd, size = 3.5) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = CI_low, xmax = CI_high),
      position = pd, height = 0.2, linewidth = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = CI_high + 0.05, label = hr_label),
      position = pd, hjust = 0, size = 3, show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("Unadjusted"     = "black",
                 "Adj + ESTIMATE"  = "#E7298A",
                 "Adj + Stromal"   = "#2E9FDF")
    ) +
    ggplot2::scale_shape_manual(
      values = c("Unadjusted"     = 16,
                 "Adj + ESTIMATE"  = 17,
                 "Adj + Stromal"   = 15)
    ) +
    ggplot2::scale_x_continuous(
      trans  = "log2",
      breaks = c(0.25, 0.5, 1, 2, 4),
      labels = c("0.25", "0.5", "1", "2", "4")
    ) +
    ggplot2::labs(
      title    = "TCGA-PAAD: Cox HR (High vs Low, Q25 vs Q75)",
      subtitle = "LRT p-value | * p<0.05, ** p<0.01, *** p<0.001",
      x        = "Hazard Ratio (log2 scale)",
      y        = "",
      color    = "Model",
      shape    = "Model"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      panel.border    = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      plot.title      = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle   = ggplot2::element_text(size = 10, color = "grey30"),
      axis.text.y     = ggplot2::element_text(face = "bold", size = 12)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 80, 10, 10))

  grDevices::png(out_png, width = 10, height = 6, units = "in", res = dpi)
  print(p)
  grDevices::dev.off()

  cat("Saved:", out_png, "\n")
}

# =============================================================================
# 10) Stratified Cox: strata(estimate_score High/Low) for each gene
#     Uses survival::strata() to allow different baseline hazards by TME group
# =============================================================================
run_stratified_cox <- function(genes,
                               q          = c(0.25, 0.75),
                               tme_var    = "estimate_score",
                               tme_split  = "median",
                               out_csv    = "PAAD_stratified_cox.csv") {

  genes <- unique(genes)
  genes_present <- genes[genes %in% rownames(expr)]
  if (length(genes_present) == 0) stop("None of the genes found in expr.")

  results <- vector("list", length(genes_present))
  names(results) <- genes_present

  for (g in genes_present) {

    # Build per-gene data
    dat <- surv_df %>%
      dplyr::mutate(gene_expr = as.numeric(expr[g, sample])) %>%
      dplyr::filter(!is.na(gene_expr), !is.na(.data[[tme_var]]))

    # Quantile stratification for gene expression
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

    # Split TME variable into High/Low
    if (tme_split == "median") {
      tme_cut <- stats::median(dat[[tme_var]], na.rm = TRUE)
    } else if (tme_split == "mean") {
      tme_cut <- mean(dat[[tme_var]], na.rm = TRUE)
    } else {
      tme_cut <- stats::median(dat[[tme_var]], na.rm = TRUE)
    }

    dat <- dat %>%
      dplyr::mutate(
        tme_group = ifelse(.data[[tme_var]] >= tme_cut, "High", "Low")
      )
    dat$tme_group <- factor(dat$tme_group, levels = c("Low", "High"))

    n_total  <- nrow(dat)
    n_events <- sum(dat$OS_event == 1, na.rm = TRUE)

    # --- Model 1: Unadjusted (no TME) ---
    cox_unadj <- tryCatch(
      survival::coxph(survival::Surv(OS_time, OS_event) ~ group,
                      data = dat),
      error = function(e) NULL
    )

    # --- Model 2: Adjusted (TME as covariate) ---
    fml_adj <- stats::as.formula(
      paste0("survival::Surv(OS_time, OS_event) ~ group + ", tme_var)
    )
    cox_adj <- tryCatch(
      survival::coxph(fml_adj, data = dat),
      error = function(e) NULL
    )

    # --- Model 3: Stratified (strata on TME group) ---
    cox_strata <- tryCatch(
      survival::coxph(
        survival::Surv(OS_time, OS_event) ~ group + survival::strata(tme_group),
        data = dat
      ),
      error = function(e) NULL
    )

    # --- Model 4: Interaction (group * TME continuous) ---
    fml_int <- stats::as.formula(
      paste0("survival::Surv(OS_time, OS_event) ~ group * ", tme_var)
    )
    cox_interact <- tryCatch(
      survival::coxph(fml_int, data = dat),
      error = function(e) NULL
    )

    # Extract stats from each model
    extract_stats <- function(fit, model_label) {
      if (is.null(fit)) {
        return(data.frame(
          Gene = g, Model = model_label,
          N = n_total, Events = n_events,
          HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
          p_lrt = NA_real_, p_wald = NA_real_,
          stringsAsFactors = FALSE
        ))
      }

      s  <- summary(fit)
      rn <- rownames(s$coefficients)
      idx <- which(rn == "groupHigh")
      if (length(idx) == 0) idx <- grep("^group", rn)
      if (length(idx) == 0) {
        return(data.frame(
          Gene = g, Model = model_label,
          N = n_total, Events = n_events,
          HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
          p_lrt = NA_real_, p_wald = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      idx <- idx[1]

      beta <- s$coefficients[idx, "coef"]
      se   <- s$coefficients[idx, "se(coef)"]
      HR   <- exp(beta)
      L    <- exp(beta - 1.96 * se)
      U    <- exp(beta + 1.96 * se)

      p_wald <- as.numeric(s$coefficients[idx, "Pr(>|z|)"])

      p_lrt <- tryCatch({
        d1 <- drop1(fit, test = "Chisq")
        gr <- grep("^group$", rownames(d1))
        if (length(gr) > 0) as.numeric(d1[gr, "Pr(>Chi)"])
        else as.numeric(s$logtest["pvalue"])
      }, error = function(e) as.numeric(s$logtest["pvalue"]))

      data.frame(
        Gene = g, Model = model_label,
        N = n_total, Events = n_events,
        HR = HR, CI_low = L, CI_high = U,
        p_lrt = p_lrt, p_wald = p_wald,
        stringsAsFactors = FALSE
      )
    }

    # Extract interaction p-value separately
    p_interact <- NA_real_
    if (!is.null(cox_interact)) {
      s_int <- summary(cox_interact)
      rn_int <- rownames(s_int$coefficients)
      int_idx <- grep("^group.*:", rn_int)
      if (length(int_idx) == 0) int_idx <- grep(":", rn_int)
      if (length(int_idx) > 0) {
        p_interact <- as.numeric(s_int$coefficients[int_idx[1], "Pr(>|z|)"])
      }
    }

    row_unadj  <- extract_stats(cox_unadj,  "Unadjusted")
    row_adj    <- extract_stats(cox_adj,     paste0("Adjusted + ", tme_var))
    row_strata <- extract_stats(cox_strata,  paste0("Stratified by ", tme_var))
    row_int    <- extract_stats(cox_interact, paste0("Interaction (group*", tme_var, ")"))

    # Add interaction p-value to interaction row
    row_int$p_interaction <- p_interact
    row_unadj$p_interaction  <- NA_real_
    row_adj$p_interaction    <- NA_real_
    row_strata$p_interaction <- NA_real_

    # Add TME cutoff info
    row_unadj$TME_cutoff  <- NA_real_
    row_adj$TME_cutoff    <- NA_real_
    row_strata$TME_cutoff <- tme_cut
    row_int$TME_cutoff    <- NA_real_

    results[[g]] <- dplyr::bind_rows(row_unadj, row_adj, row_strata, row_int)
  }

  out <- dplyr::bind_rows(results)

  # Add HR formatted column
  out <- out %>%
    dplyr::mutate(
      HR_fmt = ifelse(
        is.na(HR), "NA",
        paste0(sprintf("%.2f", HR),
               " (", sprintf("%.2f", CI_low), "\u2013", sprintf("%.2f", CI_high), ")")
      )
    )

  readr::write_csv(out, out_csv)
  cat("\n=== Stratified Cox Summary ===\n")
  print(out %>% dplyr::select(Gene, Model, N, Events, HR_fmt, p_lrt, p_wald, p_interaction))
  cat("Saved:", out_csv, "\n")
  out
}

# =============================================================================
# 11) NOMOGRAM SYSTEM: Survival Prediction, C-index, Calibration, ROC
# =============================================================================

# ---------------------------
# 10a) Build nomogram data: merge gene expression + clinical + TME
# ---------------------------
build_nomogram_data <- function(genes, time_days_max = NULL) {

  genes <- genes[genes %in% rownames(expr)]
  if (length(genes) == 0) stop("No genes found in expression matrix.")

  # Start with surv_df, add gene expression columns (continuous)
  nom_df <- surv_df

  for (g in genes) {
    nom_df[[g]] <- as.numeric(expr[g, nom_df$sample])
  }

  # Convert OS_time to months for clinical interpretability
  nom_df <- nom_df %>%
    dplyr::mutate(
      OS_months = OS_time / 30.44,
      age_years = ifelse(is.na(age_years), NA_real_, age_years)
    )

  # Remove rows with missing critical variables
  required_cols <- c("OS_months", "OS_event", "age_years", "stage_simple",
                     "estimate_score", "stromal_score", genes)
  nom_df <- nom_df %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(required_cols), ~ !is.na(.)))

  # Optional: cap follow-up time
  if (!is.null(time_days_max)) {
    time_months_max <- time_days_max / 30.44
    nom_df <- nom_df %>%
      dplyr::mutate(
        OS_event  = ifelse(OS_months > time_months_max, 0L, OS_event),
        OS_months = pmin(OS_months, time_months_max)
      )
  }

  cat("Nomogram data: N =", nrow(nom_df), "samples,",
      sum(nom_df$OS_event == 1), "events\n")
  nom_df
}

# ---------------------------
# 10b) Fit rms::cph model for nomogram
# ---------------------------
fit_nomogram_model <- function(nom_df, genes,
                               include_clinical = TRUE,
                               include_tme      = TRUE) {

  # Set rms datadist
  dd <- rms::datadist(nom_df)
  options(datadist = "dd")
  assign("dd", dd, envir = .GlobalEnv)

  # Build formula
  gene_terms <- paste(genes, collapse = " + ")
  clin_terms <- if (include_clinical) "age_years + stage_simple" else NULL
  tme_terms  <- if (include_tme) "estimate_score + stromal_score" else NULL

  all_terms <- paste(c(gene_terms, clin_terms, tme_terms), collapse = " + ")
  fml <- stats::as.formula(paste0("Surv(OS_months, OS_event) ~ ", all_terms))

  cat("Nomogram formula:", deparse(fml), "\n")

  # Fit cph model
  fit <- rms::cph(fml, data = nom_df, x = TRUE, y = TRUE, surv = TRUE,
                  time.inc = 12)  # 12 months = 1 year

  cat("Model fitted. N =", fit$stats["Obs"], ", Events =", fit$stats["Events"], "\n")
  fit
}

# ---------------------------
# 10c) Plot nomogram
# ---------------------------
plot_nomogram <- function(fit,
                          time_points = c(12, 24, 36),
                          out_png     = "Nomogram_PAAD.png",
                          dpi         = 600) {

  surv_obj <- rms::Survival(fit)

  nom <- rms::nomogram(
    fit,
    fun = lapply(time_points, function(t) {
      function(x) surv_obj(t, x)
    }),
    fun.at    = seq(0.1, 0.9, by = 0.1),
    funlabel  = paste0(time_points / 12, "-Year OS"),
    lp        = TRUE,
    maxscale  = 100
  )

  grDevices::png(out_png, width = 12, height = 10, units = "in", res = dpi)
  plot(nom, xfrac = 0.3, cex.axis = 0.8, cex.var = 1.0,
       col.grid = gray(c(0.8, 0.95)),
       main = "TCGA-PAAD Prognostic Nomogram")
  grDevices::dev.off()

  cat("Saved:", out_png, "\n")
  nom
}

# ---------------------------
# 10d) C-index with bootstrap 95% CI
# ---------------------------
compute_cindex <- function(fit, nom_df, n_boot = 1000,
                           out_csv = "Nomogram_Cindex_PAAD.csv") {

  # Apparent C-index from model
  c_apparent <- as.numeric(fit$stats["C"])

  # Optimism-corrected C-index via rms::validate
  set.seed(42)
  v <- rms::validate(fit, method = "boot", B = n_boot)

  # Dxy (Somers' D) row: corrected = original - optimism
  dxy_row    <- v["Dxy", ]
  dxy_orig   <- dxy_row["index.orig"]
  dxy_corrected <- dxy_row["index.corrected"]
  optimism   <- dxy_row["optimism"]

  c_corrected <- 0.5 + dxy_corrected / 2

  # Bootstrap CI for C-index
  boot_c <- numeric(n_boot)
  set.seed(42)
  for (i in seq_len(n_boot)) {
    idx <- sample(nrow(nom_df), replace = TRUE)
    boot_dat <- nom_df[idx, ]
    bc <- tryCatch({
      dd_b <- rms::datadist(boot_dat)
      options(datadist = "dd_b")
      fit_b <- rms::cph(fit$sformula, data = boot_dat, x = TRUE, y = TRUE, surv = TRUE)
      as.numeric(fit_b$stats["C"])
    }, error = function(e) NA_real_)
    boot_c[i] <- bc
  }
  # Restore original datadist
  options(datadist = "dd")

  boot_c <- boot_c[!is.na(boot_c)]
  ci_low  <- as.numeric(quantile(boot_c, 0.025))
  ci_high <- as.numeric(quantile(boot_c, 0.975))

  cindex_df <- data.frame(
    Metric           = c("C-index (apparent)", "C-index (optimism-corrected)",
                         "Bootstrap 95% CI lower", "Bootstrap 95% CI upper",
                         "Dxy (original)", "Dxy (corrected)", "Optimism"),
    Value            = c(c_apparent, c_corrected, ci_low, ci_high,
                         dxy_orig, dxy_corrected, optimism),
    stringsAsFactors = FALSE
  )

  readr::write_csv(cindex_df, out_csv)
  cat("\n=== C-index Results ===\n")
  print(cindex_df)
  cat("Saved:", out_csv, "\n")

  cindex_df
}

# ---------------------------
# 10e) Calibration curves (1-yr, 2-yr, 3-yr)
# ---------------------------
plot_calibration <- function(fit, nom_df,
                             time_points = c(12, 24, 36),
                             n_boot      = 200,
                             out_png     = "Calibration_PAAD.png",
                             dpi         = 600) {

  n_tp <- length(time_points)

  grDevices::png(out_png, width = 5 * n_tp, height = 5.5, units = "in", res = dpi)
  par(mfrow = c(1, n_tp), mar = c(5, 5, 3, 1))

  for (tp in time_points) {

    # Refit cph with the correct time.inc
    dd_cal <- rms::datadist(nom_df)
    options(datadist = "dd_cal")
    assign("dd_cal", dd_cal, envir = .GlobalEnv)

    fit_cal <- rms::cph(fit$sformula, data = nom_df,
                        x = TRUE, y = TRUE, surv = TRUE,
                        time.inc = tp)

    set.seed(42)
    cal <- tryCatch(
      rms::calibrate(fit_cal, u = tp, B = n_boot, cmethod = "KM"),
      error = function(e) {
        cat("Calibration failed for time =", tp, "months:", e$message, "\n")
        NULL
      }
    )

    if (!is.null(cal)) {
      plot(cal,
           xlab = "Nomogram-Predicted Probability",
           ylab = "Observed Probability (KM)",
           main = paste0(round(tp / 12, 1), "-Year OS Calibration"),
           subtitles = TRUE,
           cex.subtitle = 0.7)
      abline(0, 1, lty = 2, col = "red", lwd = 1.5)
    } else {
      plot.new()
      title(main = paste0(round(tp / 12, 1), "-Year OS Calibration\n(failed)"))
    }
  }

  grDevices::dev.off()

  # Restore original datadist
  options(datadist = "dd")
  cat("Saved:", out_png, "\n")
}

# ---------------------------
# 10f) Time-dependent ROC/AUC
# ---------------------------
plot_time_roc <- function(fit, nom_df,
                          time_points = c(12, 24, 36),
                          out_png     = "TimeROC_PAAD.png",
                          out_csv     = "TimeROC_AUC_PAAD.csv",
                          dpi         = 600) {

  # Compute linear predictor
  lp <- predict(fit, type = "lp")

  troc <- timeROC::timeROC(
    T          = nom_df$OS_months,
    delta      = nom_df$OS_event,
    marker     = lp,
    cause      = 1,
    weighting  = "marginal",
    times      = time_points,
    iid        = TRUE
  )

  # Extract AUC and CI
  auc_vals <- troc$AUC
  ci_list  <- timeROC::confint(troc, level = 0.95)

  auc_df <- data.frame(
    Time_months = time_points,
    Time_years  = time_points / 12,
    AUC         = as.numeric(auc_vals),
    CI_low      = as.numeric(ci_list$CI_AUC[, 1]),
    CI_high     = as.numeric(ci_list$CI_AUC[, 2]),
    stringsAsFactors = FALSE
  )

  readr::write_csv(auc_df, out_csv)
  cat("\n=== Time-dependent AUC ===\n")
  print(auc_df)
  cat("Saved:", out_csv, "\n")

  # Plot ROC curves
  colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")

  grDevices::png(out_png, width = 7, height = 7, units = "in", res = dpi)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "1 - Specificity (FPR)", ylab = "Sensitivity (TPR)",
       main = "TCGA-PAAD: Time-Dependent ROC Curves")
  abline(0, 1, lty = 2, col = "grey50")

  for (i in seq_along(time_points)) {
    tp_idx <- which(troc$times == time_points[i])
    tp     <- troc$TP[, tp_idx]
    fp     <- troc$FP[, tp_idx]

    lines(fp, tp, col = colors[i], lwd = 2.5)
  }

  legend("bottomright",
         legend = paste0(time_points / 12, "-yr AUC = ",
                         sprintf("%.3f", auc_df$AUC),
                         " (", sprintf("%.3f", auc_df$CI_low),
                         "\u2013", sprintf("%.3f", auc_df$CI_high), ")"),
         col = colors[seq_along(time_points)],
         lwd = 2.5, cex = 0.9, bty = "n")

  grDevices::dev.off()

  cat("Saved:", out_png, "\n")
  auc_df
}

# ---------------------------
# 10g) Master nomogram runner
# ---------------------------
run_nomogram_system <- function(genes          = GENES_TO_TEST,
                                time_points    = c(12, 24, 36),
                                n_boot_cindex  = 1000,
                                n_boot_calib   = 200,
                                include_clinical = TRUE,
                                include_tme      = TRUE,
                                prefix         = "PAAD") {

  cat("\n", strrep("=", 60), "\n")
  cat("  NOMOGRAM SYSTEM: TCGA-PAAD\n")
  cat(strrep("=", 60), "\n\n")

  # Step 1: Build data
  cat("--- Step 1: Building nomogram data ---\n")
  nom_df <- build_nomogram_data(genes)

  # Step 2: Fit model
  cat("\n--- Step 2: Fitting Cox model (rms::cph) ---\n")
  fit <- fit_nomogram_model(nom_df, genes,
                            include_clinical = include_clinical,
                            include_tme      = include_tme)

  cat("\n--- Model Summary ---\n")
  print(fit)

  # Step 3: Nomogram
  cat("\n--- Step 3: Plotting nomogram ---\n")
  nom <- plot_nomogram(fit,
                       time_points = time_points,
                       out_png     = paste0("Nomogram_", prefix, ".png"))

  # Step 4: C-index
  cat("\n--- Step 4: C-index with bootstrap ---\n")
  cindex <- compute_cindex(fit, nom_df,
                           n_boot  = n_boot_cindex,
                           out_csv = paste0("Nomogram_Cindex_", prefix, ".csv"))

  # Step 5: Calibration
  cat("\n--- Step 5: Calibration curves ---\n")
  plot_calibration(fit, nom_df,
                   time_points = time_points,
                   n_boot      = n_boot_calib,
                   out_png     = paste0("Calibration_", prefix, ".png"))

  # Step 6: Time-dependent ROC
  cat("\n--- Step 6: Time-dependent ROC/AUC ---\n")
  auc <- plot_time_roc(fit, nom_df,
                       time_points = time_points,
                       out_png     = paste0("TimeROC_", prefix, ".png"),
                       out_csv     = paste0("TimeROC_AUC_", prefix, ".csv"))

  cat("\n", strrep("=", 60), "\n")
  cat("  NOMOGRAM SYSTEM COMPLETE\n")
  cat(strrep("=", 60), "\n")

  # Return all results
  invisible(list(
    nom_df   = nom_df,
    fit      = fit,
    nomogram = nom,
    cindex   = cindex,
    auc      = auc
  ))
}

# =============================================================================
# 11) Run everything
# =============================================================================
gene_list <- GENES_TO_TEST

# Build comparison table
tme_table <- build_gene_TME_compare_table(
  genes   = gene_list,
  q       = Q_CUTS,
  out_csv = OUT_CSV
)
print(tme_table)

# Generate forest plot
plot_forest(tme_table, out_png = "Forest_PAAD_unadj_vs_TMEadj.png")

# Generate KM plots for each gene
for (g in gene_list) {
  plot_gene_km(gene = g, q = Q_CUTS)
}

# Stratified Cox: estimate_score
strata_est <- run_stratified_cox(
  genes     = gene_list,
  q         = Q_CUTS,
  tme_var   = "estimate_score",
  tme_split = "median",
  out_csv   = "PAAD_stratified_cox_estimate.csv"
)

# Stratified Cox: stromal_score
strata_str <- run_stratified_cox(
  genes     = gene_list,
  q         = Q_CUTS,
  tme_var   = "stromal_score",
  tme_split = "median",
  out_csv   = "PAAD_stratified_cox_stromal.csv"
)

# Run nomogram system (genes + clinical + TME)
nom_results <- run_nomogram_system(
  genes            = gene_list,
  time_points      = c(12, 24, 36),   # 1-yr, 2-yr, 3-yr
  n_boot_cindex    = 1000,
  n_boot_calib     = 200,
  include_clinical = TRUE,             # age_years + stage_simple
  include_tme      = TRUE,             # estimate_score + stromal_score
  prefix           = "PAAD"
)

cat("\n=== All done! ===\n")
cat("Table saved to:", OUT_CSV, "\n")
cat("Forest plot:    Forest_PAAD_unadj_vs_TMEadj.png\n")
cat("KM plots:       KM_<gene>_PAAD.png\n")
cat("Strata Cox:     PAAD_stratified_cox_estimate.csv\n")
cat("Strata Cox:     PAAD_stratified_cox_stromal.csv\n")
cat("Nomogram:       Nomogram_PAAD.png\n")
cat("C-index:        Nomogram_Cindex_PAAD.csv\n")
cat("Calibration:    Calibration_PAAD.png\n")
cat("Time-dep ROC:   TimeROC_PAAD.png\n")
cat("Time-dep AUC:   TimeROC_AUC_PAAD.csv\n")
