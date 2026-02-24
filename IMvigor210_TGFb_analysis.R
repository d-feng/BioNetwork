# ============================================================
# IMvigor210 TGF-β / PD-L1 Blockade Analysis
#
# Source: Mariathasan et al. (2018) Nature 554:544-548
#   "TGF-β attenuates tumour response to PD-L1 blockade
#    by contributing to exclusion of T cells"
#
# Data/code package: IMvigor210CoreBiologies
#   http://research-pub.gene.com/IMvigor210CoreBiologies
#   https://github.com/SiYangming/IMvigor210CoreBiologies
#
# Cohort: 368 metastatic urothelial carcinoma (mUC) patients
#   treated with atezolizumab (anti-PD-L1) in IMvigor210 trial.
#   RNA-seq: TruSeq RNA Access (Illumina)
#   Mutations: whole exome sequencing (WES)
# ============================================================

# ── 0. Install & load ────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "limma", "edgeR",
  "ComplexHeatmap", "circlize",
  "survival", "ggplot2", "dplyr",
  "reshape2", "knitr"
))

# Install IMvigor210CoreBiologies from Genentech
# Download from: http://research-pub.gene.com/IMvigor210CoreBiologies
install.packages("IMvigor210CoreBiologies_1.0.0.tar.gz", repos = NULL)

library(IMvigor210CoreBiologies)
library(DESeq2); library(limma); library(edgeR)
library(ComplexHeatmap); library(circlize)
library(ggplot2); library(dplyr); library(reshape2)
library(survival)

# ── 1. Load datasets ─────────────────────────────────────────
data(cds)                  # CountDataSet: RNA-seq counts + pData (sample metadata)
data(fmone)                # NChannelSet: somatic mutations from WES
data(human_gene_signatures)# Named list of gene sets (TGF-β, Teff, etc.)
data(color_palettes)       # Consistent colour palettes used across figures

# Key metadata columns in pData(cds):
#   "Best Confirmed Overall Response" : CR / PR / SD / PD / NE
#   "binaryResponse"                  : TRUE = CR or PR
#   "FMOne mutation burden per MB"    : tumour mutation burden (TMB)
#   "Lund2"                           : Lund molecular subtype (5 groups)
#   "TCGA Subtype"                    : TCGA bladder cancer subtype
#   "IC Level" / "TC Level"           : PD-L1 IHC expression (0/1/2/3)

irf      <- "Best Confirmed Overall Response"
ml       <- "FMOne mutation burden per MB"
tablesDir <- system.file("tables", package = "IMvigor210CoreBiologies")


# ── 2. Gene signature scoring ─────────────────────────────────
# Pipeline: raw counts → voom normalisation → per-gene z-score → mean(z) per signature
#
# human_gene_signatures contains:
#   Pan-F-TBRS  : pan-fibroblast TGF-β response signature (stromal exclusion)
#   Teff        : effector T cell signature
#   Ma-TBRS     : macrophage TGF-β response signature
#   T-TBRS      : T cell TGF-β response signature
#   EMT1/2/3    : epithelial-mesenchymal transition signatures

goi   <- names(human_gene_signatures)

# Normalise with voom (keeps genes expressed in ≥10% of samples at CPM ≥0.25)
voomD <- filterNvoom(counts(cds),
                     minSamples = ncol(counts(cds)) / 10,
                     minCpm     = 0.25)
m     <- voomD$E
m     <- t(scale(t(m), center = TRUE, scale = TRUE))   # per-gene z-score across patients

# Map Ensembl IDs → gene symbols
geneNames <- setNames(fData(cds)$Symbol, rownames(fData(cds)))
m2        <- m
rownames(m2) <- geneNames[rownames(m2)]

# Score each patient for every signature (mean z-score of signature genes)
for (sig in goi) {
  genes <- human_gene_signatures[[sig]]
  genes <- genes[genes %in% rownames(m2)]
  pData(cds)[, sig] <- gsScore(m2[genes, , drop = FALSE])
}


# ── 3. Figure 3 — Heatmap: immune subtypes × response ────────
# Ordered by Lund subtype → response; annotated with TMB, mutations, PD-L1

pdata  <- pData(cds)

# Map Lund subtypes to 5-group labels
pdata$Lund3 <- factor(pdata$Lund2,
  labels = c("UroA", "GU", "Inf", "UroB", "SCCL"))

# Drop "NE" (not evaluable) samples for response analysis
pdata2 <- pdata[pdata[, irf] != "NE", ]
pdata2[, irf] <- droplevels(pdata2[, irf])

# Gene panel for heatmap rows (curated pathway genes)
subtypeHeatSig2 <- read.csv(file.path(tablesDir, "heatmap_features.csv"),
                             as.is = TRUE)
subtypeHeatSig2 <- subtypeHeatSig2[subtypeHeatSig2$genes %in% rownames(m2), ]

matT     <- m2[subtypeHeatSig2$genes, rownames(pdata2)]
matSplit <- factor(subtypeHeatSig2$pathway,
                   levels = unique(subtypeHeatSig2$pathway),
                   labels = LETTERS[seq_along(unique(subtypeHeatSig2$pathway))])

# Sort columns: Lund subtype → response
sOrder <- order(pdata2$Lund2, pdata2$Lund, pdata2[, irf])
pdata2 <- pdata2[sOrder, ]
matT   <- matT[, sOrder]

# Mutation annotation (WES)
dnaGoi <- c("TP53", "RB1", "FGFR3", "CDKN2A", "ERBB2", "PIK3CA")
mutGoi <- any_mutation(fmone[dnaGoi, ])
mutGoi <- apply(mutGoi, 2, function(x) ifelse(x, "mutant", "non-mutant"))
colnames(mutGoi) <- as.character(pData(fmone)$ANONPT_ID)
mt     <- match(as.character(pdata2$ANONPT_ID), colnames(mutGoi))
mutGoi <- mutGoi[, mt]

# Fisher test: RB1 and TP53 enrichment in GU vs SCCL
for (g in c("RB1", "TP53")) {
  tmp  <- table(mutGoi[g, ], pdata2$Lund2)
  tmp2 <- tmp[, c("Genomically unstable", "Basal/SCC-like")]
  cat(sprintf("Fisher P (%s — GU vs SCCL): %s\n", g,
              signif(fisher.test(tmp2)$p.value, 3)))
}

# Build heatmap annotation bar
ha <- HeatmapAnnotation(
  Lund     = pdata2$Lund3,
  TCGA     = pdata2$"TCGA Subtype",
  IC       = pdata2$"IC Level",
  TC       = pdata2$"TC Level",
  Response = pdata2[, irf],
  TMB      = anno_barplot(as.numeric(pdata2[, ml]),
               border = FALSE, gp = gpar(fill = "black")),
  TP53     = mutGoi["TP53", ],
  RB1      = mutGoi["RB1",  ],
  FGFR3    = mutGoi["FGFR3", ],
  CDKN2A   = mutGoi["CDKN2A", ],
  ERBB2    = mutGoi["ERBB2",  ],
  PIK3CA   = mutGoi["PIK3CA", ],
  col = list(
    Response = color_palettes$irf_palette,
    Lund     = color_palettes$lund_palette3,
    TCGA     = color_palettes$tcga_palette,
    IC       = color_palettes$ic_palette,
    TC       = color_palettes$tc_palette,
    TP53     = c(mutant = "black", "non-mutant" = "white"),
    RB1      = c(mutant = "black", "non-mutant" = "white"),
    FGFR3    = c(mutant = "black", "non-mutant" = "white"),
    CDKN2A   = c(mutant = "black", "non-mutant" = "white"),
    ERBB2    = c(mutant = "black", "non-mutant" = "white"),
    PIK3CA   = c(mutant = "black", "non-mutant" = "white")
  )
)

# Draw heatmap
heat_colors <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))

ht <- Heatmap(limitRange(matT),
  name             = "Expression",
  top_annotation   = ha,
  col              = heat_colors,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  show_column_names = FALSE,
  row_names_gp     = gpar(fontsize = 5),
  split            = matSplit,
  gap              = unit(1, "mm"))

draw(ht, padding = unit(c(2, 0, 7, 0), "mm"))

# Add subtype boundary lines
for (xpos in c(0.30, 0.48, 0.765, 0.815)) {
  decorate_heatmap_body("Expression", {
    grid.lines(c(xpos, xpos), unit(c(-24.8, 7.4), "native"),
               gp = gpar(col = "black", lwd = 1))
  })
}


# ── 4. Figure 3 inset — Response rate by Lund subtype ────────
# Logistic regression: does subtype predict binary response?

tmpDat <- pdata[pdata[, irf] != "NE", ]
tmpDat[, irf] <- droplevels(tmpDat[, irf])

fit0 <- glm(binaryResponse ~ 1,     data = tmpDat, family = "binomial")
fit1 <- glm(binaryResponse ~ Lund2, data = tmpDat, family = "binomial")
print(anova(fit0, fit1, test = "Chisq"))    # likelihood ratio test

# GU vs all others (Fisher)
tmpDat$group <- ifelse(tmpDat$Lund2 == "Genomically unstable", "GU", "Others")
a <- table(tmpDat$binaryResponse, tmpDat$group)
cat(sprintf("Fisher P (binary response, GU vs other): %s\n",
            signif(fisher.test(a)$p.value, 2)))

# Stacked barplot: response fraction per subtype
tmpDat$Lund3 <- factor(tmpDat$Lund2,
  labels = c("UroA", "GU", "Inf", "UroB", "SCCL"))
ic <- prop.table(table(tmpDat[, irf], tmpDat$Lund3), margin = 2)
barplot(ic,
  ylab   = "Fraction of patients",
  col    = color_palettes$irf_palette,
  legend.text = rownames(ic),
  args.legend = list(bty = "n", x = "topright"))


# ── 5. Figure 4 — Mouse model: TGF-β blockade + anti-PD-L1 ──
# Four treatment arms: Vehicle, aPDL1, 1D11 (anti-TGFβ), 1D11+aPDL1
# Key question: does combo (1D11+aPDL1) > either monotherapy?

scores_raw <- readRDS(system.file("mouse", "scores.rds",
                                  package = "IMvigor210CoreBiologies"))
pheno      <- readRDS(system.file("mouse", "pheno.rds",
                                  package = "IMvigor210CoreBiologies"))

scores <- tibble::as_tibble(t(scores_raw))
colnames(scores) <- c("Pan-F-TBRS", "Teff", "Ma-TBRS", "T-TBRS",
                       "EMT1", "EMT2", "EMT3")

# Key derived biomarker: Teff minus TGF-β fibroblast score
# High Teff-FTBRS = T cells present AND TGF-β stromal signal low → good prognosis
scores <- scores %>%
  mutate(`Teff-FTBRS` = Teff - `Pan-F-TBRS`,
         Treatment    = factor(pheno$Drug,
           levels = c("Vehicle", "aPDL1", "1D11", "1D11+aPDL1")))

scoresMelt <- melt(scores, id.vars = "Treatment")

# Pairwise Wilcoxon tests for each signature across all treatment pairs
pairs <- list(
  c("Vehicle",  "1D11"),
  c("Vehicle",  "aPDL1"),
  c("Vehicle",  "1D11+aPDL1"),
  c("1D11",     "aPDL1"),
  c("1D11",     "1D11+aPDL1"),
  c("aPDL1",    "1D11+aPDL1")
)

sigs <- setdiff(colnames(scores), "Treatment")

for (sig in sigs) {
  cat(sprintf("\n── %s ──\n", sig))
  dat <- dplyr::filter(scoresMelt, variable == sig)
  for (pair in pairs) {
    tmp <- dplyr::filter(dat, Treatment %in% pair)
    p   <- wilcox.test(value ~ Treatment, data = tmp)$p.value
    cat(sprintf("  %s vs %s : p = %s\n", pair[1], pair[2], signif(p, 3)))
  }

  # Boxplot per signature
  g <- ggplot(dat, aes(x = Treatment, y = value, fill = Treatment)) +
    geom_boxplot(outlier.size = 0.8) +
    theme_bw(base_size = 11) +
    labs(title = sig, x = NULL, y = "Score") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
  print(g)
}


# ── 6. Key biomarker summary ──────────────────────────────────
# The Teff-FTBRS ratio is the central finding:
#
#   High Teff + Low Pan-F-TBRS → responds to aPDL1 alone
#   Low  Teff + High Pan-F-TBRS → T cell exclusion phenotype;
#                                  requires TGFβ blockade (1D11) + aPDL1
#
# This motivates the combination therapy (bintrafusp alfa: TGFβtrap + PDL1 Ab)

cat("\nDone. Key result: combo 1D11 + aPDL1 rescues TGF-β-excluded tumours.\n")
