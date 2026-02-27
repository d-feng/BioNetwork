# ============================================================
# IMvigor210 TGF-β / PD-L1 Blockade Analysis
#
# Source: Mariathasan et al. (2018) Nature 554:544-548
#   "TGF-β attenuates tumour response to PD-L1 blockade
#    by contributing to exclusion of T cells"
#
# Original Rmd authors:
#   Figure 3: Dorothee Nickles
#   Figure 4: Yasin Senbabaoglu
#
# Local package path (update if package is elsewhere):
#   PKG_DIR below
# ============================================================

# ── 0. Paths ──────────────────────────────────────────────────
PKG_DIR   <- "C:/Users/difen/Rcode/IMvigor210CoreBiologies_1.0.1/IMvigor210CoreBiologies"
tablesDir <- file.path(PKG_DIR, "inst", "tables")
mouseDir  <- file.path(PKG_DIR, "inst", "mouse")

OUT_DIR <- file.path(PKG_DIR, "output_figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("Figures will be saved to: ", OUT_DIR)

# ── 1. Install dependencies (run once) ────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "limma", "edgeR",
                       "ComplexHeatmap", "circlize", "biomaRt"))
install.packages(c("ggplot2", "dplyr", "reshape2", "knitr", "DT", "spatstat"))

# Install IMvigor210CoreBiologies from local source
# (lsmeans removed from DESCRIPTION — no longer on CRAN)
setwd(PKG_DIR)
install.packages(".", repos = NULL, type = "source")

# ── 2. CountDataSet stub (run once per R installation) ────────
# data(cds) loads a CountDataSet — S4 class from the retired DESeq
# package (removed from Bioconductor 3.13). Biobase's .requirePackage()
# checks for a namespace called "DESeq", so we build a minimal stub
# package that defines CountDataSet. Installed once; persists across
# sessions.
if (!requireNamespace("DESeq", quietly = TRUE)) {
  message("Building DESeq compatibility stub...")
  stub <- file.path(tempdir(), "DESeq")
  dir.create(file.path(stub, "R"), recursive = TRUE, showWarnings = FALSE)

  writeLines(c(
    "Package: DESeq",
    "Version: 1.36.0",
    "Title: CountDataSet Compatibility Stub",
    "Description: Minimal stub defining the CountDataSet S4 class so",
    "    that IMvigor210CoreBiologies data(cds) loads without the",
    "    retired DESeq package.",
    "Depends: R (>= 3.3), Biobase",
    "License: LGPL-3",
    "NeedsCompilation: no"
  ), con = file.path(stub, "DESCRIPTION"))

  writeLines(c(
    "import(Biobase)",
    "exportClasses(CountDataSet)"
  ), con = file.path(stub, "NAMESPACE"))

  writeLines(c(
    "# CountDataSet: eSet subclass from the legacy DESeq package.",
    "# Extends eSet with two extra slots; sizeFactors live in pData.",
    "setClass('CountDataSet',",
    "  contains = 'eSet',",
    "  representation(",
    "    fitInfo   = 'environment',",
    "    dispTable = 'character'",
    "  )",
    ")"
  ), con = file.path(stub, "R", "CountDataSet.R"))

  install.packages(stub, repos = NULL, type = "source")
  message("DESeq stub installed.")
}

# ── 3. Libraries ──────────────────────────────────────────────
library(IMvigor210CoreBiologies)
library(Biobase)
library(DESeq)       # CountDataSet class; satisfies .requirePackage("DESeq")
library(DESeq2); library(limma); library(edgeR)
library(ComplexHeatmap); library(circlize)
library(ggplot2); library(dplyr); library(reshape2)
library(knitr)


# ============================================================
# FIGURE 3  (original author: Dorothee Nickles)
# Heatmap: immune subtypes × clinical response
# + Extended Data Figure 3
# ============================================================

# ── Settings ──────────────────────────────────────────────────
data(human_gene_signatures)
ind_genes <- human_gene_signatures

data(color_palettes)

irf       <- "Best Confirmed Overall Response"
ml        <- "FMOne mutation burden per MB"
goi       <- names(ind_genes)

labCex    <- 0.9
namesCex  <- 0.9
legendCex <- 0.9
titleCex  <- 1
axisCex   <- 0.9
titleF    <- 1

# ── Data preparation ──────────────────────────────────────────
data(cds)
cds2 <- cds
data(fmone)
fmi  <- fmone

geneNames <- setNames(fData(cds2)$Symbol,
                      as.character(rownames(fData(cds2))))

# assayData()[["counts"]] replaces DESeq's counts() accessor
cnt   <- assayData(cds2)[["counts"]]
voomD <- filterNvoom(cnt,
                     minSamples = ncol(cnt) / 10,
                     minCpm     = 0.25)
m  <- voomD$E
m  <- t(scale(t(m), center = TRUE, scale = TRUE))

m2           <- m
rownames(m2) <- geneNames[rownames(m2)]

# Calculate gene set scores
for (sig in goi) {
  pData(cds2)[, sig] <- NA
  genes <- ind_genes[[sig]]
  genes <- genes[genes %in% rownames(m2)]
  tmp   <- m2[genes, , drop = FALSE]
  pData(cds2)[, sig] <- gsScore(tmp)
}


# ── Figure 3 heatmap ──────────────────────────────────────────
dropNE  <- TRUE
pdata   <- pData(cds2)

pdata$Lund <- factor(pdata$Lund,
  levels = c("MS1a","MS1b","MS2a1","MS2a2","MS2b1","MS2b2.1","MS2b2.2"))
pdata2       <- pdata
pdata2$Lund3 <- factor(pdata2$Lund2,
  labels = c("UroA","GU","Inf","UroB","SCCL"))

if (dropNE) {
  ind           <- pdata2[, irf] != "NE"
  pdata2        <- pdata2[ind, ]
  pdata2[, irf] <- droplevels(pdata2[, irf])
}

subtypeHeatSig2 <- read.csv(file.path(tablesDir, "heatmap_features.csv"),
                             as.is = TRUE)
subtypeHeatSig2 <- subtypeHeatSig2[subtypeHeatSig2$genes %in% rownames(m2), ]

matT     <- m2[subtypeHeatSig2$genes, rownames(pdata2)]
matSplit <- subtypeHeatSig2$pathway
matSplit <- factor(matSplit,
  levels = unique(matSplit),
  labels = LETTERS[1:length(unique(matSplit))])

sOrder <- order(pdata2$"Lund2", pdata2$"Lund", pdata2[, irf])
pdata2 <- pdata2[sOrder, ]
matT   <- matT[, sOrder]

# Mutation annotation
dnaGoi           <- c("TP53","RB1","FGFR3","CDKN2A","ERBB2","PIK3CA")
fmiGoi           <- fmi[dnaGoi, ]
mutGoi           <- any_mutation(fmiGoi)
mutGoi           <- apply(mutGoi, 2, function(x) ifelse(x, "mutant", "non-mutant"))
colnames(mutGoi) <- as.character(pData(fmi)$"ANONPT_ID")
mt               <- match(as.character(pdata2$"ANONPT_ID"), colnames(mutGoi))
mutGoi           <- mutGoi[, mt]

# Fisher tests: RB1 and TP53 enrichment in GU vs SCCL
for (gene in c("RB1", "TP53")) {
  tmp  <- table(mutGoi[gene, ], pdata2$Lund2)
  tmp2 <- tmp[, c("Genomically unstable", "Basal/SCC-like")]
  print(paste("Fisher P-value for", gene, "mutants GU vs SCCL:",
              signif(fisher.test(tmp2)$p.value)))
  print(prop.table(tmp, 2)); print(prop.table(tmp, 1))
}

# Sample annotations
ha <- HeatmapAnnotation(
  Lund     = pdata2$Lund3,
  TCGA     = pdata2$"TCGA Subtype",
  IC       = pdata2$"IC Level",
  TC       = pdata2$"TC Level",
  Response = pdata2[, irf],
  TMB      = anno_barplot(as.numeric(pdata2[, ml]),
               border = FALSE, gp = gpar(fill = "black")),
  TP53     = mutGoi["TP53",   ],
  RB1      = mutGoi["RB1",    ],
  FGFR3    = mutGoi["FGFR3",  ],
  CDKN2A   = mutGoi["CDKN2A", ],
  ERBB2    = mutGoi["ERBB2",  ],
  PIK3CA   = mutGoi["PIK3CA", ],
  annotation_height = unit.c(
    rep(unit(0.3, "cm"), 5),
    unit(0.6, "cm"),
    unit.c(rep(unit(0.3, "cm"), 6))),
  annotation_legend_param = list(
    labels_gp = gpar(fontsize = 9),
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    ncol      = 2),
  gap = unit(c(1, rep(0, 4), 1, rep(0, 5)), "mm"),
  col = list(
    IC       = color_palettes$ic_palette,
    TC       = color_palettes$tc_palette,
    Response = color_palettes$irf_palette,
    Lund     = color_palettes$lund_palette3,
    TCGA     = color_palettes$tcga_palette,
    TP53     = c(mutant = "black", "non-mutant" = "white"),
    RB1      = c(mutant = "black", "non-mutant" = "white"),
    FGFR3    = c(mutant = "black", "non-mutant" = "white"),
    CDKN2A   = c(mutant = "black", "non-mutant" = "white"),
    ERBB2    = c(mutant = "black", "non-mutant" = "white"),
    PIK3CA   = c(mutant = "black", "non-mutant" = "white")),
  show_annotation_name = TRUE,
  annotation_name_gp   = gpar(fontsize = 7))

heat_colors2 <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))

ht1org <- Heatmap(limitRange(matT),
  name                 = "Expression",
  top_annotation       = ha,
  cluster_rows         = FALSE,
  col                  = heat_colors2,
  color_space          = "RGB",
  cluster_columns      = FALSE,
  row_order            = NULL,
  column_order         = NULL,
  show_column_names    = FALSE,
  row_names_gp         = gpar(fontsize = 5),
  split                = matSplit,
  gap                  = unit(1, "mm"),
  column_title         = "",
  column_title_gp      = gpar(fontsize = 5),
  width                = unit(8, "cm"),
  show_heatmap_legend  = TRUE,
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 9),
    title_gp  = gpar(fontsize = 9, fontface = "bold")))

pdf(file.path(OUT_DIR, "Fig3_heatmap.pdf"), width = 12, height = 10)
draw(ht1org, padding = unit(c(2, 0, 7, 0), "mm"))
uLim <- 7.4; lLim <- -24.8
for (xpos in c(0.30, 0.48, 0.765, 0.815)) {
  decorate_heatmap_body("Expression", {
    grid.lines(c(xpos, xpos), unit(c(lLim, uLim), "native"),
               gp = gpar(col = "black", lwd = 1))
  })
}
decorate_heatmap_body("Expression", {
  grid.text("UroA",  0.15,   8.7, default.units = "native", gp = gpar(fontsize = 9))
  grid.text("GU",    0.39,   8.7, default.units = "native", gp = gpar(fontsize = 9))
  grid.text("Inf",   0.6225, 8.7, default.units = "native", gp = gpar(fontsize = 9))
  grid.text("SCCL",  0.9075, 8.7, default.units = "native", gp = gpar(fontsize = 9))
})
dev.off()
message("Saved: Fig3_heatmap.pdf")


# ── Figure 3 inset — response rate by Lund subtype ───────────
tmpDat        <- pdata[pdata[, irf] != "NE", ]
tmpDat[, irf] <- droplevels(tmpDat[, irf])

fit_resp_0 <- glm(tmpDat$binaryResponse ~ 1,           family = "binomial")
fit_resp_l <- glm(tmpDat$binaryResponse ~ tmpDat$Lund2, family = "binomial")
aFit       <- anova(fit_resp_0, fit_resp_l, test = "Chisq")
print(aFit)
pvalAll <- signif(aFit$"Pr(>Chi)"[2])
print(summary(fit_resp_l))

tmpDat$group <- ifelse(tmpDat$Lund2 == "Genomically unstable", "GU", "Others")
a            <- table(tmpDat$binaryResponse, tmpDat$group)
print(prop.table(a, 2))
pvalF <- signif(fisher.test(a)$p.value, 2)
print(paste("Fisher P binary response GU versus other Lund:", pvalF))

for (rr in c("CR", "PR")) {
  tmpDat$group2 <- factor(ifelse(tmpDat[, irf] == rr, rr, "Others"))
  fit_resp_0    <- glm(tmpDat$group2 ~ 1,           family = "binomial")
  fit_resp_l    <- glm(tmpDat$group2 ~ tmpDat$Lund2, family = "binomial")
  print(paste(rr, "vs. others"))
  print(anova(fit_resp_0, fit_resp_l, test = "Chisq"))
  print(summary(fit_resp_l))
}

pdf(file.path(OUT_DIR, "Fig3_inset_response_Lund.pdf"), width = 7, height = 5)
oldMar <- par()$mar
par(mar = c(5.5, 4.1, 2, 0.5))

tmpDat$Lund3 <- factor(tmpDat$Lund2,
  labels = c("UroA", "GU", "Inf", "UroB", "SCCL"))
ic       <- table(tmpDat$Lund3, tmpDat[, irf])
nSamples <- rowSums(ic)
ic       <- prop.table(t(ic), margin = 2)
bWidth   <- 0.1

b    <- barplot(ic,
  cex.axis = axisCex, cex.names = namesCex, cex.lab = labCex,
  ylab = "fraction of patients",
  legend.text = rownames(ic),
  col  = color_palettes$irf_palette,
  width = bWidth, xlim = c(0, 1),
  args.legend = list(bty = "n", cex = legendCex, x = 0, y = 1),
  plot = FALSE)
xLim <- b[1] * 2 + b[5] + bWidth * 2

a <- barplot(ic,
  cex.axis = axisCex, cex.names = namesCex, cex.lab = labCex,
  ylab = "fraction of patients",
  legend.text = rownames(ic),
  col  = color_palettes$irf_palette,
  width = bWidth, xlim = c(0, 1),
  args.legend = list(bty = "n", cex = legendCex, x = xLim, y = 1),
  xaxt = "n")
axis(1, at = a, labels = FALSE)
text(x = a, par("usr")[3] - 0.06,
     labels = levels(tmpDat$Lund3),
     srt = -45, xpd = TRUE, adj = 0, cex = namesCex)
mtext(nSamples, side = 3, at = a, line = 0, cex = 0.6)
mtext("Response, Lund", side = 3, at = a[3], line = 1,
      cex = titleCex, font = titleF)
par(mar = oldMar)
dev.off()
message("Saved: Fig3_inset_response_Lund.pdf")


# ── ED Figure 3b — signature scores by Lund subtype ──────────
tmps       <- pData(cds2)
tmps$PPARG <- scale(m2["PPARG", ], center = TRUE, scale = TRUE)

for (sig in c("FGFR3 related", "WNT target", "PPARG")) {
  sig_safe <- gsub("[^A-Za-z0-9]", "_", sig)
  pdf(file.path(OUT_DIR, paste0("EDFig3b_", sig_safe, ".pdf")), width = 6, height = 5)
  par(mar = c(5.5, 4.1, 2, 2))
  feat   <- "Lund2"
  tmpDat <- tmps
  tmpDat[, feat] <- factor(tmpDat[, feat],
    levels = levels(tmpDat[, feat]),
    labels = c("UroA", "GU", "Inf", "UroB", "SCCL"))
  a <- boxplot(tmpDat[, sig] ~ tmpDat[, feat],
    ylab     = paste(sig, "score"),
    col      = color_palettes[["lund_palette3"]][levels(tmpDat[, feat])],
    cex.axis = axisCex, cex.lab = labCex, cex.names = namesCex,
    whisklty = 1, xaxt = "n")
  axis(1, at = 1:nlevels(tmpDat[, feat]),
       labels = rep("", nlevels(tmpDat[, feat])))
  yrange <- par("usr")[4] - par("usr")[3]
  yunit  <- yrange / 60
  text(x = 1:nlevels(tmpDat[, feat]),
       y = par("usr")[3] - yunit * 4,
       labels = levels(tmpDat[, feat]),
       srt = -45, xpd = TRUE, adj = 0, cex = namesCex)
  pval <- signif(getPfromAnova(tmpDat[, sig], tmpDat[, feat]), 2)
  print(paste("P for", sig, "by Lund subtype:", pval))
  mtext(a$n, side = 3, at = 1:nlevels(tmpDat[, feat]), line = 0, cex = 0.8)
  mtext(paste(sig, "subtype", sep = ", "),
        side = 3, at = 3, line = 1, font = titleF, cex = titleCex)
  dev.off()
  message("Saved: EDFig3b_", sig_safe, ".pdf")
}


# ── ED Figure 3c — immune phenotype by Lund subtype ──────────
pdf(file.path(OUT_DIR, "EDFig3c_immunophenotype_Lund.pdf"), width = 7, height = 5)
feat   <- "Lund2"
oldMar <- par()$mar
par(mar = c(5.5, 4.1, 2, 0))
cols   <- color_palettes[["lund_palette3"]]
tmpDat <- pData(cds2)[!is.na(pData(cds2)$"Immune phenotype"), ]
tmpDat[, feat] <- factor(tmpDat[, feat],
  levels = rev(levels(tmpDat[, feat])))
ic <- table(tmpDat$"Immune phenotype", tmpDat[, feat])

print(prop.table(ic, 1)); print(prop.table(ic, 2))
pval <- signif(chisq.test(ic)$p.value, 2)
print(paste("Chisquared P Immunophenotype by Lund:", pval))

nSamples     <- rowSums(ic)
ic           <- prop.table(t(ic), margin = 2)
rownames(ic) <- sub("Basal/SCC-like",      "SCCL", rownames(ic))
rownames(ic) <- sub("iltrated",            "",     rownames(ic))
rownames(ic) <- sub("Genomically unstable","GU",   rownames(ic))

a <- barplot(ic,
  ylab = "fraction of patients",
  legend.text = rownames(ic),
  col  = cols[rownames(ic)],
  width = 0.16, xlim = c(0, 1),
  args.legend = list(bty = "n", cex = legendCex, x = "topright"),
  xaxt = "n",
  cex.axis = axisCex, cex.names = namesCex, cex.lab = labCex)
axis(1, at = a, labels = FALSE)
text(x = a, par("usr")[3] - 0.06,
     labels = levels(tmpDat$"Immune phenotype"),
     srt = -45, xpd = TRUE, adj = 0, cex = namesCex)
mtext(nSamples, side = 3, at = a, line = 0, cex = 0.75)
mtext("Lund, phenotype", side = 3, at = a[2] + 0.02,
      line = 0.7, cex = titleCex, font = titleF)
par(mar = oldMar)
dev.off()
message("Saved: EDFig3c_immunophenotype_Lund.pdf")


# ============================================================
# FIGURE 4  (original author: Yasin Senbabaoglu)
# Mouse model: EMT6 tumours, 4 treatment arms
# Panels k,l + Extended Data Figure A11 panels f,g,h
# ============================================================

x     <- readRDS(file.path(mouseDir, "scores.rds"))
pheno <- readRDS(file.path(mouseDir, "pheno.rds"))

scores <- tbl_df(t(x))
colnames(scores) <- c("Pan-F-TBRS", "Teff", "Ma-TBRS",
                      "T-TBRS", "EMT1", "EMT2", "EMT3")

# Derived biomarker: Teff minus TGF-β fibroblast score
# High Teff-FTBRS = T cells present + TGFb stromal signal low  -> responds to aPDL1 alone
# Low  Teff-FTBRS = T cell exclusion                           -> requires 1D11 + aPDL1
scores <- scores %>%
  dplyr::select(Teff, everything()) %>%
  mutate(`Teff-FTBRS` = Teff - `Pan-F-TBRS`) %>%
  mutate(Treatment    = pheno$Drug) %>%
  dplyr::select(Treatment, everything())

scoresMelt <- melt(scores)
sname      <- colnames(scores)[-1]

for (i in seq_along(sname)) {
  dat <- filter(scoresMelt, variable == sname[i])

  pvec    <- rep(NA, 6)
  pvec[1] <- signif(wilcox.test(value ~ Treatment,
    data = filter(dat, Treatment %in% c("Vehicle",  "1D11")))$p.value,       3)
  pvec[2] <- signif(wilcox.test(value ~ Treatment,
    data = filter(dat, Treatment %in% c("Vehicle",  "aPDL1")))$p.value,      3)
  pvec[3] <- signif(wilcox.test(value ~ Treatment,
    data = filter(dat, Treatment %in% c("Vehicle",  "1D11+aPDL1")))$p.value, 3)
  pvec[4] <- signif(wilcox.test(value ~ Treatment,
    data = filter(dat, Treatment %in% c("1D11",     "aPDL1")))$p.value,      3)
  pvec[5] <- signif(wilcox.test(value ~ Treatment,
    data = filter(dat, Treatment %in% c("1D11",     "1D11+aPDL1")))$p.value, 3)
  pvec[6] <- signif(wilcox.test(value ~ Treatment,
    data = filter(dat, Treatment %in% c("aPDL1",    "1D11+aPDL1")))$p.value, 3)

  cat("\n###", sname[i], "\n")
  dat$Treatment <- factor(dat$Treatment,
    levels = c("Vehicle", "aPDL1", "1D11", "1D11+aPDL1"))

  g <- ggplot(data = dat, aes(x = Treatment, y = value, fill = Treatment)) +
    geom_boxplot() + theme_bw() + labs(x = "", y = "Score", title = sname[i])
  print(g)
  sig_safe <- gsub("[^A-Za-z0-9]", "_", sname[i])
  ggsave(file.path(OUT_DIR, paste0("Fig4_", sig_safe, ".pdf")),
         plot = g, width = 6, height = 4)
  message("Saved: Fig4_", sig_safe, ".pdf")

  cat("Vehicle vs 1D11 p-value =",         pvec[1], "\n\n")
  cat("Vehicle vs aPDL1 p-value =",         pvec[2], "\n\n")
  cat("Vehicle vs 1D11+aPDL1 p-value =",    pvec[3], "\n\n")
  cat("1D11 vs aPDL1 p-value =",            pvec[4], "\n\n")
  cat("1D11 vs 1D11+aPDL1 p-value =",       pvec[5], "\n\n")
  cat("aPDL1 vs 1D11+aPDL1 p-value =",      pvec[6], "\n")
}
