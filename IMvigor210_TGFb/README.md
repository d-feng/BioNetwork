# IMvigor210 TGFβ Analysis

Reproduction of key figures from:

> Mariathasan et al. (2018). *TGF-β attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells.* **Nature**, 554, 544–548.
> https://doi.org/10.1038/nature25501

---

## Overview

This script re-creates **Figure 3** (gene-signature heatmap, Lund subtype inset, extended data heatmaps) and **Figure 4** (mouse cohort Wilcoxon / boxplots) from the IMvigor210 bladder-cancer immunotherapy trial, using the companion R data package `IMvigor210CoreBiologies`.

All plots are saved as PDF files to `output_figures/`.

---

## Prerequisites

### R version
R ≥ 3.5 recommended (tested on R 4.x).

### Required packages

The script auto-installs everything it needs. The main dependencies are:

| Package | Source | Notes |
|---------|--------|-------|
| `IMvigor210CoreBiologies` | Local install (see below) | Provides `cds`, `pData`, mouse data |
| `DESeq` | Built as stub at runtime | Retired package — stub satisfies namespace check for `data(cds)` |
| `Biobase`, `DESeq2`, `edgeR`, `limma` | Bioconductor | Expression analysis |
| `ComplexHeatmap`, `circlize` | Bioconductor | Heatmap rendering |
| `ggplot2`, `dplyr`, `reshape2`, `plyr` | CRAN | Plotting & data wrangling |
| `survival`, `corrplot`, `spatstat`, `DT` | CRAN / Bioconductor | Survival & misc |
| `biomaRt` | Bioconductor | Gene ID mapping |

### Install the data package

Download `IMvigor210CoreBiologies_1.0.1.tar.gz` from the Nature supplementary materials or Genentech/Roche GitHub, then install it locally:

```r
# Remove lsmeans from DESCRIPTION (archived from CRAN) before installing:
# Edit IMvigor210CoreBiologies/DESCRIPTION and remove "lsmeans," from Imports.
install.packages(
  "path/to/IMvigor210CoreBiologies_1.0.1",
  repos = NULL, type = "source"
)
```

---

## Usage

1. **Set `PKG_DIR`** at the top of `IMvigor210_TGFb_analysis.R` to the path of the installed package source (needed to locate the `data/` folder):

   ```r
   PKG_DIR <- "C:/Users/difen/Rcode/IMvigor210CoreBiologies_1.0.1/IMvigor210CoreBiologies"
   ```

2. **Run the script** in its entirety from an R session or RStudio:

   ```r
   source("IMvigor210_TGFb_analysis.R")
   ```

3. **Output** — all figures are written to `output_figures/` (created automatically):

   | File | Description |
   |------|-------------|
   | `Fig3_heatmap.pdf` | Main Figure 3 signature heatmap |
   | `Fig3_inset_response_Lund.pdf` | Response & Lund-subtype inset |
   | `EDFig3b_FGFR3_related.pdf` | Extended Data Fig 3b – FGFR3-related |
   | `EDFig3b_WNT_target.pdf` | Extended Data Fig 3b – WNT target |
   | `EDFig3b_PPARG.pdf` | Extended Data Fig 3b – PPARG |
   | `EDFig3c_immunophenotype_Lund.pdf` | Extended Data Fig 3c |
   | `Fig4_<signature>.pdf` × 8 | Figure 4 mouse Wilcoxon boxplots (one per signature) |

---

## Key Technical Notes

- **`DESeq` stub** — `data(cds)` from `IMvigor210CoreBiologies` serialised the `CountDataSet` S4 class from the retired `DESeq` package (removed from Bioconductor 3.13, 2021). The script builds and installs a minimal stub package named `DESeq` at runtime so that `Biobase::.requirePackage("DESeq")` succeeds without the original package.
- **`lsmeans` removed** — `lsmeans` was archived from CRAN in 2020 and was listed in the data package's `DESCRIPTION` but is not used anywhere in the code. Remove it from `Imports` before installing `IMvigor210CoreBiologies` from source.
- **Plot saving** — `ComplexHeatmap` plots use `pdf()/dev.off()`; all `draw()` and `decorate_heatmap_body()` calls must share the same graphics device. `ggplot2` plots use `ggsave()`.

---

## Repository Structure

```
IMvigor210_TGFb/
├── IMvigor210_TGFb_analysis.R   # Main analysis script
├── output_figures/              # Generated PDFs (git-ignored)
├── README.md
└── .gitignore
```

---

## License / Citation

The analysis code in this repository is adapted from the `IMvigor210CoreBiologies` R package authored by Dorothee Nickles and Richard Bourgon (Genentech/Roche) and released under the package's own license. Please cite the original Nature paper when using these results.
