# Blood-based Gene Expression Biomarkers for Alzheimer's Disease

Identification of candidate blood biomarkers for early Alzheimer's Disease (AD) diagnosis using LASSO regularisation, Elastic Net, and differential expression analysis on publicly available whole-blood transcriptomic data.

---

## Background

Early diagnosis of Alzheimer's Disease typically relies on expensive or invasive procedures such as PET imaging or cerebrospinal fluid (CSF) analysis. Identifying robust gene expression signatures in peripheral blood offers a non-invasive and scalable diagnostic alternative.

This project uses **GSE63060** (Lunnon et al., 2012 — *Molecular Psychiatry*), a publicly available GEO dataset containing whole-blood Illumina microarray profiles from AD patients and age-matched healthy controls.

---

## Methods

| Step | Method | R package |
|---|---|---|
| Data acquisition | GEO direct download | `GEOquery` |
| Normalisation | Log2 + quantile normalisation | `limma` |
| Differential expression | Linear model + eBayes | `limma` |
| Multiple testing correction | FDR (Benjamini-Hochberg) | base R |
| Dimensionality reduction | PCA | base R |
| Biomarker selection | LASSO (α=1) | `glmnet` |
| Comparison | Elastic Net (α=0.5) | `glmnet` |
| Model evaluation | 10-fold cross-validation (AUC) | `glmnet` |
| Visualisation | Volcano plot, PCA, heatmap | `ggplot2`, `pheatmap` |

### Why LASSO?

In transcriptomic studies, the number of features (genes/probes) far exceeds the number of samples — a classic **p >> n** problem. Standard regression overfits in this setting. LASSO addresses this by adding an L1 penalty that shrinks irrelevant coefficients exactly to zero, performing simultaneous **variable selection** and regularisation.

### Why FDR instead of Bonferroni?

Testing 2,000+ probes simultaneously with Bonferroni correction (FWER) would require p < 0.000025 per test — eliminating most true signals. FDR (Benjamini-Hochberg) instead controls the *proportion* of false positives among significant results, retaining statistical power while limiting false discovery rate to 5%.

---

## Dataset

| Field | Value |
|---|---|
| GEO accession | [GSE63060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63060) |
| Platform | Illumina HumanHT-12 V3.0 |
| Tissue | Whole blood |
| Groups | AD patients vs healthy controls |
| Reference | Lunnon et al. (2012), *Molecular Psychiatry* |

Data is downloaded programmatically within the `.Rmd` — no manual download required.

---

## Repository structure

```
blood-biomarkers-alzheimer-lasso/
│
├── analysis.Rmd          # Main analysis notebook (R Markdown)
├── analysis.html         # Rendered output (knit from .Rmd)
└── README.md
```

---

## How to run

### Requirements

```r
install.packages(c("glmnet", "ggplot2", "dplyr", "tidyr", "pheatmap", "caret"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma"))
```

### Execute

Open `analysis.Rmd` in RStudio and click **Knit**, or run from the R console:

```r
rmarkdown::render("analysis.Rmd")
```

The script downloads GSE63060 automatically on first run (requires internet). Results are cached by knitr after the first download.

---

## Key outputs

- **Volcano plot** — differential expression AD vs Control with FDR correction
- **PCA plot** — sample separation by diagnosis in gene expression space  
- **LASSO cross-validation curve** — AUC vs regularisation parameter λ
- **Selected biomarker plot** — probes retained by LASSO with their coefficients
- **Heatmap** — expression of LASSO-selected probes across all samples

---

## Relevance to multi-omics and exposome research

This project demonstrates core skills directly applicable to integrative exposome and multi-omics research:

- Handling **high-dimensional biological data** (p >> n)
- **Regularisation methods** for feature selection in omics contexts
- **Multiple testing correction** in genome-scale analyses
- Working with **population-based cohort data** from public repositories
- Familiarity with the **brain–body axis**: the dataset links peripheral blood signatures to a neurodegenerative phenotype, illustrating the potential of blood-based proxies for brain health

---

## References

Lunnon K et al. (2012). Mitochondrial dysfunction and immune activation are detectable in early Alzheimer's disease blood. *Journal of Alzheimer's Disease*, 30(3), 685–710.

Barrett T et al. (2013). NCBI GEO: archive for functional genomics data sets — update. *Nucleic Acids Research*, 41(D1), D991–D995.

Tibshirani R (1996). Regression shrinkage and selection via the lasso. *Journal of the Royal Statistical Society: Series B*, 58(1), 267–288.

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300.
