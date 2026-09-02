# Software

The analyses were run with **R 4.5.0** on Linux (x86_64) against
**Bioconductor 3.22**. The versions below are those of the environment in
which the thesis results were produced, queried from that environment for
every package the code actually loads.

## Package versions

### Data import, filtering and normalisation

| Package | Version | Used for |
|---|---|---|
| GEOquery | 2.78.0 | downloading GEO series |
| minfi | 1.56.0 | IDAT import, detection p-values |
| ChAMP | 2.40.0 | probe filtering, BMIQ normalisation, KNN imputation, ComBat |
| sva | 3.58.0 | ComBat implementation, called through `ChAMP::champ.runCombat` |
| impute | 1.84.0 | KNN imputation backend used by ChAMP |
| wateRmelon | 2.16.0 | BMIQ backend used by ChAMP |
| readxl | 1.5.0 | Excel metadata from the source publications |

### Epigenetic clocks

| Package | Version | Used for |
|---|---|---|
| methylCIPHER | 0.2.0 | epigenetic clock estimation |
| CTSclocks | 0.0.1 | liver-specific clocks |
| DunedinPoAm38 | 0.1.0 | pace-of-ageing estimator required by methylCIPHER |
| prcPhenoAge | 0.1.0 | principal-component clock required by methylCIPHER |
| HiBED | 1.8.0 | cell-type deconvolution required by CTSclocks |
| qs2 | 0.2.2 | serialisation format of the methylCIPHER reference objects |

### Annotation

| Package | Version | Used for |
|---|---|---|
| IlluminaHumanMethylation450kanno.ilmn12.hg19 | 0.6.1 | 450K probe annotation |
| IlluminaHumanMethylationEPICanno.ilm10b4.hg19 | 0.6.0 | EPIC probe annotation |
| org.Hs.eg.db | 3.22.0 | gene annotation |
| TxDb.Hsapiens.UCSC.hg19.knownGene | 3.22.1 | transcript coordinates for CpG context |
| GenomicFeatures | 1.62.0 | promoter and gene-body regions |
| GenomicRanges / IRanges / S4Vectors | 1.62.1 / 2.44.0 / 0.48.1 | interval arithmetic |
| AnnotationDbi | 1.72.0 | annotation database queries |
| Biobase | 2.70.0 | expression-set containers used by upstream Bioconductor packages |

### Models and statistics

| Package | Version | Used for |
|---|---|---|
| glmnet | 4.1.10 | Elastic Net |
| ranger | 0.18.0 | Random Forest |
| xgboost | 3.2.1.1 | boosted trees |
| caret | 7.0.1 | resampling, tuning |
| yardstick | 1.4.0 | classification metrics |
| pROC | 1.19.0.1 | ROC and AUC |
| limma | 3.66.0 | differential methylation in feature selection |
| broom | 1.0.12 | tidying model output |
| matrixStats | 1.5.0 | row-wise variance for CpG pre-selection |
| Matrix | 1.7.5 | sparse design matrices for glmnet |

### Exploratory structure

| Package | Version | Used for |
|---|---|---|
| uwot | 0.2.4 | UMAP embedding |
| Rtsne | 0.17 | t-SNE embedding |

### Figures

| Package | Version | Used for |
|---|---|---|
| ggplot2 | 4.0.3 | figures |
| patchwork / viridis / GGally | 1.3.2 / 0.6.5 / 2.4.0 | figure composition, palettes, pair plots |
| ggrepel | 0.9.8 | non-overlapping point labels |
| scales | 1.4.0 | axis formatting |
| png | 0.1.9 | raster panels assembled into composite figures |

### Data handling and infrastructure

| Package | Version | Used for |
|---|---|---|
| dplyr / tidyr / readr / tibble / stringr / forcats / purrr | 1.2.1 / 1.3.2 / 2.2.0 / 3.3.1 / 1.6.0 / 1.0.1 / 1.2.2 | data handling |
| janitor | 2.2.1 | column-name cleaning of source metadata |
| rlang | 1.3.0 | tidy evaluation in the helper functions |
| yaml / fs / glue | 2.3.12 / 2.1.0 / 1.8.1 | configuration and paths |
| R.utils | 2.13.0 | gzip handling of downloaded matrices |
| BiocManager / remotes | 1.30.27 / 2.5.0 | package installation |

The random-effects meta-analysis in
`classifier/scripts/11_meta_analyze_candidate_cpgs.R` is implemented directly
(DerSimonian–Laird estimator of τ²), so no meta-analysis package is required.

The random-effects meta-analysis in
`classifier/scripts/11_meta_analyze_candidate_cpgs.R` is implemented directly
(DerSimonian–Laird estimator of τ²), so no meta-analysis package is required.

A machine-readable version including the exact runtime state is in
[`../thesis_results/supplementary/tables/T34_software_versions_runtime.csv`](../thesis_results/supplementary/tables/T34_software_versions_runtime.csv).

## Installation

All packages above are installed by

```bash
Rscript workflow/install_dependencies.R
```

which installs whatever is missing from the current library — CRAN and
Bioconductor packages by name, the four clock packages from GitHub — and then
prints the resulting versions so a fresh environment can be compared against
the table above. It exits with an error if any package is still missing.

The two clock packages are not on CRAN or Bioconductor and are pinned in the
script to the commits used for the thesis results:

| Package | Source | Commit |
|---|---|---|
| methylCIPHER | `HigginsChenLab/methylCIPHER` | `73f2d7d` |
| CTSclocks | `HGT-UwU/CTSclocks` | `7c242cf` |

`methylCIPHER` needs the PC-clock reference object, which is distributed
separately by the package authors and is several hundred megabytes; it is not
included in this repository.
