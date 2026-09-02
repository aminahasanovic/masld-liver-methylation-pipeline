# DNA methylation and epigenetic age acceleration in MASLD liver tissue

Analysis code for my master's thesis, *Epigenetic Age Acceleration and DNA
Methylation Signatures Associated with Disease Severity in Metabolic
Dysfunction-associated Steatotic Liver Disease*.

Six published liver DNA methylation studies were harmonised into a single
cohort of **621 samples and 377,857 shared CpGs** (450K and EPIC arrays,
640 array profiles before technical replicates were resolved to one profile
per sample) and used for three analyses:

1. epigenetic age acceleration (EAA) across disease stages, using established
   methylation clocks;
2. supervised CpG-based classifiers of histological disease stage, validated
   both by random train/test splits and leave-one-study-out (LOSO);
3. candidate CpGs associated with disease severity, with cross-study effect
   modelling and random-effects meta-analysis.

The repository contains the code and the result tables reported in the thesis.
Methylation matrices are not redistributed here — all source datasets are
public and listed in [`data/source_manifest.csv`](data/source_manifest.csv).

## Repository layout

```
data_preprocessing/   loading, harmonisation, filtering, normalisation, PCA
epigenetic_age/       clock estimation and EAA models
classifier/           supervised models, candidate CpGs, sensitivity analyses
R/                    shared I/O and plotting helpers
workflow/             pipeline entry point
docs/                 pipeline, data and software documentation
thesis_results/       result tables as reported in the thesis
```

## Running the pipeline

Requires R 4.5. Install the packages listed in
[`docs/SOFTWARE.md`](docs/SOFTWARE.md) with

```bash
Rscript workflow/install_dependencies.R
```

download the source data as described in [`docs/DATA.md`](docs/DATA.md), then:

```bash
bash workflow/run_pipeline.sh core
```

`core` runs preprocessing, PCA, epigenetic age, and the primary four-class
Elastic Net under both validation schemes. `full` additionally runs the
Healthy-obese exclusion sensitivity analysis. Individual analyses can be run
script by script; see [`docs/PIPELINE.md`](docs/PIPELINE.md).

Preprocessing is the expensive step (large beta matrices, KNN imputation,
ComBat) and needs a machine with substantial memory; the classifier scripts
run on a laptop once the matrices exist.

## Main findings

EAA results are clock-specific rather than uniform. Classifier accuracy is
high under random train/test splits and drops sharply under LOSO, where the
binary Healthy-vs-Disease contrast is the only reasonably transferable signal
(`thesis_results/tables/T06_classifier_model_comparison.csv`). Disease stage
and source study are partly
confounded in the integrated cohort, which limits how far the stage-related
signal can be separated from study background.

## Data availability

Six of the seven processed datasets come from GEO; the ITEN cohort
(K-BDS, `KAP240571`) may require an access request. The Kurokawa
transfer analysis additionally needs a sample-level age/sex table that is not
redistributed here — see [`docs/DATA.md`](docs/DATA.md).

## License

Code is released under the MIT License (see [`LICENSE`](LICENSE)). The
published datasets and any third-party supplementary tables remain subject to
their original terms.
