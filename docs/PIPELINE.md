# Pipeline

Two stages. The preprocessing/EAA stage runs from the repository root and
writes to `data/processed/` and `outputs/`; the classifier stage runs from
`classifier/` and writes to `classifier/results/`.

`workflow/run_pipeline.sh core` runs preprocessing, PCA, epigenetic age, and
the primary four-class Elastic Net under both validation schemes.
`workflow/run_pipeline.sh full` additionally runs the exploratory UMAP/t-SNE
embeddings, the EAA audit and BMI sensitivity analysis, the candidate CpG chain
up to the final thesis outputs, and the Healthy-obese sensitivity analysis. It
pauses at the one manual step of the candidate chain described below: if the
curated top-30 literature review is not present, the chain stops before script
13 with a message instead of writing incomplete biological support tiers. The
comparison models and the transfer analyses are not part of either mode.
Everything below can also be run script by script.

## Preprocessing and exploratory structure

`data_preprocessing/00_config.R` holds paths, colour and shape mappings, and
the list of clocks. It is sourced by the other scripts and is not run on its
own.

`01_preprocessing_liver.R` loads each dataset (GEO signal intensities or
IDATs via `minfi`), harmonises the phenotype tables into the six-group
`DiseaseGroup` label, applies ChAMP probe filtering and BMIQ normalisation per
dataset, intersects to the CpGs shared by all arrays, resolves technical
replicates, and writes the filtered matrix plus a KNN-imputed and a
ComBat-corrected version. The filtered, non-imputed matrix is the classifier
input; the ComBat version is used for exploratory analyses only.

`02_pca_liver.R` runs PCA on the 50,000 most variable CpGs, raw-imputed
against ComBat-corrected, plus within-study PCA and beta-distribution
summaries.

`03_umap_tsne_liver.R` adds the exploratory nonlinear embeddings: UMAP
(`uwot`) and t-SNE (`Rtsne`) on the first 30 principal components of the same
50,000 CpGs of the ComBat matrix, with a fixed random seed. It writes the
embedding coordinates and the parameter settings alongside the figure. Only
run in `full` mode, since it needs two packages the rest of the pipeline does
not use.

## Epigenetic age

`epigenetic_age/01_epi_age_calculation.R` estimates the clocks via
`methylCIPHER` and `CTSclocks`, computes EAA as the residual of clock age on
chronological age, and fits group models adjusted for sex and dataset with BH
correction applied within each result family. Twelve clocks are reported —
Hannum, Horvath1, Horvath2, PhenoAge, GrimAgeV1, their principal-component
versions (PCHannum, PCHorvath1, PCHorvath2, PCPhenoAge, PCGrimAge) and the two
liver-specific clocks HepClock and LiverClock. The models run on the 565
samples with a recorded chronological age, out of the 621 samples in the final
cohort. Clock references are collected
in [`references/epigenetic_clock_citations.bib`](references/epigenetic_clock_citations.bib).

`02_bmi_audit_eaa_sensitivity.R` documents BMI availability per dataset and
runs the corresponding sensitivity analysis.

`03_eaa_audit_outputs.R` reads the two effect-size tables written by
`01_epi_age_calculation.R` (`DiseaseGroup` and the ordinal `Progression3`
staging variable) and assembles the reported EAA outputs: the combined
contrast table with BH-adjusted q-values, the multiple-testing summary, the
clock probe-coverage and residualisation audits, and the β ± 95 % CI figure.
It checks the expected number of contrasts per staging variable and matrix, so
it fails rather than writing a partial table if an upstream model is missing.

## Classifier

Run from `classifier/`. Three environment variables select the run
configuration; each combination gets its own output directory:

| Variable | Values used in the thesis |
|---|---|
| `CLASSIFIER_OUTCOME` | `disease_group_4` (primary), `disease_group_5_obesity_split`, `disease_group_3`, `binary_healthy_vs_disease`, and the `_no_healthy_obese` variants |
| `CLASSIFIER_COVARIATE_SET` | `cpg_only` (primary), `clinical`, `study_array`, `all_covariates` |
| `CLASSIFIER_CV_STRATEGY` | `train_test`, `loso` |

Outcome definitions and their mapping onto the six harmonised disease groups
are in `config/classifier_config.yaml`.

Primary model, random 70/30 split:

```bash
export CLASSIFIER_OUTCOME=disease_group_4 CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=train_test
Rscript scripts/01_prepare_classifier_data.R
Rscript scripts/03_train_elastic_net.R
Rscript scripts/06_evaluate_models.R
```

Same model, leave-one-study-out:

```bash
export CLASSIFIER_OUTCOME=disease_group_4 CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=loso
Rscript scripts/01_prepare_classifier_data.R
Rscript scripts/07_evaluate_loso_elastic_net.R
```

Feature selection (`02_feature_selection.R`) is sourced by the training
scripts rather than run separately, so it is refitted inside each training set
and inside each LOSO fold.

`08_plot_elastic_net_results.R` summarises the LOSO run: the per-fold and
pooled metric tables, the number of non-zero CpGs per fold, and the confusion
heatmaps reported in the thesis.

Comparison models: `04_train_random_forest.R`, `05_train_boosted_trees.R`,
`16_evaluate_loso_random_forest.R`, and `17_compare_classifier_models.R` for
the combined comparison table.

## Candidate CpGs

The candidate CpGs reported in the thesis are the Elastic Net features that
remain non-zero across LOSO folds of several outcome definitions. This
selection runs under a second configuration,
`config/classifier_config_sparse_stability.yaml`, selected with
`CLASSIFIER_CONFIG`. That configuration drops ridge (alpha = 0) from the alpha
grid, since a ridge fit never sets coefficients to zero and therefore carries
no sparse selection information, and it writes to separate
`*_sparse_stability` output directories so the primary four-class results are
never overwritten. The LOSO runs entering the aggregation are listed in it
under `stability_selection.eligible_run_ids`.

The whole chain — the eligible LOSO runs, the aggregation, and the candidate
figures and tables — is run by

```bash
bash scripts/run_sparse_stability_candidate_chain.sh
```

which reads the eligible run ids from the configuration and stops at the
manual curation step described below. `workflow/run_pipeline.sh full` calls it.
The individual steps, in order:

1. `09_annotate_elastic_net_cpgs.R` — aggregate stability across the eligible LOSO runs and annotate the selected CpGs (Illumina 450K/EPIC manifests)
2. `10_plot_candidate_cpgs.R` — stability across LOSO folds, priority table
3. `11_meta_analyze_candidate_cpgs.R` — per-study effects and random-effects meta-analysis
4. `12_build_candidate_cpg_plots.R` — forest and effect plots
5. `13a_prepare_top30_biological_curation.R` — curation sheet for the top-30 literature review
6. `13_add_biological_context_top_cpgs.R` — gene and regulatory context for the top candidates
7. `28_build_final_candidate_thesis_outputs.R` — final candidate figures and tables

Each script stops with an explicit message if its input from the previous step
is missing.

Between steps 5 and 6 there is one manual step. `13a` writes
`elastic_net_top30_biological_curation_input.csv`, one row per top-30 CpG with
empty columns for the literature review. Those columns were filled in by hand
from the primary literature and the result saved next to it as
`elastic_net_top30_biological_curation_final.csv`, which `13` reads. The
curated file is therefore an input to the pipeline rather than an output of it;
it is the basis of the biological support tiers in main table T08 and of
supplementary table S09.

The curated table used in the thesis is tracked in the repository at
`classifier/curation/elastic_net_top30_biological_curation_final.csv`, so the
chain also completes on a fresh clone. `13` prefers a file of that name in the
results directory and falls back to the tracked copy, which means a repeated
review round only requires saving the new version next to the other candidate
results.

## Sensitivity and transfer analyses

- `21_validate_holdout_samples.R`, `23_summarize_johnson_all_holdout.R`,
  `24_plot_johnson_holdout_validation.R` — Johnson holdout, activated with
  `CLASSIFIER_HOLDOUT=johnson_nonfibrotic` or `johnson_all`
- `22_summarize_obesity_split_sensitivity.R`,
  `26_summarize_no_healthy_obese_sensitivity.R` and
  `run_elastic_net_no_healthy_obese_sensitivity.sh` — Healthy-obese handling
- `25_kurokawa_transfer_analysis.R` — applying the model to the Kurokawa
  cross-etiology cohort (needs the metadata file described in `DATA.md`)

## Runtime

Preprocessing dominates: loading and normalising the six datasets, the
common-CpG intersection, KNN imputation and ComBat over a 377,857 × 621 matrix
need a large-memory machine and several hours. The classifier scripts run in
minutes to tens of minutes on the prepared matrices, LOSO being the slowest
because feature selection is repeated per fold.
