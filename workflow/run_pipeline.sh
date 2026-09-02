#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-core}"
if [[ "${MODE}" != "core" && "${MODE}" != "full" ]]; then
  echo "Usage: bash workflow/run_pipeline.sh [core|full]" >&2
  exit 2
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"
export PROJECT_ROOT="${PROJECT_ROOT:-${ROOT_DIR}}"

echo "Running preprocessing"
Rscript data_preprocessing/01_preprocessing_liver.R

echo "Running PCA"
Rscript data_preprocessing/02_pca_liver.R

if [[ "${MODE}" == "full" ]]; then
  echo "Running exploratory UMAP/t-SNE embeddings"
  Rscript data_preprocessing/03_umap_tsne_liver.R
fi

echo "Running epigenetic-age analysis"
Rscript epigenetic_age/01_epi_age_calculation.R

if [[ "${MODE}" == "full" ]]; then
  echo "Running BMI availability audit and EAA sensitivity analysis"
  Rscript epigenetic_age/02_bmi_audit_eaa_sensitivity.R

  echo "Assembling reported EAA outputs"
  Rscript epigenetic_age/03_eaa_audit_outputs.R
fi

echo "Running main classifier analyses"
(
  cd classifier

  env CLASSIFIER_OUTCOME=disease_group_4 CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=train_test \
    Rscript scripts/01_prepare_classifier_data.R
  env CLASSIFIER_OUTCOME=disease_group_4 CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=train_test \
    Rscript scripts/03_train_elastic_net.R
  env CLASSIFIER_OUTCOME=disease_group_4 CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=train_test \
    Rscript scripts/06_evaluate_models.R

  env CLASSIFIER_OUTCOME=disease_group_4 CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=loso \
    Rscript scripts/01_prepare_classifier_data.R
  env CLASSIFIER_OUTCOME=disease_group_4 CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=loso \
    Rscript scripts/07_evaluate_loso_elastic_net.R

  if [[ "${MODE}" == "full" ]]; then
    export CLASSIFIER_OUTCOME=disease_group_4
    export CLASSIFIER_COVARIATE_SET=cpg_only
    export CLASSIFIER_CV_STRATEGY=loso

    echo "Summarising the LOSO Elastic Net run"
    Rscript scripts/08_plot_elastic_net_results.R

    echo "Running the Healthy-obese sensitivity analysis"
    bash scripts/run_elastic_net_no_healthy_obese_sensitivity.sh

    echo "Running the candidate CpG chain under the sparse-stability config"
    bash scripts/run_sparse_stability_candidate_chain.sh
  fi
)

echo "Pipeline complete"
