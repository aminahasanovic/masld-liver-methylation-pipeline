#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p logs

run_step() {
  local log_name="$1"
  shift

  echo
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${log_name}"
  "$@" 2>&1 | tee "logs/${log_name}.log"
}

run_train_test() {
  local outcome="$1"

  run_step "prepare_${outcome}_cpg_only_train_test" \
    env CLASSIFIER_OUTCOME="${outcome}" CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=train_test \
    Rscript scripts/01_prepare_classifier_data.R

  run_step "elastic_net_${outcome}_cpg_only_train_test" \
    env CLASSIFIER_OUTCOME="${outcome}" CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=train_test \
    Rscript scripts/03_train_elastic_net.R

  run_step "evaluate_${outcome}_cpg_only_train_test" \
    env CLASSIFIER_OUTCOME="${outcome}" CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=train_test \
    Rscript scripts/06_evaluate_models.R
}

run_loso() {
  local outcome="$1"

  run_step "prepare_${outcome}_cpg_only_loso" \
    env CLASSIFIER_OUTCOME="${outcome}" CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=loso \
    Rscript scripts/01_prepare_classifier_data.R

  run_step "elastic_net_loso_${outcome}_cpg_only_loso" \
    env CLASSIFIER_OUTCOME="${outcome}" CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=loso \
    Rscript scripts/07_evaluate_loso_elastic_net.R
}

for outcome in \
  disease_group_4_no_healthy_obese \
  disease_group_3_no_healthy_obese \
  binary_healthy_vs_disease_no_healthy_obese
do
  run_train_test "${outcome}"
  run_loso "${outcome}"
done



run_step "summarize_no_healthy_obese_sensitivity" \
  Rscript scripts/26_summarize_no_healthy_obese_sensitivity.R
