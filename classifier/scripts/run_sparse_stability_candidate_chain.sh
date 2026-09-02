#!/usr/bin/env bash
# Candidate CpG chain of the thesis.
#
# The candidate CpGs are the Elastic Net features that stay non-zero across
# LOSO folds of several outcome definitions. This selection runs under its own
# configuration (`config/classifier_config_sparse_stability.yaml`), which drops
# ridge from the alpha grid and writes to separate `*_sparse_stability` output
# directories, so it never overwrites the primary four-class results.
#
# The LOSO runs that feed the stability aggregation are listed in that config
# under `stability_selection.eligible_run_ids` and are read from there below,
# so this script cannot drift away from the configuration.

set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p logs_sparse_stability

export CLASSIFIER_CONFIG=config/classifier_config_sparse_stability.yaml

run_step() {
  local log_name="$1"
  shift
  echo
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${log_name}"
  "$@" 2>&1 | tee "logs_sparse_stability/${log_name}.log"
}

mapfile -t run_ids < <(
  Rscript -e 'cat(paste(yaml::read_yaml("config/classifier_config_sparse_stability.yaml")$stability_selection$eligible_run_ids, collapse = "\n"), "\n")'
)

if [[ "${#run_ids[@]}" -eq 0 ]]; then
  echo "No eligible run ids found in the sparse-stability config." >&2
  exit 1
fi

for run_id in "${run_ids[@]}"; do
  run_id="${run_id//[[:space:]]/}"
  [[ -z "${run_id}" ]] && continue

  # Run ids are built as <outcome>_<covariate_set>_<cv_strategy>.
  outcome="${run_id%_cpg_only_loso}"
  if [[ "${outcome}" == "${run_id}" ]]; then
    echo "Run id ${run_id} is not a cpg_only LOSO run; skipping." >&2
    continue
  fi

  run_step "prepare_${run_id}" \
    env CLASSIFIER_OUTCOME="${outcome}" CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=loso \
    Rscript scripts/01_prepare_classifier_data.R

  run_step "elastic_net_loso_${run_id}" \
    env CLASSIFIER_OUTCOME="${outcome}" CLASSIFIER_COVARIATE_SET=cpg_only CLASSIFIER_CV_STRATEGY=loso \
    Rscript scripts/07_evaluate_loso_elastic_net.R
done

# Stability aggregation, annotation, meta-analysis and candidate figures.
export CLASSIFIER_OUTCOME=disease_group_4
export CLASSIFIER_COVARIATE_SET=cpg_only
export CLASSIFIER_CV_STRATEGY=loso

run_step "annotate_candidates"   Rscript scripts/09_annotate_elastic_net_cpgs.R
run_step "plot_candidates"       Rscript scripts/10_plot_candidate_cpgs.R
run_step "meta_analyze"          Rscript scripts/11_meta_analyze_candidate_cpgs.R
run_step "candidate_plots"       Rscript scripts/12_build_candidate_cpg_plots.R
run_step "curation_sheet"        Rscript scripts/13a_prepare_top30_biological_curation.R

curation_file="results/selected_features_sparse_stability/elastic_net_top30_biological_curation_final.csv"
tracked_curation_file="curation/elastic_net_top30_biological_curation_final.csv"
if [[ -f "${curation_file}" || -f "${tracked_curation_file}" ]]; then
  run_step "biological_context"  Rscript scripts/13_add_biological_context_top_cpgs.R
  run_step "thesis_outputs"      Rscript scripts/28_build_final_candidate_thesis_outputs.R
else
  echo
  echo "Manual curation step pending."
  echo "Fill in the review columns of"
  echo "  results/selected_features_sparse_stability/elastic_net_top30_biological_curation_input.csv,"
  echo "save it as elastic_net_top30_biological_curation_final.csv in the same directory,"
  echo "then re-run this script to finish steps 13 and 28."
fi
