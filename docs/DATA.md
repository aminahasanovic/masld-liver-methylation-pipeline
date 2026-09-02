# Data

No methylation data is stored in this repository. All datasets used in the
thesis are public; this file describes where they come from and where the
pipeline expects them.

## Source datasets

| Dataset | Accession | Repository | Array | Input | Role |
|---|---|---|---|---|---|
| Ahrens | GSE48325 | GEO | 450K | signal intensities | main cohort |
| Horvath | GSE61258 | GEO | 450K | signal intensities | main cohort |
| Johnson | GSE180474 | GEO | EPIC | signal intensities | main cohort, holdout validation |
| Murphy | GSE49542 | GEO | 450K | signal intensities | main cohort |
| VanDijck | GSE294806 | GEO | EPIC | IDAT | main cohort |
| ITEN | KAP240571 | K-BDS | EPIC | IDAT | main cohort |
| Kurokawa | GSE60753 | GEO | 450K | signal intensities | cross-etiology transfer analysis |

Machine-readable version: [`../data/source_manifest.csv`](../data/source_manifest.csv).

The ITEN cohort is hosted at the Korea BioData Station under study accession
`KAP240571` and may require an access request.

## Expected directory layout

```
data/
  raw/
    Ahrens_GSE48325/
    Horvath_GSE61258/
    Johnson_GSE180474/
    Murphy_GSE49542/
    VanDijck_GSE294806/
    Kim_KSE102917/
      idat/                      212 IDAT files
      BioSample_metadata-2.xlsx  K-BDS BioSample export
      BioSample_idat_map.csv     BioSample name -> Sentrix basename
    Kurokawa_GSE60753/
  processed/                     written by the preprocessing script
  metadata/
  references/
```

GEO series are downloaded with `GEOquery`; paths are configured in
`data_preprocessing/00_config.R`.

## Files not redistributed here

**Kurokawa age and sex table** — `data/metadata/S_table1_27_12_2024.xlsx`,
read with `skip = 13` in `01_preprocessing_liver.R`. Sample-level table with
the columns `Sample`, `Group`, `Etiology`, `Simple Group`, `Age`, `Sex`. It is
not included because it contains individual-level covariates that were not
published alongside the GEO series but directly retrieved from the authors of the paper. Without it, the Kurokawa transfer analysis
(`classifier/scripts/25_kurokawa_transfer_analysis.R`) cannot be reproduced;
the rest of the pipeline is unaffected. In case you want to reproduce this step, reach out to me privately.

## Cohort after preprocessing

621 samples and 377,857 shared CpGs. Sample-level counts per study, disease
group, array and metadata completeness are in
[`../thesis_results/tables/T00_study_overview_final_cohort.csv`](../thesis_results/tables/T00_study_overview_final_cohort.csv)
and `T00E_thesis_dataset_overview.csv`; the exclusion steps from 650
harmonised profiles to the final 621 samples are in
`../thesis_results/supplementary/tables/T25_sample_exclusion_flow_650_to_621.csv`.

Disease labels from the source studies are harmonised into six groups
(`Healthy`, `Healthy obese`, `MASL`, `MASH`, `Mild fibrosis`,
`Advanced fibrosis`); the per-study mapping is documented in
`T00C_disease_label_harmonization_summary.csv`. Coarser outcomes used by the
classifier are derived from these six groups in
`classifier/config/classifier_config.yaml`.
