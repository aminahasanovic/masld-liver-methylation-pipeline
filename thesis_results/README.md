# Result tables

Tables as reported in the thesis, copied from the pipeline output. Figures are
not included; they are in the thesis document.

- `tables/` — main tables (T00–T23)
- `supplementary/tables/` — supplementary tables (S04–S09, T24–T43)
- `diagnostics/tables/` — internal validation metrics

## Numbering and provenance

The `F…`/`S…`/`T…` identifiers are a curated layer on top of the pipeline
output: they order the results as they appear in the thesis, and only a subset
of them is assigned inside the analysis scripts themselves. The manifests make
the mapping explicit, so every reported result can be traced back to the code
that produced it:

- `tables/thesis_figure_manifest.csv` — identifier, file stem, description,
  the thesis figure number (`thesis_figure`, empty for outputs that were
  produced but not shown in the thesis), the generating script
  (`source_script`) and that script's own output filename
  (`source_output_file`).
- `tables/thesis_table_manifest.csv` and the two manifests in
  `supplementary/tables/` — the same `source_script` / `source_output_file`
  columns.

Where `source_output_file` reads "curated summary …", the numbered table is a
formatted summary of that script's outputs rather than a file the script
writes under this name.

Rerunning the pipeline writes to `outputs/` and `classifier/results/`, not
here, so these files stay as the reference version of the reported numbers.
