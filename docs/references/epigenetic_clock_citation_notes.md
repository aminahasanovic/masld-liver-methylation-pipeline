# Epigenetic clock citations

Maps the clocks used in the EAA analysis to the BibTeX keys in
`epigenetic_clock_citations.bib`. DOIs are listed in
`epigenetic_clock_dois.txt`.

## Clocks reported in the thesis

| Clock column | Citation key |
|---|---|
| `Horvath1` | `Horvath2013DNAmAge` |
| `Horvath2` | `Horvath2018SkinBloodClock` |
| `Hannum` | `Hannum2013GenomeWideMethylation` |
| `PhenoAge` | `Levine2018PhenoAge` |
| `GrimAgeV1` | `Lu2019GrimAge` |
| `PCHorvath1`, `PCHorvath2`, `PCHannum`, `PCPhenoAge`, `PCGrimAge` | `HigginsChen2022PCClocks`, plus the corresponding base-clock paper |
| `HepClock`, `LiverClock` | `Tong2024CTSClocksSubmitted`, plus `Tong2026CTSclocksPackage` for the software |

## Software

| Software | Citation key |
|---|---|
| methylCIPHER | `Thrush2026MethylCIPHER` |
| CTSclocks | `Tong2026CTSclocksPackage` |

## Estimated but not reported

`GrimAgeV2` and HRS-InCH PhenoAge are computed and kept in the clock output
tables, but are not part of the reported EAA results. For `GrimAgeV2`, cite
`Lu2022GrimAgeV2`.

`citeMyClocks("calcHRSInChPhenoAge")` returns no publication, and the
installed methylCIPHER version exposes only the CpG weights for that clock, so
no source is given here.
