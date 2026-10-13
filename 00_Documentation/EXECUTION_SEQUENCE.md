# Suggested execution and verification order

The code is a **historical modular pipeline** requiring original data and stage-specific saved outputs. Do not execute all files automatically; inspect their configuration first.

**Early validation / history**: Steps 07, 09–17b are early chr8 validation and feasibility work. The reconstructed early steps 01–06/08 are stored separately as provenance, not as verified originals.

**Primary main pipeline**:

1. Step18: audit GWAS summary data.
2. Step19: setup LDSC resources and references.
3. Step20: construct LDSC inputs.
4. Step21: estimate SNP-based heritability.
5. Step22: estimate 28 genetic correlations and FDR.
6. Step23: exclude chromosome 8, recompute genetic correlation.
7. Step24: assess genome-wide locus overlap.
8. Step25 FINAL v3: strict LAVA local h2 / local rg screening.
9. Step25b: exploratory LAVA nominal local h2 sensitivity screen.
10. Step26 FINAL v4: regional colocalisation; Step26b recovery scripts are optional remediation of the historical pipeline and should **not** be rerun without confirming whether the base run needs recovery.
11. Step27 / 27b: functional annotation and PEPD eQTL colocalisation.
12. Step28: cross-organ local architecture.
13. Step29: main univariable MR and historical MVMR.
14. Step29b: corrected **final** MVMR after global union clumping.
15. Step30a: four glycaemic traits × two pancreatic traits × two loci, 16 ABF tests.
16. Step30b: six Step30a candidate tests, SuSiE-RSS and coloc.susie. The script accepts only an upstream 16/16-complete run and can pin that run by env var.

**Important**: Steps rely on historical output folders, external LD reference data and often automatically discover runs. Confirm stage dependencies and the selected input run. For a clean deployment, consider pinned config paths and a version-locked R environment. No automated one-command launch is provided because it would falsely imply a tested end-to-end pipeline.
