# Historical Step31 renamed to submission Step30a/30b

Historical Step30 was a manuscript table/figure assembly step; it has been omitted from the analytical sequence as requested by the author. The subsequent glycaemic analyses were originally Step31 and Step31b. The curated scripts are now:

| Published script | Historical original |
|---|---|
| `Step30a_glycaemic_ABF_colocalisation.R` | `Step30_prioritized_locus_glycemic_followup_v7_ORIENTATION_FIX.R` (internal `STEP 31`) |
| `Step30b_glycaemic_SuSiE_colocalisation.R` | `Step30_glycemic_SuSiE_coloc_COMPLETE.R` (internal `STEP 31b`) |

Both historical original files are preserved unchanged in `02_Provenance/original_glycaemic_scripts/` and identified by SHA256 in the manifest.

In the submission copies: Step31/31b text identifiers, output subdirectories, report filenames and completion markers have been converted consistently to Step30a/30b; the two codes continue to use the same GWAS, variants, priors, LD reference and statistical procedures. The project root has an optional environment override, without changing the original default. The second script now accepts **only** the exact successful Step30a status (16/16 comparisons OK), rather than the old broad prefix matching that also accepted `WITH_REVIEW_FLAGS` runs. The candidate comparisons remain explicitly restricted to six. `PANCREAS_STEP30A_RUN_DIR` optionally pins the exact upstream run; otherwise the newest strictly completed run is selected. This is a workflow-safety change, not an alteration of estimators.

The **historical** reference checkpoints in `03_Reference_checkpoints` retain original Step31 file headers and output names, and are not falsely represented as fresh Step30 re-runs.
