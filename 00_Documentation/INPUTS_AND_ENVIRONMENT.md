# Inputs and environment

## GWAS accession IDs (GRCh37)

Imaging: GCST90016666 (liver volume), GCST90016669 (pancreatic volume), GCST90016671 (visceral adipose tissue), GCST90016672 (abdominal subcutaneous adipose tissue), GCST90016673 (liver fat), GCST90016674 (liver iron), GCST90016675 (pancreatic fat), GCST90016676 (pancreatic iron). The order follows the project's downloaded-file mapping; verify against the exact source accession records before independently re-running.

Glycaemic: GCST90002232 (fasting glucose), GCST90002238 (fasting insulin), GCST90002244 (HbA1c), GCST006867 (type 2 diabetes). Access via GWAS Catalog/OpenGWAS source links. For T2D, use the **actual version-specific** case/control metadata, not assumed numbers from other Xue 2018 distributions.

eQTL: eQTL Catalogue, GTEx pancreas `QTD000296`, van de Bunt islets `QTD000554`, PISA islets `QTD000574`. The main analysed regions use genome build GRCh37; variant harmonisation follows rsID/allele matching as stated in the code.

## External dependencies

- R and packages used in the scripts including data.table, coloc, susieR, LAVA, GenomicRanges, Rsamtools, IRanges, TwoSampleMR, httr, jsonlite, and helpers. Some stages also call R packages conditionally.
- PLINK 1.9 executable; 1000 Genomes European LD reference panel and required bed/bim/fam files.
- LDSC executable and EUR LD scores, HapMap3 variants, and LAVA European reference blocks.
- Linux/Unix utilities `curl`, `python3`, `zip`; for some stages `tabix`, `gzip`, `wget` or other common command-line tools.
- Network is needed for particular GWAS/eQTL dataset downloads and package installation. Large source data are intentionally excluded from this archive.

## Expected original working directory

Most archived scripts default to `~/胰腺` (for the data and outputs) and to `~/mr_tools` (for the command-line tools). For a faithful rerun, create/provision these paths or carefully update the local path configuration while recording your changes. **Do not** change paths inside a scientific workflow without verifying each downstream artifact.

For the two published glycaemic follow-ups only, you may set:

```sh
export PANCREAS_PROJECT_ROOT="$HOME/胰腺"
export PANCREAS_STEP30A_RUN_DIR="$HOME/胰腺/30a_glycaemic_locus_followup/glycemic_followup_<specific-id>"
```

`PANCREAS_STEP30A_RUN_DIR` is optional and is read only by Step30b. It must point to an output with a fully successful completion marker. It is safer to pin a specific run when multiple attempts exist.

## Restore environment provenance

Before rerunning, record `R.version.string`, `sessionInfo()`, `packageVersion()` for each R package, `plink --version`, `python3 --version`, and checksums for each downloaded GWAS/LD/eQTL file. Historical scripts may reject mismatched R package or PLINK versions (especially the reconstructed early steps). No `renv.lock` is supplied because a verified common R environment is not available in this archive.
