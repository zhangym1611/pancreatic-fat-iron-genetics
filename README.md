# Pancreatic fat and pancreatic iron — analysis code (Diabetologia manuscript)

This repository-style distribution contains the R analysis scripts used in the genome-wide, local, cross-organ, glycaemic-colocalisation and Mendelian randomisation work on MRI-derived pancreatic fat and pancreatic iron. It is a **curated submission candidate**, not a claim of one-click reproducibility or fresh end-to-end validation.

## Directory layout

- `01_Analysis_scripts/`: 29 analytic scripts arranged by scientific stage, including final glycaemic follow-up Step30a and Step30b.
- `02_Provenance/`: the reconstructed early-stage script (not an original historical source) and byte-identical historical Step31/31b scripts.
- `03_Reference_checkpoints/`: small, historical published-analysis output excerpts for comparison; intentionally retains historical filenames.
- `00_Documentation/`: requirements, execution sequence, analysis-to-manuscript mapping, limitations, release checklist, file inventory and SHA-256 checksums.
- `04_Utilities/verify_package.py`: validates file hashes and key curated dependency markers without running R.

## Code naming and analytical provenance

The historical manuscript-output Step30 was omitted by the author. The glycaemic ABF and SuSiE stages previously numbered Step31/31b have been renamed in the curated source to **Step30a/30b**. Only identifiers, output paths, safe upstream-run selection and optional project-root settings were changed, not statistical models. See `00_Documentation/RENAME_AND_PROVENANCE.md` and the originals under `02_Provenance/`.

For independent execution you must supply the publicly available GWAS/eQTL summary data and the European LD/LDSC/LAVA reference resources. Most of the earlier R scripts expect historical input/output paths in `~/胰腺` and tool binaries under `~/mr_tools`. See `00_Documentation/INPUTS_AND_ENVIRONMENT.md` for accession IDs and software requirements, and `00_Documentation/EXECUTION_SEQUENCE.md` for the run order.

## Main result attribution

- **Primary LAVA:** strict screen did not yield testable dual-trait loci; chromosome 19 is exploratory under nominal local heritability thresholds.
- **Colocalisation:** primary chr8 and exploratory chr19 ABF and SuSiE checks use different models; high PP.H4 does not identify a specific effector gene.
- **MR:** Step29 is used for univariable MR, **Step29b is the final MVMR**. Step30a/30b support the glycaemic colocalisation analyses.
- **Study limitations:** overlapping UK Biobank-based data; ancestry restrictions; MRI iron is an imaging estimate; gene–environment equivalence not established by MR. These interpretation caveats do not disappear through code curation.

## Status and integrity verification

Run `python3 04_Utilities/verify_package.py` from the package root to check hashes and the Step30a/30b hand-off. **No R statistical scripts were executed during this packaging task**. The results in `03_Reference_checkpoints` are copies of historical outputs, not newly generated independent replications.

Please read `00_Documentation/KNOWN_LIMITATIONS.md` before depositing this archive in a public repository. Do not write 'fully reproducible' in the paper or code repository before the environment/data test is completed.

## Suggested repository citation

Assign a permanent repository/Zenodo identifier, version number, authors and license after the authors approve the release. No DOI or licence is implied by this draft.
