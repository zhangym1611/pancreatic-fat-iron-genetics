# To do before public repository submission

- [ ] Confirm the author's permission to publish the source code; select an explicit code licence (e.g. MIT, BSD-3-Clause, or as permitted by all contributors).
- [ ] Check no hidden tokens, passwords, access credentials, private paths or unpublished proprietary source data are present.
- [ ] Verify all original GWAS IDs against downloaded data headers and the GRCh37 build.
- [ ] Verify Step30a T2D `N`, `ncase`, `ncontrol`, and per-variant sample size against the exact used version of GCST006867.
- [ ] Verify the Step25b exploratory threshold timing; change manuscript `prespecified` if unsupported.
- [ ] Confirm Step29b MVMR estimates are used in Tables/Figures (not preliminary Step29 estimates).
- [ ] Re-run the curated Step30a, Step30b together in the original R/LD environment; check 16/16 and 6/6 statuses respectively.
- [ ] Record R `sessionInfo()`, package versions, PLINK and LDSC versions, checksums for all external input files.
- [ ] Test a fresh run with a second user/environment, or document precise prerequisites and run blockers.
- [ ] Supply complete manuscript figures/tables/ESM from the separate final figure package (this repository does not regenerate all layout outputs).
- [ ] Assign repository URL/Zenodo DOI and update manuscript **Code availability** accordingly.
- [ ] Check disclosures about use of AI tooling against actual use in analysis/code/manuscript preparation.
