# Manuscript analysis ↔ script/result map

| Manuscript component | Principal R script(s) | Reference output shipped | Notes |
|---|---|---|---|
| Imaging GWAS QC | Step18 | no checkpoint | eight source traits and QC checks |
| Heritabilities | Step20–21 | no checkpoint | GWAS LDSC pipeline |
| 28 pairwise genetic correlations | Step22 | no checkpoint | genetic correlation/FDR |
| Chromosome-8 exclusion | Step23 | `Step23__full_vs_nochr8_rg_comparison.tsv` | whole-chromosome LDSC sensitivity |
| Genome-wide significant locus overlap | Step24 | no checkpoint | LD clumping and exact SNP overlap |
| Primary local-h2 LAVA screen | Step25 FINAL v3 | no checkpoint | no region passed strict h2 screen for both traits |
| Exploratory LAVA | Step25b | `Step25b__Step25b_LAVA_local_rg_all_tested.tsv` | chromosome 19 exploratory; correction for all 2,495 blocks discussed |
| chr8/chr19 ABF + multi-signal | Step26 FINAL v4 with Step26b fix/recovery | `Step26__Step26_core_summary.tsv` | Step26b comprises repair/recovery, not independent analysis |
| Pancreatic functional annotation / eQTL | Step27 FINAL v2, Step27b FINAL v6 | no checkpoint | MYC, PVT1, PEPD; 3 eQTL datasets |
| Cross-organ, two regions × six pairs | Step28 FINAL v2 | `Step28__Step28_cross_organ_coloc_main.tsv` | main and sensitivity |
| Forward/reverse MR | Step29 FINAL v6 | `Step29_UVMR__Step29_UVMR_results.tsv` | primary forward eight IV after harmonisation |
| FINAL multivariable MR | Step29b GLOBAL CLUMPING | `Step29b_final_MVMR__Step29b_corrected_MVMR_results.tsv` | supersedes the earlier uncorrected MVMR |
| Glycaemic ABF, 16 comparisons | Step30a (renamed from historical Step31) | `Step30a_historical_Step31__Step31_glycemic_coloc_main.tsv` | T2D source sample metadata needs verification |
| Glycaemic SuSiE, six candidate comparisons | Step30b (renamed from historical Step31b) | `Step30b_historical_Step31b__Step31b_integrated_results.tsv` | read only successful Step30a run |
| Manuscript tables/figures | not bundled as an analysis step | separate prepared-figure package | historic graphical Step30 intentionally omitted |

The archived reference outputs are **checkpoints from past analyses**, not proof that the publication scripts have been newly re-executed in the current environment.
