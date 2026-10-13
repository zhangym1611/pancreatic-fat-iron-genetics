# ============================================================
# STEP 31b
# Signal-level glycemic colocalization at prioritized loci
# SuSiE-RSS + coloc.susie
#
# PURPOSE
# -------
# Follow up ONLY the six Step31 ABF candidate comparisons:
#
# chr8 primary locus:
#   1) pancreatic fat  vs T2D
#   2) pancreatic iron vs T2D
#
# chr19 secondary locus:
#   3) pancreatic fat  vs fasting glucose
#   4) pancreatic fat  vs fasting insulin
#   5) pancreatic fat  vs T2D
#   6) pancreatic iron vs fasting insulin
#
# This step addresses multi-signal architecture. A strong ABF PP.H4 alone
# is not sufficient to conclude that the glycemic signal is the same local
# signal already shared by pancreatic fat and pancreatic iron.
#
# INPUT
# -----
# Automatically finds the latest completed Step31 directory under:
#   ~/胰腺/31_glycemic_locus_followup/
#
# Uses Step31 harmonized files and:
#   ~/胰腺/00_ld_reference/EUR
# or
#   ~/胰腺/ld_reference/EUR
#
# PLINK:
#   ~/mr_tools/plink
#   ~/mr_tools/plink/plink
#   or system plink
#
# METHODS
# -------
# - Retain unique rsID SNVs
# - MAF >= 0.01 already enforced in Step31; rechecked here
# - Exact GRCh37 position match to 1000G EUR reference
# - Allele-set match allowing strand complement
# - PLINK reference subset + signed LD matrix
# - Align both pancreatic and glycemic effects to PLINK A1
# - SuSiE-RSS via coloc::runsusie()
# - coloc::coloc.susie()
# - priors: p1=p2=1e-4, p12=1e-5
# - L=10, coverage=0.95
#
# IMPORTANT INTERPRETATION
# ------------------------
# - Signal-pair PP.H4 supports a shared local signal under the model.
# - It does NOT prove that pancreatic fat/iron causes the glycemic trait.
# - chr19 remains an exploratory secondary locus.
# - T2D and pancreatic MRI GWAS may share UK Biobank participants; this
#   does not invalidate locus colocalization but should be noted.
#
# OUTPUT
# ------
# ~/胰腺/31b_glycemic_susie/<run_id>/
#   susie/<pair_id>/
#   reports/
#   31b_glycemic_susie_review.zip
# ============================================================

(function() {

  options(stringsAsFactors = FALSE)

  base <- path.expand("~/胰腺")

  # ------------------------------------------------------------
  # 1. R packages
  # ------------------------------------------------------------
  lib <- path.expand("~/R_pancreas_step31b")
  dir.create(lib, recursive=TRUE, showWarnings=FALSE)
  .libPaths(c(lib, .libPaths()))

  need_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
      install.packages(pkg, lib=lib, repos="https://cloud.r-project.org")
    }
    if (!requireNamespace(pkg, quietly=TRUE)) {
      stop("Cannot install/load R package: ", pkg)
    }
  }

  need_pkg("data.table")
  need_pkg("coloc")
  need_pkg("susieR")

  suppressPackageStartupMessages({
    library(data.table)
    library(coloc)
    library(susieR)
  })

  cat("data.table: ", as.character(packageVersion("data.table")), "\n", sep="")
  cat("coloc:      ", as.character(packageVersion("coloc")), "\n", sep="")
  cat("susieR:     ", as.character(packageVersion("susieR")), "\n", sep="")

  # ------------------------------------------------------------
  # 2. Find latest completed Step31 run
  # ------------------------------------------------------------
  step31_parent <- file.path(base, "31_glycemic_locus_followup")

  if (!dir.exists(step31_parent)) {
    stop("Step31 parent directory not found: ", step31_parent)
  }

  status_files <- list.files(
    step31_parent,
    pattern="completion_status\\.txt$",
    recursive=TRUE,
    full.names=TRUE
  )

  if (!length(status_files)) {
    stop("No Step31 completion_status.txt found.")
  }

  completed <- status_files[
    vapply(status_files, function(f) {
      x <- tryCatch(readLines(f, warn=FALSE), error=function(e) "")
      any(grepl("^COMPLETED_STEP31", x))
    }, logical(1))
  ]

  if (!length(completed)) {
    stop("No completed Step31 run found.")
  }

  # choose newest completed Step31 directory
  completed <- completed[
    order(file.info(completed)$mtime, decreasing=TRUE)
  ]

  step31_dir <- dirname(dirname(completed[1]))

  candidate_file <- file.path(
    step31_dir,
    "reports",
    "Step31_SuSiE_followup_candidates.tsv"
  )

  if (!file.exists(candidate_file)) {
    stop("Step31 candidate file missing: ", candidate_file)
  }

  cand_all <- fread(candidate_file)

  # ------------------------------------------------------------
  # 3. Freeze the six intended comparisons
  # ------------------------------------------------------------
  target <- data.table(
    region_id = c(
      "chr8_primary_shared_locus",
      "chr8_primary_shared_locus",
      "chr19_secondary_shared_locus",
      "chr19_secondary_shared_locus",
      "chr19_secondary_shared_locus",
      "chr19_secondary_shared_locus"
    ),
    pancreatic_trait = c(
      "pancreas_fat",
      "pancreas_iron",
      "pancreas_fat",
      "pancreas_fat",
      "pancreas_fat",
      "pancreas_iron"
    ),
    glycemic_trait = c(
      "T2D",
      "T2D",
      "fasting_glucose",
      "fasting_insulin",
      "T2D",
      "fasting_insulin"
    )
  )

  target[, pair_id := paste(
    region_id,
    pancreatic_trait,
    glycemic_trait,
    sep="__"
  )]

  cand_all[, pair_id := paste(
    region_id,
    pancreatic_trait,
    glycemic_trait,
    sep="__"
  )]

  # Verify that Step31 actually selected each target.
  missing_candidates <- setdiff(target$pair_id, cand_all$pair_id)

  if (length(missing_candidates)) {
    stop(
      "The following prespecified Step31b targets were not present in ",
      "Step31_SuSiE_followup_candidates.tsv:\n",
      paste(missing_candidates, collapse="\n")
    )
  }

  target <- merge(
    target,
    cand_all[, .(
      pair_id,
      ABF_PP.H3 = PP.H3,
      ABF_PP.H4 = PP.H4
    )],
    by="pair_id",
    all.x=TRUE
  )

  # ------------------------------------------------------------
  # 4. Glycemic metadata
  # ------------------------------------------------------------
  meta_file <- file.path(
    step31_dir,
    "reports",
    "Step31_glycemic_GWAS_metadata.tsv"
  )

  if (!file.exists(meta_file)) {
    stop("Step31 glycemic metadata missing.")
  }

  gly_meta <- fread(meta_file)

  # pancreatic sample sizes from the original MRI GWAS
  pan_N <- c(
    pancreas_fat = 25617,
    pancreas_iron = 25617
  )

  # ------------------------------------------------------------
  # 5. PLINK and EUR reference
  # ------------------------------------------------------------
  plink_candidates <- unique(c(
    path.expand("~/mr_tools/plink"),
    path.expand("~/mr_tools/plink/plink"),
    Sys.which("plink")
  ))

  plink_candidates <- plink_candidates[
    nzchar(plink_candidates) &
    file.exists(plink_candidates)
  ]

  if (!length(plink_candidates)) {
    stop(
      "PLINK not found. Checked:\n",
      paste(c(
        path.expand("~/mr_tools/plink"),
        path.expand("~/mr_tools/plink/plink"),
        "system PATH"
      ), collapse="\n")
    )
  }

  plink_bin <- plink_candidates[1]

  ref_candidates <- c(
    file.path(base, "00_ld_reference", "EUR"),
    file.path(base, "ld_reference", "EUR")
  )

  is_plink_ref <- function(pref) {
    all(file.exists(paste0(pref, c(".bed",".bim",".fam"))))
  }

  good_ref <- ref_candidates[
    vapply(ref_candidates, is_plink_ref, logical(1))
  ]

  if (!length(good_ref)) {
    stop(
      "No EUR PLINK reference found. Checked:\n",
      paste(ref_candidates, collapse="\n")
    )
  }

  ref_prefix <- good_ref[1]

  cat("Step31 input directory: ", step31_dir, "\n", sep="")
  cat("PLINK: ", plink_bin, "\n", sep="")
  cat("EUR reference: ", ref_prefix, "\n", sep="")

  # ------------------------------------------------------------
  # 6. Output
  # ------------------------------------------------------------
  parent <- file.path(base, "31b_glycemic_susie")
  dir.create(parent, recursive=TRUE, showWarnings=FALSE)

  out <- tempfile("glycemic_susie_", tmpdir=parent)
  dir.create(out)

  susie_dir <- file.path(out, "susie")
  reports <- file.path(out, "reports")
  logs <- file.path(out, "logs")

  for (d in c(susie_dir, reports, logs)) {
    dir.create(d, recursive=TRUE, showWarnings=FALSE)
  }

  writeLines("INCOMPLETE", file.path(reports, "completion_status.txt"))

  fwrite(
    target,
    file.path(reports, "Step31b_targets.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 7. Utility functions
  # ------------------------------------------------------------
  comp <- function(a) {
    chartr("ACGT","TGCA",toupper(a))
  }

  max_asymmetry <- function(R) {
    max(abs(R - t(R)), na.rm=TRUE)
  }

  safe_get_cs <- function(x) {
    tryCatch(
      susieR::susie_get_cs(
        x,
        coverage=0.95,
        min_abs_corr=0
      ),
      error=function(e) NULL
    )
  }

  count_cs <- function(x) {
    if (is.null(x) || is.null(x$cs)) return(0L)
    length(x$cs)
  }

  credible_set_table <- function(susie_obj, cs_obj, label) {

    if (is.null(cs_obj) || is.null(cs_obj$cs) || !length(cs_obj$cs)) {
      return(data.table())
    }

    outlist <- list()

    for (j in seq_along(cs_obj$cs)) {

      idx <- cs_obj$cs[[j]]
      if (!length(idx)) next

      alpha_col <- j
      if (ncol(susie_obj$alpha) < alpha_col) next

      dt <- data.table(
        phenotype=label,
        credible_set=j,
        SNP=susie_obj$variable_name[idx],
        PIP=susie_obj$pip[idx],
        alpha=susie_obj$alpha[alpha_col, idx]
      )

      setorder(dt, -PIP)
      outlist[[length(outlist)+1]] <- dt
    }

    if (!length(outlist)) return(data.table())
    rbindlist(outlist, fill=TRUE)
  }

  # ------------------------------------------------------------
  # 8. Read EUR BIM once
  # ------------------------------------------------------------
  bim <- fread(
    paste0(ref_prefix, ".bim"),
    col.names=c(
      "REF_CHR","SNP","REF_CM","REF_BP","REF_A1","REF_A2"
    ),
    showProgress=FALSE
  )

  # unique rsIDs only
  bim <- bim[
    !duplicated(SNP) &
    !duplicated(SNP, fromLast=TRUE)
  ]

  setkey(bim, SNP)

  # ------------------------------------------------------------
  # 9. One-pair runner
  # ------------------------------------------------------------
  run_pair <- function(region_id, pan_trait, gly_trait, abf_h3, abf_h4) {

    pair_id <- paste(
      region_id,
      pan_trait,
      gly_trait,
      sep="__"
    )

    cat("\n====================================================\n")
    cat("Step31b: ", pair_id, "\n", sep="")
    cat("====================================================\n")

    pair_dir <- file.path(susie_dir, pair_id)
    dir.create(pair_dir, recursive=TRUE, showWarnings=FALSE)

    infile <- file.path(
      step31_dir,
      "harmonized",
      paste0(pair_id, ".tsv")
    )

    if (!file.exists(infile)) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="HARMONIZED_INPUT_MISSING",
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    d <- fread(infile, showProgress=FALSE)

    required <- c(
      "CHR","BP","SNP",
      "REF_pan","ALT_pan",
      "BETA_pan","SE_pan","MAF_pan",
      "BETA_gly_aligned","SE_gly",
      "EAF_gly_aligned"
    )

    miss <- setdiff(required, names(d))

    if (length(miss)) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status=paste0(
          "INPUT_COLUMNS_MISSING:",
          paste(miss, collapse=",")
        ),
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    # ----------------------------------------------------------
    # 9.1 Basic QC
    # ----------------------------------------------------------
    d <- d[
      is.finite(BETA_pan) &
      is.finite(SE_pan) &
      SE_pan > 0 &
      is.finite(BETA_gly_aligned) &
      is.finite(SE_gly) &
      SE_gly > 0 &
      is.finite(MAF_pan) &
      MAF_pan >= 0.01 &
      MAF_pan <= 0.5
    ]

    d <- d[
      grepl("^rs[0-9]+$", SNP)
    ]

    # unique rsIDs only
    d <- d[
      !duplicated(SNP) &
      !duplicated(SNP, fromLast=TRUE)
    ]

    n_harmonized <- nrow(d)

    if (n_harmonized < 100) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="TOO_FEW_HARMONIZED_RSID_VARIANTS",
        n_harmonized=n_harmonized,
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    # ----------------------------------------------------------
    # 9.2 Merge with EUR reference by rsID and exact GRCh37 pos
    # ----------------------------------------------------------
    z <- merge(
      d,
      bim,
      by="SNP",
      all=FALSE
    )

    z <- z[
      CHR == REF_CHR &
      BP == REF_BP
    ]

    n_position_match <- nrow(z)

    if (n_position_match < 100) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="LOW_REFERENCE_POSITION_COVERAGE",
        n_harmonized=n_harmonized,
        n_position_match=n_position_match,
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    # allele set match allowing complement
    allele_ok <- (
      (z$REF_pan == z$REF_A1 & z$ALT_pan == z$REF_A2) |
      (z$REF_pan == z$REF_A2 & z$ALT_pan == z$REF_A1) |
      (comp(z$REF_pan) == z$REF_A1 & comp(z$ALT_pan) == z$REF_A2) |
      (comp(z$REF_pan) == z$REF_A2 & comp(z$ALT_pan) == z$REF_A1)
    )

    z <- z[allele_ok]

    n_allele_match <- nrow(z)

    if (n_allele_match < 100) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="LOW_REFERENCE_ALLELE_COVERAGE",
        n_harmonized=n_harmonized,
        n_position_match=n_position_match,
        n_allele_match=n_allele_match,
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    # ----------------------------------------------------------
    # 9.3 PLINK subset
    # ----------------------------------------------------------
    extract_file <- file.path(pair_dir, "extract_rsids.txt")

    fwrite(
      z[,.(SNP)],
      extract_file,
      col.names=FALSE
    )

    subpref <- file.path(pair_dir, "ref_subset")
    makebed_log <- file.path(pair_dir, "plink_makebed.log")

    rc1 <- system2(
      plink_bin,
      c(
        "--bfile", shQuote(ref_prefix),
        "--extract", shQuote(extract_file),
        "--make-bed",
        "--keep-allele-order",
        "--allow-no-sex",
        "--out", shQuote(subpref)
      ),
      stdout=makebed_log,
      stderr=makebed_log
    )

    if (!identical(as.integer(rc1),0L)) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="PLINK_MAKEBED_FAILED",
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    sbim <- fread(
      paste0(subpref,".bim"),
      col.names=c("CHR","SNP","CM","BP","A1","A2"),
      showProgress=FALSE
    )

    if (nrow(sbim) < 100 || anyDuplicated(sbim$SNP)) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="FINAL_REFERENCE_VARIANTS_INVALID",
        n_final_reference=nrow(sbim),
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    # ----------------------------------------------------------
    # 9.4 Reorder summary stats exactly to subset BIM
    # ----------------------------------------------------------
    setkey(z, SNP)
    zz <- z[sbim$SNP]

    if (
      nrow(zz) != nrow(sbim) ||
      anyNA(zz$SNP) ||
      !identical(as.character(zz$SNP), as.character(sbim$SNP))
    ) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="FINAL_SNP_ORDER_MATCH_FAILED",
        n_final_reference=nrow(sbim),
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    # ----------------------------------------------------------
    # 9.5 Signed LD matrix
    # ----------------------------------------------------------
    ldpref <- file.path(pair_dir, "ld")
    ld_log <- file.path(pair_dir, "plink_ld.log")

    rc2 <- system2(
      plink_bin,
      c(
        "--bfile", shQuote(subpref),
        "--r", "square", "gz",
        "--allow-no-sex",
        "--out", shQuote(ldpref)
      ),
      stdout=ld_log,
      stderr=ld_log
    )

    ldfile <- paste0(ldpref, ".ld.gz")

    if (!identical(as.integer(rc2),0L) || !file.exists(ldfile)) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="PLINK_LD_FAILED",
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    R <- as.matrix(
      fread(
        ldfile,
        header=FALSE,
        showProgress=FALSE
      )
    )

    storage.mode(R) <- "double"

    if (
      nrow(R) != nrow(sbim) ||
      ncol(R) != nrow(sbim) ||
      any(!is.finite(R))
    ) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="LD_MATRIX_INVALID",
        ld_nrow=nrow(R),
        ld_ncol=ncol(R),
        bim_n=nrow(sbim),
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    R <- (R + t(R))/2
    diag(R) <- 1

    rownames(R) <- as.character(sbim$SNP)
    colnames(R) <- as.character(sbim$SNP)

    asym <- max_asymmetry(R)

    if (!is.finite(asym) || asym > 1e-8) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="LD_NOT_SYMMETRIC",
        max_asymmetry=asym,
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    # ----------------------------------------------------------
    # 9.6 Align BOTH beta vectors to PLINK A1
    #
    # Step31 BETA_gly_aligned has already been aligned to pancreatic
    # ALT_pan. Therefore the same sign factor can be applied to both
    # pancreatic and glycemic betas.
    # ----------------------------------------------------------
    alt <- toupper(zz$ALT_pan)
    A1 <- toupper(sbim$A1)
    A2 <- toupper(sbim$A2)

    alt_to_A1 <- (
      alt == A1 |
      comp(alt) == A1
    )

    alt_to_A2 <- (
      alt == A2 |
      comp(alt) == A2
    )

    # ambiguous palindromic orientation relative to LD reference
    ambiguous_align <- alt_to_A1 & alt_to_A2

    if (any(ambiguous_align)) {

      keep <- !ambiguous_align

      zz <- zz[keep]
      sbim <- sbim[keep]
      R <- R[keep,keep,drop=FALSE]

      alt_to_A1 <- alt_to_A1[keep]
      alt_to_A2 <- alt_to_A2[keep]

      rownames(R) <- as.character(sbim$SNP)
      colnames(R) <- as.character(sbim$SNP)
    }

    if (nrow(sbim) < 100) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="TOO_FEW_AFTER_LD_ALLELE_ALIGNMENT",
        n_final=nrow(sbim),
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    bad_align <- !(alt_to_A1 | alt_to_A2)

    if (any(bad_align)) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="ALLELE_TO_LD_ALIGNMENT_FAILED",
        n_bad=sum(bad_align),
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    sign_to_A1 <- ifelse(alt_to_A1, 1, -1)

    beta_pan_ld <- zz$BETA_pan * sign_to_A1
    beta_gly_ld <- zz$BETA_gly_aligned * sign_to_A1

    varbeta_pan <- zz$SE_pan^2
    varbeta_gly <- zz$SE_gly^2

    # glycemic MAF: prefer glycemic EAF, fallback to pancreatic MAF
    gly_maf <- ifelse(
      is.finite(zz$EAF_gly_aligned),
      pmin(
        zz$EAF_gly_aligned,
        1-zz$EAF_gly_aligned
      ),
      zz$MAF_pan
    )

    if (
      any(!is.finite(beta_pan_ld)) ||
      any(!is.finite(beta_gly_ld)) ||
      any(!is.finite(varbeta_pan)) ||
      any(!is.finite(varbeta_gly)) ||
      any(varbeta_pan <= 0) ||
      any(varbeta_gly <= 0) ||
      any(!is.finite(gly_maf)) ||
      any(gly_maf <= 0 | gly_maf > 0.5)
    ) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="NONFINITE_SUMMARY_STATS_AFTER_LD_ALIGNMENT",
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    final_snps <- as.character(sbim$SNP)

    # save final alignment QC
    qc <- data.table(
      order=seq_along(final_snps),
      SNP=final_snps,
      CHR=sbim$CHR,
      BP=sbim$BP,
      LD_A1=sbim$A1,
      LD_A2=sbim$A2,
      pancreatic_REF=zz$REF_pan,
      pancreatic_ALT=zz$ALT_pan,
      beta_pancreatic_A1=beta_pan_ld,
      beta_glycemic_A1=beta_gly_ld,
      varbeta_pancreatic=varbeta_pan,
      varbeta_glycemic=varbeta_gly,
      MAF_pancreatic=zz$MAF_pan,
      MAF_glycemic=gly_maf
    )

    fwrite(
      qc,
      file.path(pair_dir, "final_SNP_order_and_alignment.tsv.gz"),
      sep="\t"
    )

    # ----------------------------------------------------------
    # 9.7 Build SuSiE datasets
    # ----------------------------------------------------------
    D1 <- list(
      beta=beta_pan_ld,
      varbeta=varbeta_pan,
      snp=final_snps,
      position=as.numeric(sbim$BP),
      N=as.numeric(pan_N[[pan_trait]]),
      MAF=zz$MAF_pan,
      type="quant",
      sdY=1,
      LD=R
    )

    gm <- gly_meta[trait==gly_trait]

    if (nrow(gm) != 1) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="GLYCEMIC_METADATA_NOT_UNIQUE",
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    if (gm$type[1] == "quant") {

      # Do not force sdY=1 for MAGIC glycemic traits.
      # coloc::runsusie/check_dataset estimates sdY from N, MAF,
      # beta and varbeta when needed.
      D2 <- list(
        beta=beta_gly_ld,
        varbeta=varbeta_gly,
        snp=final_snps,
        position=as.numeric(sbim$BP),
        N=as.numeric(gm$N[1]),
        MAF=gly_maf,
        type="quant",
        LD=R
      )

    } else {

      D2 <- list(
        beta=beta_gly_ld,
        varbeta=varbeta_gly,
        snp=final_snps,
        position=as.numeric(sbim$BP),
        N=as.numeric(gm$N[1]),
        s=as.numeric(gm$s[1]),
        MAF=gly_maf,
        type="cc",
        LD=R
      )
    }

    if (
      !identical(D1$snp, rownames(D1$LD)) ||
      !identical(D2$snp, rownames(D2$LD))
    ) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status="DATASET_SNP_LD_NAME_MISMATCH",
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    # ----------------------------------------------------------
    # 9.8 SuSiE-RSS
    # ----------------------------------------------------------
    cat("Running pancreatic SuSiE-RSS ...\n")

    s_pan <- tryCatch(
      coloc::runsusie(
        D1,
        L=10,
        maxit=10000,
        coverage=0.95
      ),
      error=function(e) e
    )

    if (inherits(s_pan,"error")) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status=paste0(
          "RUNSUSIE_PAN_ERROR:",
          conditionMessage(s_pan)
        ),
        n_final=length(final_snps),
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    cat("Running glycemic SuSiE-RSS ...\n")

    s_gly <- tryCatch(
      coloc::runsusie(
        D2,
        L=10,
        maxit=10000,
        coverage=0.95
      ),
      error=function(e) e
    )

    if (inherits(s_gly,"error")) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status=paste0(
          "RUNSUSIE_GLY_ERROR:",
          conditionMessage(s_gly)
        ),
        n_final=length(final_snps),
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    cs_pan <- safe_get_cs(s_pan)
    cs_gly <- safe_get_cs(s_gly)

    n_cs_pan <- count_cs(cs_pan)
    n_cs_gly <- count_cs(cs_gly)

    # save CS tables
    cs_pan_tab <- credible_set_table(
      s_pan,
      cs_pan,
      pan_trait
    )

    cs_gly_tab <- credible_set_table(
      s_gly,
      cs_gly,
      gly_trait
    )

    if (nrow(cs_pan_tab)) {
      fwrite(
        cs_pan_tab,
        file.path(pair_dir, "credible_sets_pancreatic.tsv"),
        sep="\t"
      )
    }

    if (nrow(cs_gly_tab)) {
      fwrite(
        cs_gly_tab,
        file.path(pair_dir, "credible_sets_glycemic.tsv"),
        sep="\t"
      )
    }

    # ----------------------------------------------------------
    # 9.9 coloc.susie
    # ----------------------------------------------------------
    cat("Running coloc.susie ...\n")

    csu <- tryCatch(
      coloc::coloc.susie(
        s_pan,
        s_gly,
        p1=1e-4,
        p2=1e-4,
        p12=1e-5
      ),
      error=function(e) e
    )

    if (inherits(csu,"error")) {
      return(data.table(
        pair_id=pair_id,
        region_id=region_id,
        pancreatic_trait=pan_trait,
        glycemic_trait=gly_trait,
        status=paste0(
          "COLOC_SUSIE_ERROR:",
          conditionMessage(csu)
        ),
        n_final=length(final_snps),
        n_cs_pancreatic=n_cs_pan,
        n_cs_glycemic=n_cs_gly,
        ABF_PP.H3=abf_h3,
        ABF_PP.H4=abf_h4
      ))
    }

    # save full objects
    saveRDS(
      list(
        pancreatic=s_pan,
        glycemic=s_gly,
        pancreatic_cs=cs_pan,
        glycemic_cs=cs_gly,
        coloc_susie=csu
      ),
      file.path(pair_dir, "Step31b_susie_objects.rds")
    )

    # ----------------------------------------------------------
    # 9.10 Signal-pair summary
    # ----------------------------------------------------------
    signal_summary <- data.table()

    if (
      !is.null(csu$summary) &&
      nrow(as.data.frame(csu$summary)) > 0
    ) {

      signal_summary <- as.data.table(csu$summary)

      fwrite(
        signal_summary,
        file.path(pair_dir, "coloc_susie_summary.tsv"),
        sep="\t"
      )
    } else {

      fwrite(
        data.table(
          note="coloc.susie completed but returned no signal-pair summary rows"
        ),
        file.path(pair_dir, "coloc_susie_summary.tsv"),
        sep="\t"
      )
    }

    # SNP-level results if present
    if (!is.null(csu$results)) {
      try(
        fwrite(
          as.data.table(csu$results),
          file.path(pair_dir, "coloc_susie_SNP_results.tsv.gz"),
          sep="\t"
        ),
        silent=TRUE
      )
    }

    # ----------------------------------------------------------
    # 9.11 Extract maximum H4 and top shared signal pair
    # ----------------------------------------------------------
    max_h4 <- NA_real_
    top_hit1 <- NA_character_
    top_hit2 <- NA_character_
    top_h3 <- NA_real_

    if (nrow(signal_summary)) {

      h4col <- intersect(
        c("PP.H4.abf","PP.H4"),
        names(signal_summary)
      )

      h3col <- intersect(
        c("PP.H3.abf","PP.H3"),
        names(signal_summary)
      )

      if (length(h4col)) {

        hv <- as.numeric(signal_summary[[h4col[1]]])

        if (any(is.finite(hv))) {

          imax <- which.max(hv)
          max_h4 <- hv[imax]

          if ("hit1" %in% names(signal_summary)) {
            top_hit1 <- as.character(signal_summary$hit1[imax])
          }

          if ("hit2" %in% names(signal_summary)) {
            top_hit2 <- as.character(signal_summary$hit2[imax])
          }

          if (length(h3col)) {
            top_h3 <- as.numeric(
              signal_summary[[h3col[1]]][imax]
            )
          }
        }
      }
    }

    # classify signal-level evidence
    susie_label <- ifelse(
      is.na(max_h4),
      "NO_SIGNAL_PAIR_RESULT",
      ifelse(
        max_h4 >= 0.95,
        "VERY_STRONG_SIGNAL_LEVEL_SHARING",
        ifelse(
          max_h4 >= 0.80,
          "STRONG_SIGNAL_LEVEL_SHARING",
          ifelse(
            max_h4 >= 0.50,
            "MODERATE_SIGNAL_LEVEL_SHARING",
            "LITTLE_SIGNAL_LEVEL_SHARING"
          )
        )
      )
    )

    data.table(
      pair_id=pair_id,
      region_id=region_id,
      pancreatic_trait=pan_trait,
      glycemic_trait=gly_trait,
      status="OK",
      n_step31_harmonized=n_harmonized,
      n_reference_position_match=n_position_match,
      n_reference_allele_match=n_allele_match,
      n_final_LD=length(final_snps),
      n_cs_pancreatic=n_cs_pan,
      n_cs_glycemic=n_cs_gly,
      n_coloc_susie_signal_pairs=nrow(signal_summary),
      max_PP.H4_susie=max_h4,
      top_pancreatic_signal=top_hit1,
      top_glycemic_signal=top_hit2,
      top_PP.H3_susie=top_h3,
      ABF_PP.H3=abf_h3,
      ABF_PP.H4=abf_h4,
      signal_level_interpretation=susie_label
    )
  }

  # ------------------------------------------------------------
  # 10. Run six candidate pairs
  # ------------------------------------------------------------
  status_list <- vector("list", nrow(target))

  for (i in seq_len(nrow(target))) {

    status_list[[i]] <- tryCatch(
      run_pair(
        region_id=target$region_id[i],
        pan_trait=target$pancreatic_trait[i],
        gly_trait=target$glycemic_trait[i],
        abf_h3=target$ABF_PP.H3[i],
        abf_h4=target$ABF_PP.H4[i]
      ),
      error=function(e) {
        data.table(
          pair_id=target$pair_id[i],
          region_id=target$region_id[i],
          pancreatic_trait=target$pancreatic_trait[i],
          glycemic_trait=target$glycemic_trait[i],
          status=paste0(
            "R_EXCEPTION:",
            conditionMessage(e)
          ),
          ABF_PP.H3=target$ABF_PP.H3[i],
          ABF_PP.H4=target$ABF_PP.H4[i]
        )
      }
    )
  }

  status <- rbindlist(status_list, fill=TRUE)

  fwrite(
    status,
    file.path(reports, "Step31b_susie_status.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 11. Build compact cross-pair results table
  # ------------------------------------------------------------
  result <- copy(status)

  if (!"max_PP.H4_susie" %in% names(result)) {
    result[, max_PP.H4_susie := NA_real_]
  }

  result[, passes_signal_level_0.80 := (
    status=="OK" &
    is.finite(max_PP.H4_susie) &
    max_PP.H4_susie >= 0.80
  )]

  result[, passes_signal_level_0.95 := (
    status=="OK" &
    is.finite(max_PP.H4_susie) &
    max_PP.H4_susie >= 0.95
  )]

  fwrite(
    result,
    file.path(reports, "Step31b_integrated_results.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 12. Key-pattern summary
  # ------------------------------------------------------------
  key_summary <- data.table(
    metric=c(
      "n_target_pairs",
      "n_completed_OK",
      "n_max_susie_H4_ge_0.80",
      "n_max_susie_H4_ge_0.95",
      "chr8_T2D_pairs_completed",
      "chr19_fasting_insulin_pairs_completed"
    ),
    value=c(
      nrow(target),
      sum(status$status=="OK"),
      sum(
        status$status=="OK" &
        is.finite(status$max_PP.H4_susie) &
        status$max_PP.H4_susie>=0.80,
        na.rm=TRUE
      ),
      sum(
        status$status=="OK" &
        is.finite(status$max_PP.H4_susie) &
        status$max_PP.H4_susie>=0.95,
        na.rm=TRUE
      ),
      sum(
        status$status=="OK" &
        status$region_id=="chr8_primary_shared_locus" &
        status$glycemic_trait=="T2D"
      ),
      sum(
        status$status=="OK" &
        status$region_id=="chr19_secondary_shared_locus" &
        status$glycemic_trait=="fasting_insulin"
      )
    )
  )

  fwrite(
    key_summary,
    file.path(reports, "Step31b_summary.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 13. Interpretation guardrails
  # ------------------------------------------------------------
  writeLines(
    c(
      "Step31b performs signal-level SuSiE-RSS + coloc.susie only for the six preselected Step31 ABF candidates.",
      "The primary question is whether a glycemic association maps to a specific local signal rather than merely showing region-wide ABF colocalization.",
      "A high signal-pair PP.H4 supports shared local genetic signal architecture under the model; it does not establish causal direction.",
      "The chromosome 8 T2D results should be interpreted in the context of the previously established pancreatic fat–iron multi-signal locus.",
      "The chromosome 19 fasting-insulin results should be interpreted together with the pancreatic fat–VAT sharing identified in Step28.",
      "Chromosome 19 remains exploratory because it was identified after nominal local-heritability screening in Step25b.",
      "Fasting glucose and fasting insulin MAGIC GWAS are BMI-adjusted.",
      "HbA1c was not carried into Step31b because its Step31 ABF evidence was weak.",
      "No new glycemic traits are added in Step31b."
    ),
    file.path(reports, "Step31b_interpretation_guardrails.txt")
  )

  # ------------------------------------------------------------
  # 14. Completion status
  # ------------------------------------------------------------
  n_ok <- sum(status$status=="OK")

  if (n_ok == 6) {
    completion <- "COMPLETED_STEP31B_SUSIE_GLYCEMIC_COLOCALIZATION"
  } else {
    completion <- "COMPLETED_STEP31B_WITH_REVIEW_FLAGS"
  }

  writeLines(
    completion,
    file.path(reports, "completion_status.txt")
  )

  # ------------------------------------------------------------
  # 15. Small review ZIP
  # Exclude LD matrices, BED files and large RDS objects.
  # ------------------------------------------------------------
  review_zip <- file.path(
    out,
    "31b_glycemic_susie_review.zip"
  )

  oldwd <- getwd()
  on.exit(setwd(oldwd), add=TRUE)
  setwd(out)

  review_files <- c(
    list.files(
      "reports",
      recursive=TRUE,
      full.names=TRUE
    ),
    list.files(
      "susie",
      pattern="coloc_susie_summary\\.tsv$|credible_sets_.*\\.tsv$|final_SNP_order_and_alignment\\.tsv\\.gz$",
      recursive=TRUE,
      full.names=TRUE
    )
  )

  review_files <- review_files[
    file.exists(review_files) &
    !dir.exists(review_files)
  ]

  utils::zip(
    review_zip,
    files=review_files,
    flags="-q"
  )

  # ------------------------------------------------------------
  # 16. Console summary
  # ------------------------------------------------------------
  cat("\n====================================================\n")
  cat("STEP31b finished\n")
  cat("====================================================\n\n")

  print(
    status[, .(
      region_id,
      pancreatic_trait,
      glycemic_trait,
      status,
      ABF_PP.H4,
      n_cs_pancreatic,
      n_cs_glycemic,
      n_coloc_susie_signal_pairs,
      max_PP.H4_susie,
      top_pancreatic_signal,
      top_glycemic_signal,
      signal_level_interpretation
    )]
  )

  cat("\nSummary:\n")
  print(key_summary)

  cat("\nCompletion status: ", completion, "\n", sep="")
  cat("Output directory:\n", out, "\n", sep="")
  cat("Review ZIP:\n", review_zip, "\n", sep="")

})()
