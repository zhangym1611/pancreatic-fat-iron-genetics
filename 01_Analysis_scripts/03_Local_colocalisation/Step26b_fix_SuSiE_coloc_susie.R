# 第26b步：修复并完成chr19与chr8的SuSiE-RSS + coloc.susie
# pancreas fat vs pancreas iron
#
# 背景：
# 第26步ABF已完成；SuSiE阶段因LD矩阵缺少rownames/colnames而失败：
#   RUNSUSIE_FAT_ERROR: LD rownames != colnames
#
# 本步只修复并重跑：
#   1) chr19 LOC2340 primary region
#   2) chr8 positive-control region
#
# 不重复：
#   - 原始VCF全区域提取
#   - ABF coloc主分析及敏感性分析
#
# 核心修复：
#   - 以最终PLINK subset BIM顺序为唯一SNP顺序
#   - 设置 rownames(LD) <- sbim$SNP
#   - 设置 colnames(LD) <- sbim$SNP
#   - summary statistics严格重排到同一SNP顺序
#   - effect allele严格对齐到PLINK A1
#   - 检查LD维度、对称性、有限值、SNP命名一致性
#
# 解释边界：
#   coloc.susie用于允许区域内存在多个因果信号时的共定位分析。
#   PP.H4支持共享局部信号模型，但不表示胰腺脂肪导致胰腺铁或相反。

(function() {

  options(stringsAsFactors = FALSE)

  # ============================================================
  # 0. 基础路径与依赖
  # ============================================================
  base <- path.expand("~/胰腺")

  cran_lib <- path.expand("~/R_pancreas_step26")
  dir.create(cran_lib, recursive=TRUE, showWarnings=FALSE)
  .libPaths(c(cran_lib, .libPaths()))

  need_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
      install.packages(pkg, lib=cran_lib, repos="https://cloud.r-project.org")
    }
    if (!requireNamespace(pkg, quietly=TRUE)) {
      stop("无法加载R包：", pkg)
    }
  }

  need_pkg("data.table")
  need_pkg("coloc")
  need_pkg("susieR")

  suppressPackageStartupMessages(library(data.table))

  cat("coloc version: ", as.character(packageVersion("coloc")), "\n", sep="")
  cat("susieR version: ", as.character(packageVersion("susieR")), "\n", sep="")

  # ============================================================
  # 1. 自动定位最新的第26步结果目录
  # ============================================================
  step26_parent <- file.path(base, "26_cross_trait_colocalization")

  if (!dir.exists(step26_parent)) {
    stop("未找到第26步目录：", step26_parent)
  }

  dirs26 <- list.dirs(
    step26_parent,
    recursive=FALSE,
    full.names=TRUE
  )

  if (!length(dirs26)) {
    stop("第26步目录下没有结果批次。")
  }

  # 优先要求存在第26步核心ABF输出和区域harmonized文件
  valid26 <- dirs26[vapply(dirs26, function(d) {

    core1 <- file.path(
      d,
      "coloc_abf",
      "Step26_chr19_primary_coloc.tsv"
    )

    core2 <- file.path(
      d,
      "coloc_abf",
      "Step26_chr8_positive_control_coloc.tsv"
    )

    h19 <- file.path(
      d,
      "regional_sumstats",
      "harmonized__chr19_LOC2340_primary.tsv.gz"
    )

    h8 <- file.path(
      d,
      "regional_sumstats",
      "harmonized__chr8_positive_control.tsv.gz"
    )

    all(file.exists(c(core1, core2, h19, h8)))

  }, logical(1))]

  if (!length(valid26)) {
    stop(
      "没有找到同时包含ABF结果与harmonized区域数据的第26步批次。"
    )
  }

  step26 <- valid26[
    which.max(file.info(valid26)$mtime)
  ]

  cat("\n使用第26步结果目录：\n", step26, "\n")

  # ============================================================
  # 2. 读取第26步ABF主结果，作为26b背景记录
  # ============================================================
  abf_chr19_file <- file.path(
    step26,
    "coloc_abf",
    "Step26_chr19_primary_coloc.tsv"
  )

  abf_chr8_file <- file.path(
    step26,
    "coloc_abf",
    "Step26_chr8_positive_control_coloc.tsv"
  )

  abf_chr19 <- fread(abf_chr19_file)
  abf_chr8 <- fread(abf_chr8_file)

  # ============================================================
  # 3. 输出目录
  # ============================================================
  parent <- file.path(base, "26b_coloc_susie_fix")
  dir.create(parent, recursive=TRUE, showWarnings=FALSE)

  out <- tempfile(
    "fat_iron_susie_fix_",
    tmpdir=parent
  )

  if (!dir.create(out)) {
    stop("无法创建第26b输出目录：", out)
  }

  reports <- file.path(out, "reports")
  susiedir <- file.path(out, "susie")
  logs <- file.path(out, "logs")

  dir.create(reports)
  dir.create(susiedir)
  dir.create(logs)

  writeLines(
    "INCOMPLETE_CHECK_LOGS",
    file.path(reports, "completion_status.txt")
  )

  writeLines(
    c(
      paste0("Step26_source_dir\t", step26),
      "Step26b_scope\tSuSiE-RSS and coloc.susie only",
      "Fix\tLD rownames/colnames explicitly set to final PLINK subset BIM SNP order",
      "Trait1\tpancreas_fat",
      "Trait2\tpancreas_iron",
      "Primary_p1\t1e-4",
      "Primary_p2\t1e-4",
      "Primary_p12\t1e-5",
      "sdY\t1",
      "L\t10"
    ),
    file.path(reports, "analysis_context.txt")
  )

  # ============================================================
  # 4. 固定研究参数
  # ============================================================
  N <- c(
    pancreas_fat = 25617,
    pancreas_iron = 25617
  )

  target_regions <- data.table(
    region_id=c(
      "chr19_LOC2340_primary",
      "chr8_positive_control"
    ),
    role=c(
      "STEP25B_EXPLORATORY_CANDIDATE",
      "POSITIVE_CONTROL"
    )
  )

  # ============================================================
  # 5. PLINK / 1000G EUR reference
  # ============================================================
  plink_candidates <- c(
    path.expand("~/mr_tools/plink"),
    path.expand("~/mr_tools/plink/plink"),
    Sys.which("plink")
  )

  plink_candidates <- unique(
    plink_candidates[nzchar(plink_candidates)]
  )

  plink_bin <- plink_candidates[
    file.exists(plink_candidates)
  ][1]

  if (length(plink_bin)==0 || is.na(plink_bin)) {
    stop("未找到PLINK可执行文件。")
  }

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
      "未找到可用1000G EUR PLINK reference。\n检查路径：\n",
      paste(ref_candidates, collapse="\n")
    )
  }

  ref_prefix <- good_ref[1]

  cat("PLINK: ", plink_bin, "\n", sep="")
  cat("EUR reference prefix: ", ref_prefix, "\n", sep="")

  # ============================================================
  # 6. 工具函数
  # ============================================================
  comp <- function(a) {
    chartr("ACGT", "TGCA", toupper(a))
  }

  safe_cor_symmetry <- function(R) {
    max(abs(R - t(R)), na.rm=TRUE)
  }

  extract_pp_cols <- function(x) {
    nms <- names(x)
    keep <- grep(
      "^PP\\.H[0-4]",
      nms,
      value=TRUE
    )
    keep
  }

  # ============================================================
  # 7. 单区域SuSiE + coloc.susie
  # ============================================================
  run_region <- function(region_id) {

    cat("\n========================================\n")
    cat("处理区域：", region_id, "\n")
    cat("========================================\n")

    region_dir <- file.path(susiedir, region_id)
    dir.create(region_dir, recursive=TRUE, showWarnings=FALSE)

    harm_file <- file.path(
      step26,
      "regional_sumstats",
      paste0("harmonized__", region_id, ".tsv.gz")
    )

    if (!file.exists(harm_file)) {
      return(list(
        status="HARMONIZED_FILE_MISSING",
        region_id=region_id
      ))
    }

    d <- fread(harm_file)

    required <- c(
      "CHR","BP",
      "SNP_fat","SNP_iron",
      "REF_fat","ALT_fat",
      "BETA_fat","VARBETA_fat","MAF_fat",
      "BETA_iron_aligned","VARBETA_iron","MAF_iron_aligned"
    )

    missing_cols <- setdiff(required, names(d))

    if (length(missing_cols)) {
      return(list(
        status=paste0(
          "HARMONIZED_COLUMNS_MISSING:",
          paste(missing_cols, collapse=",")
        ),
        region_id=region_id
      ))
    }

    # ----------------------------------------------------------
    # 7.1 与第26步主ABF一致：MAF >= 0.01
    # ----------------------------------------------------------
    d <- d[
      is.finite(BETA_fat) &
      is.finite(VARBETA_fat) &
      VARBETA_fat > 0 &
      is.finite(BETA_iron_aligned) &
      is.finite(VARBETA_iron) &
      VARBETA_iron > 0 &
      is.finite(MAF_fat) &
      is.finite(MAF_iron_aligned) &
      MAF_fat >= 0.01 &
      MAF_iron_aligned >= 0.01
    ]

    n_harmonized_maf <- nrow(d)

    if (n_harmonized_maf < 100) {
      return(list(
        status="TOO_FEW_HARMONIZED_VARIANTS",
        region_id=region_id,
        n_harmonized_maf=n_harmonized_maf
      ))
    }

    # SuSiE reference匹配仅使用两个GWAS中一致的rsID
    d <- d[
      !is.na(SNP_fat) &
      SNP_fat != "." &
      SNP_fat != "" &
      SNP_fat == SNP_iron
    ]

    # rsID必须唯一
    d <- d[
      !duplicated(SNP_fat) &
      !duplicated(SNP_fat, fromLast=TRUE)
    ]

    n_unique_rsid <- nrow(d)

    if (n_unique_rsid < 100) {
      return(list(
        status="TOO_FEW_UNIQUE_RSID_VARIANTS",
        region_id=region_id,
        n_harmonized_maf=n_harmonized_maf,
        n_unique_rsid=n_unique_rsid
      ))
    }

    # ----------------------------------------------------------
    # 7.2 读取reference BIM并使用独立字段名
    # ----------------------------------------------------------
    bim <- fread(
      paste0(ref_prefix, ".bim"),
      col.names=c(
        "REF_CHR","SNP","REF_CM","REF_BP","REF_A1","REF_A2"
      ),
      showProgress=FALSE
    )

    bim <- bim[
      !duplicated(SNP) &
      !duplicated(SNP, fromLast=TRUE)
    ]

    setkey(bim, SNP)

    z <- merge(
      d,
      bim,
      by.x="SNP_fat",
      by.y="SNP",
      all=FALSE
    )

    required_merge <- c(
      "SNP_fat","CHR","BP",
      "REF_fat","ALT_fat",
      "BETA_fat","VARBETA_fat",
      "BETA_iron_aligned","VARBETA_iron",
      "REF_CHR","REF_BP","REF_A1","REF_A2"
    )

    if (!all(required_merge %in% names(z))) {
      return(list(
        status=paste0(
          "BIM_MERGE_COLUMNS_UNEXPECTED:",
          paste(
            setdiff(required_merge, names(z)),
            collapse=","
          )
        ),
        region_id=region_id
      ))
    }

    # 坐标必须完全一致
    z <- z[
      CHR == REF_CHR &
      BP == REF_BP
    ]

    n_position_match <- nrow(z)

    if (n_position_match < 100) {
      return(list(
        status="LOW_REFERENCE_POSITION_COVERAGE",
        region_id=region_id,
        n_harmonized_maf=n_harmonized_maf,
        n_unique_rsid=n_unique_rsid,
        n_position_match=n_position_match
      ))
    }

    # ----------------------------------------------------------
    # 7.3 等位基因集合必须一致（允许strand complement）
    # ----------------------------------------------------------
    allele_ok <- (
      (z$REF_fat == z$REF_A1 & z$ALT_fat == z$REF_A2) |
      (z$REF_fat == z$REF_A2 & z$ALT_fat == z$REF_A1) |
      (comp(z$REF_fat) == z$REF_A1 & comp(z$ALT_fat) == z$REF_A2) |
      (comp(z$REF_fat) == z$REF_A2 & comp(z$ALT_fat) == z$REF_A1)
    )

    z <- z[allele_ok]

    n_allele_match <- nrow(z)

    if (n_allele_match < 100) {
      return(list(
        status="LOW_REFERENCE_ALLELE_COVERAGE",
        region_id=region_id,
        n_harmonized_maf=n_harmonized_maf,
        n_unique_rsid=n_unique_rsid,
        n_position_match=n_position_match,
        n_allele_match=n_allele_match
      ))
    }

    # ----------------------------------------------------------
    # 7.4 用PLINK生成最终subset
    # ----------------------------------------------------------
    extract_file <- file.path(region_dir, "extract_rsids.txt")

    fwrite(
      z[,.(SNP_fat)],
      extract_file,
      col.names=FALSE
    )

    subpref <- file.path(region_dir, "ref_subset")
    makebed_log <- file.path(region_dir, "plink_makebed.log")

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

    if (!identical(as.integer(rc1), 0L)) {
      return(list(
        status="PLINK_MAKEBED_FAILED",
        region_id=region_id
      ))
    }

    if (!all(file.exists(paste0(subpref,c(".bed",".bim",".fam"))))) {
      return(list(
        status="PLINK_SUBSET_FILES_MISSING",
        region_id=region_id
      ))
    }

    sbim <- fread(
      paste0(subpref, ".bim"),
      col.names=c(
        "CHR","SNP","CM","BP","A1","A2"
      )
    )

    if (nrow(sbim) < 100) {
      return(list(
        status="TOO_FEW_FINAL_REFERENCE_VARIANTS",
        region_id=region_id,
        n_final_reference=nrow(sbim)
      ))
    }

    if (anyDuplicated(sbim$SNP)) {
      return(list(
        status="FINAL_REFERENCE_DUPLICATED_SNP",
        region_id=region_id
      ))
    }

    # ----------------------------------------------------------
    # 7.5 按最终subset BIM顺序重排GWAS
    # ----------------------------------------------------------
    setkey(z, SNP_fat)

    zz <- z[sbim$SNP]

    if (nrow(zz) != nrow(sbim) ||
        anyNA(zz$SNP_fat) ||
        !identical(as.character(zz$SNP_fat), as.character(sbim$SNP))) {

      return(list(
        status="FINAL_SNP_ORDER_MATCH_FAILED",
        region_id=region_id,
        n_final_reference=nrow(sbim),
        n_matched=nrow(zz)
      ))
    }

    # ----------------------------------------------------------
    # 7.6 生成LD矩阵
    # ----------------------------------------------------------
    ldpref <- file.path(region_dir, "ld")
    ld_log <- file.path(region_dir, "plink_ld.log")

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

    if (!identical(as.integer(rc2), 0L) ||
        !file.exists(ldfile)) {
      return(list(
        status="PLINK_LD_FAILED",
        region_id=region_id
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

    if (nrow(R) != nrow(sbim) ||
        ncol(R) != nrow(sbim)) {
      return(list(
        status="LD_DIMENSION_MISMATCH",
        region_id=region_id,
        ld_nrow=nrow(R),
        ld_ncol=ncol(R),
        bim_n=nrow(sbim)
      ))
    }

    if (any(!is.finite(R))) {
      return(list(
        status="LD_NONFINITE",
        region_id=region_id
      ))
    }

    # 数值对称化
    R <- (R + t(R)) / 2
    diag(R) <- 1

    # ---------------- 核心修复 ----------------
    rownames(R) <- sbim$SNP
    colnames(R) <- sbim$SNP
    # ------------------------------------------

    if (!identical(rownames(R), colnames(R))) {
      return(list(
        status="LD_ROWCOL_NAMES_NOT_IDENTICAL",
        region_id=region_id
      ))
    }

    if (!identical(rownames(R), as.character(sbim$SNP))) {
      return(list(
        status="LD_NAMES_NOT_BIM_ORDER",
        region_id=region_id
      ))
    }

    ld_asym <- safe_cor_symmetry(R)

    if (!is.finite(ld_asym) || ld_asym > 1e-8) {
      return(list(
        status="LD_NOT_SYMMETRIC_AFTER_FIX",
        region_id=region_id,
        max_asymmetry=ld_asym
      ))
    }

    # ----------------------------------------------------------
    # 7.7 effect allele对齐到PLINK A1
    # ----------------------------------------------------------
    alt <- zz$ALT_fat

    alt_to_A1 <- (
      alt == sbim$A1 |
      comp(alt) == sbim$A1
    )

    alt_to_A2 <- (
      alt == sbim$A2 |
      comp(alt) == sbim$A2
    )

    if (any(!(alt_to_A1 | alt_to_A2))) {
      badn <- sum(!(alt_to_A1 | alt_to_A2))
      return(list(
        status="ALLELE_TO_LD_ALIGNMENT_FAILED",
        region_id=region_id,
        n_bad=badn
      ))
    }

    # 对于A/T、C/G palindromic，若两种关系都可能成立，
    # 由于reference和GWAS均已按同一rsID和坐标匹配，且这里只有二等位SNV，
    # 用直接等位基因方向优先；不能唯一确定者剔除更稳妥。
    ambiguous_align <- alt_to_A1 & alt_to_A2

    if (any(ambiguous_align)) {

      keep <- !ambiguous_align

      zz <- zz[keep]
      sbim <- sbim[keep]
      R <- R[keep, keep, drop=FALSE]
      alt_to_A1 <- alt_to_A1[keep]
      alt_to_A2 <- alt_to_A2[keep]

      rownames(R) <- sbim$SNP
      colnames(R) <- sbim$SNP
    }

    if (nrow(sbim) < 100) {
      return(list(
        status="TOO_FEW_AFTER_AMBIGUOUS_ALIGNMENT_FILTER",
        region_id=region_id,
        n_final=nrow(sbim)
      ))
    }

    sign_to_A1 <- ifelse(alt_to_A1, 1, -1)

    beta_fat_ld <- zz$BETA_fat * sign_to_A1
    beta_iron_ld <- zz$BETA_iron_aligned * sign_to_A1

    varbeta_fat <- zz$VARBETA_fat
    varbeta_iron <- zz$VARBETA_iron

    if (any(!is.finite(beta_fat_ld)) ||
        any(!is.finite(beta_iron_ld)) ||
        any(!is.finite(varbeta_fat)) ||
        any(!is.finite(varbeta_iron)) ||
        any(varbeta_fat <= 0) ||
        any(varbeta_iron <= 0)) {

      return(list(
        status="NONFINITE_SUMSTATS_AFTER_ALIGNMENT",
        region_id=region_id
      ))
    }

    # ----------------------------------------------------------
    # 7.8 最终一致性QC
    # ----------------------------------------------------------
    final_snps <- as.character(sbim$SNP)

    if (!identical(rownames(R), final_snps) ||
        !identical(colnames(R), final_snps)) {
      return(list(
        status="FINAL_LD_NAME_CHECK_FAILED",
        region_id=region_id
      ))
    }

    # 保存最终LD/SNP顺序QC
    qc_snp <- data.table(
      order=seq_along(final_snps),
      SNP=final_snps,
      CHR=sbim$CHR,
      BP=sbim$BP,
      LD_A1=sbim$A1,
      LD_A2=sbim$A2,
      GWAS_REF=zz$REF_fat,
      GWAS_ALT=zz$ALT_fat,
      beta_fat_A1=beta_fat_ld,
      beta_iron_A1=beta_iron_ld,
      varbeta_fat=varbeta_fat,
      varbeta_iron=varbeta_iron
    )

    fwrite(
      qc_snp,
      file.path(region_dir, "final_SNP_order_and_alignment.tsv.gz"),
      sep="\t",
      na="NA"
    )

    # ----------------------------------------------------------
    # 7.9 构建coloc/runsusie datasets
    # ----------------------------------------------------------
    D1 <- list(
      beta=beta_fat_ld,
      varbeta=varbeta_fat,
      snp=final_snps,
      position=as.numeric(sbim$BP),
      N=as.numeric(N["pancreas_fat"]),
      type="quant",
      sdY=1,
      LD=R
    )

    D2 <- list(
      beta=beta_iron_ld,
      varbeta=varbeta_iron,
      snp=final_snps,
      position=as.numeric(sbim$BP),
      N=as.numeric(N["pancreas_iron"]),
      type="quant",
      sdY=1,
      LD=R
    )

    # 最后再检查一次
    if (!identical(D1$snp, rownames(D1$LD)) ||
        !identical(D2$snp, rownames(D2$LD))) {
      return(list(
        status="DATASET_SNP_LD_NAME_MISMATCH",
        region_id=region_id
      ))
    }

    # ----------------------------------------------------------
    # 7.10 runsusie
    # ----------------------------------------------------------
    cat("运行 pancreas fat SuSiE-RSS ...\n")

    s1 <- tryCatch(
      coloc::runsusie(
        D1,
        L=10,
        maxit=10000,
        coverage=0.95
      ),
      error=function(e) e
    )

    if (inherits(s1, "error")) {
      return(list(
        status=paste0(
          "RUNSUSIE_FAT_ERROR:",
          conditionMessage(s1)
        ),
        region_id=region_id,
        n_final=nrow(sbim)
      ))
    }

    cat("运行 pancreas iron SuSiE-RSS ...\n")

    s2 <- tryCatch(
      coloc::runsusie(
        D2,
        L=10,
        maxit=10000,
        coverage=0.95
      ),
      error=function(e) e
    )

    if (inherits(s2, "error")) {
      return(list(
        status=paste0(
          "RUNSUSIE_IRON_ERROR:",
          conditionMessage(s2)
        ),
        region_id=region_id,
        n_final=nrow(sbim)
      ))
    }

    # ----------------------------------------------------------
    # 7.11 credible sets
    # ----------------------------------------------------------
    cs1 <- tryCatch(
      susieR::susie_get_cs(
        s1,
        coverage=0.95,
        min_abs_corr=0
      ),
      error=function(e) NULL
    )

    cs2 <- tryCatch(
      susieR::susie_get_cs(
        s2,
        coverage=0.95,
        min_abs_corr=0
      ),
      error=function(e) NULL
    )

    n_cs_fat <- if (!is.null(cs1) && !is.null(cs1$cs))
      length(cs1$cs) else 0L

    n_cs_iron <- if (!is.null(cs2) && !is.null(cs2$cs))
      length(cs2$cs) else 0L

    # ----------------------------------------------------------
    # 7.12 coloc.susie
    # ----------------------------------------------------------
    cat("运行 coloc.susie ...\n")

    csu <- tryCatch(
      coloc::coloc.susie(
        s1,
        s2,
        p1=1e-4,
        p2=1e-4,
        p12=1e-5
      ),
      error=function(e) e
    )

    if (inherits(csu, "error")) {
      return(list(
        status=paste0(
          "COLOC_SUSIE_ERROR:",
          conditionMessage(csu)
        ),
        region_id=region_id,
        n_final=nrow(sbim),
        n_cs_fat=n_cs_fat,
        n_cs_iron=n_cs_iron
      ))
    }

    # ----------------------------------------------------------
    # 7.13 保存完整对象
    # ----------------------------------------------------------
    saveRDS(
      list(
        D1=D1,
        D2=D2,
        susie_fat=s1,
        susie_iron=s2,
        cs_fat=cs1,
        cs_iron=cs2,
        coloc_susie=csu
      ),
      file.path(region_dir, "Step26b_susie_objects.rds")
    )

    # SuSiE credible set index / SNPs
    if (n_cs_fat > 0) {
      fat_cs_rows <- rbindlist(
        lapply(seq_along(cs1$cs), function(i) {
          idx <- cs1$cs[[i]]
          data.table(
            trait="pancreas_fat",
            cs=paste0("CS",i),
            SNP=final_snps[idx],
            BP=sbim$BP[idx]
          )
        }),
        fill=TRUE
      )
      fwrite(
        fat_cs_rows,
        file.path(region_dir, "pancreas_fat_credible_sets.tsv"),
        sep="\t"
      )
    }

    if (n_cs_iron > 0) {
      iron_cs_rows <- rbindlist(
        lapply(seq_along(cs2$cs), function(i) {
          idx <- cs2$cs[[i]]
          data.table(
            trait="pancreas_iron",
            cs=paste0("CS",i),
            SNP=final_snps[idx],
            BP=sbim$BP[idx]
          )
        }),
        fill=TRUE
      )
      fwrite(
        iron_cs_rows,
        file.path(region_dir, "pancreas_iron_credible_sets.tsv"),
        sep="\t"
      )
    }

    # ----------------------------------------------------------
    # 7.14 coloc.susie summary
    # ----------------------------------------------------------
    if (is.null(csu$summary) || nrow(csu$summary)==0) {

      fwrite(
        data.table(
          note="coloc.susie completed but returned no signal-pair summary rows"
        ),
        file.path(region_dir, "coloc_susie_summary.tsv"),
        sep="\t"
      )

      max_h4 <- NA_real_
      n_signal_pairs <- 0L

    } else {

      ss <- as.data.table(csu$summary)

      fwrite(
        ss,
        file.path(region_dir, "coloc_susie_summary.tsv"),
        sep="\t",
        na="NA"
      )

      h4cols <- grep(
        "^PP\\.H4",
        names(ss),
        value=TRUE
      )

      max_h4 <- if (length(h4cols)) {
        max(
          as.numeric(unlist(ss[, ..h4cols])),
          na.rm=TRUE
        )
      } else {
        NA_real_
      }

      if (is.infinite(max_h4)) {
        max_h4 <- NA_real_
      }

      n_signal_pairs <- nrow(ss)
    }

    # SNP-level coloc.susie result if available
    if (!is.null(csu$results)) {
      rr <- as.data.table(csu$results)
      fwrite(
        rr,
        file.path(region_dir, "coloc_susie_SNP_results.tsv.gz"),
        sep="\t",
        na="NA"
      )
    }

    # ----------------------------------------------------------
    # 7.15 区域QC
    # ----------------------------------------------------------
    region_qc <- data.table(
      region_id=region_id,
      n_harmonized_MAF_ge_001=n_harmonized_maf,
      n_unique_matching_rsids=n_unique_rsid,
      n_reference_position_match=n_position_match,
      n_reference_allele_match=n_allele_match,
      n_final_LD_SNPs=nrow(sbim),
      LD_max_asymmetry=ld_asym,
      n_cs_pancreas_fat=n_cs_fat,
      n_cs_pancreas_iron=n_cs_iron,
      n_coloc_susie_signal_pairs=n_signal_pairs,
      max_PP_H4_susie=max_h4
    )

    fwrite(
      region_qc,
      file.path(region_dir, "region_QC.tsv"),
      sep="\t",
      na="NA"
    )

    return(as.list(
      c(
        list(
          status="OK",
          region_id=region_id
        ),
        as.list(region_qc[1])
      )
    ))
  }

  # ============================================================
  # 8. 跑两个区域
  # ============================================================
  results <- list()

  for (rid in target_regions$region_id) {

    result <- tryCatch(
      run_region(rid),
      error=function(e) {
        list(
          status=paste0(
            "UNHANDLED_R_ERROR:",
            conditionMessage(e)
          ),
          region_id=rid
        )
      }
    )

    results[[rid]] <- result

    cat(
      "\n",rid,
      " -> ",
      result$status,
      "\n",
      sep=""
    )
  }

  # ============================================================
  # 9. 汇总状态
  # ============================================================
  status_tab <- rbindlist(
    lapply(names(results), function(rid) {

      x <- results[[rid]]

      # list -> one-row data.table
      dt <- as.data.table(
        lapply(x, function(v) {
          if (length(v)==0) NA else v[1]
        })
      )

      if (!"region_id" %in% names(dt)) {
        dt[, region_id := rid]
      }

      dt
    }),
    fill=TRUE
  )

  # role
  status_tab <- merge(
    target_regions,
    status_tab,
    by="region_id",
    all.x=TRUE
  )

  fwrite(
    status_tab,
    file.path(reports, "Step26b_susie_status.tsv"),
    sep="\t",
    na="NA"
  )

  # ============================================================
  # 10. 将第26步ABF主结果一并写入总表
  # ============================================================
  get_one <- function(dt, nm) {
    if (nm %in% names(dt)) as.numeric(dt[[nm]][1]) else NA_real_
  }

  abf_context <- data.table(
    region_id=c(
      "chr19_LOC2340_primary",
      "chr8_positive_control"
    ),
    ABF_PP_H3=c(
      get_one(abf_chr19,"PP.H3"),
      get_one(abf_chr8,"PP.H3")
    ),
    ABF_PP_H4=c(
      get_one(abf_chr19,"PP.H4"),
      get_one(abf_chr8,"PP.H4")
    )
  )

  final_summary <- merge(
    abf_context,
    status_tab,
    by="region_id",
    all.x=TRUE
  )

  fwrite(
    final_summary,
    file.path(reports, "Step26b_ABF_plus_SuSiE_summary.tsv"),
    sep="\t",
    na="NA"
  )

  # ============================================================
  # 11. 完整性判定
  # ============================================================
  all_ok <- all(
    status_tab$status == "OK",
    na.rm=FALSE
  )

  if (all_ok) {

    completion <- "COMPLETED_STEP26B_SUSIE_COLOCALIZATION"

  } else {

    completion <- "COMPLETED_STEP26B_WITH_REVIEW_FLAGS"
  }

  writeLines(
    completion,
    file.path(reports, "completion_status.txt")
  )

  # ============================================================
  # 12. 方法说明
  # ============================================================
  writeLines(
    c(
      "Step26b repairs the SuSiE-RSS stage of Step26 only.",
      "The Step26 ABF analyses are not rerun.",
      "The LD matrix is generated from the same EUR PLINK reference used in Step26.",
      "The final PLINK subset BIM defines the unique SNP order used for both GWAS summary statistics and the LD matrix.",
      "The LD matrix rownames and colnames are explicitly set to the final BIM SNP IDs before runsusie.",
      "GWAS effect alleles are aligned to the PLINK A1 allele before SuSiE-RSS.",
      "Both traits use N=25,617, quantitative type, sdY=1, L=10, and 95% credible sets.",
      "coloc.susie uses p1=1e-4, p2=1e-4, and p12=1e-5.",
      "The chr19 region is an exploratory candidate selected in Step25b and is not an independent replication.",
      "The chr8 region is used as a positive-control/reproducibility region.",
      "SuSiE/coloc posterior support indicates shared local genetic signals under model assumptions and does not establish causal direction."
    ),
    file.path(reports, "methods_and_interpretation.txt")
  )

  # ============================================================
  # 13. 打包
  # ============================================================
  zipfile <- file.path(
    out,
    "26b_coloc_susie_review.zip"
  )

  oldwd <- getwd()
  on.exit(setwd(oldwd), add=TRUE)
  setwd(out)

  files_to_zip <- c(
    list.files("reports", full.names=TRUE),
    list.files(
      "susie",
      recursive=TRUE,
      full.names=TRUE
    ),
    list.files("logs", full.names=TRUE)
  )

  files_to_zip <- files_to_zip[
    file.exists(files_to_zip) &
    !dir.exists(files_to_zip)
  ]

  utils::zip(
    zipfile,
    files=files_to_zip,
    flags="-q"
  )

  # ============================================================
  # 14. 控制台输出
  # ============================================================
  cat("\n========================================\n")
  cat("第26b步完成\n")
  cat("========================================\n\n")

  print(final_summary)

  cat("\n完成状态：", completion, "\n")
  cat("结果目录：", out, "\n")
  cat("复核包：", zipfile, "\n")

  if (!all_ok) {
    cat(
      "\n注意：至少一个区域仍存在review flag。",
      "请不要继续第27步，先上传26b复核包。\n",
      sep=""
    )
  } else {
    cat(
      "\nchr19和chr8两个区域的SuSiE-RSS + coloc.susie均已成功完成。\n"
    )
  }

})()
