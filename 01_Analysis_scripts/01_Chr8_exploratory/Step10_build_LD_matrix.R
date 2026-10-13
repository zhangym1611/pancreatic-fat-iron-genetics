(function() {

  # ============================================================
  # 10：构建有符号LD矩阵及对齐的GWAS输入
  # 原文件不变；输出到新目录
  # ============================================================

  stopifnot(requireNamespace("data.table", quietly = TRUE))

  base_dir <- path.expand("~/胰腺")
  plink_bin <- path.expand("~/mr_tools/plink")

  audit_file <- file.path(
    base_dir,
    "09_LD_reference_check",
    "chr8_check_c35bc395c198a",
    "LD_reference_audit.rds"
  )

  stopifnot(file.exists(audit_file), file.exists(plink_bin))

  if (file.access(plink_bin, mode = 1) != 0L) {
    stop("PLINK没有执行权限。")
  }

  obj <- readRDS(audit_file)
  ref_prefix <- obj$region_prefix

  if (!all(file.exists(
    paste0(ref_prefix, c(".bed", ".bim", ".fam"))
  ))) {
    stop("上一步的区域参考文件不完整。")
  }

  parent_dir <- file.path(base_dir, "10_LD_matrix")
  dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)

  out_dir <- tempfile("chr8_LD_", tmpdir = parent_dir)
  if (!dir.create(out_dir)) stop("无法创建结果目录。")

  sink(file.path(out_dir, "run_log.txt"), split = TRUE)
  on.exit(sink(), add = TRUE)

  cat("结果目录：", out_dir, "\n")

  traits <- c(
    "pancreas_fat", "pancreas_iron", "pancreas_volume"
  )

  key_snps <- c(
    "rs4733612", "rs1561927", "rs7824269",
    "rs7819154", "rs12676584"
  )

  audit <- as.data.frame(obj$audit)

  a <- audit[
    !is.na(audit$eligible_for_next_check) &
      audit$eligible_for_next_check,
    , drop = FALSE
  ]

  if (!nrow(a) || anyDuplicated(a$rsid)) {
    stop("候选为空或rsID不唯一。")
  }

  stopifnot(
    all(a$sign_ALT_to_A1 %in% c(-1, 1)),
    all(key_snps %in% a$rsid)
  )

  # 1. 写出提取清单与明确的计数等位基因
  extract_file <- file.path(out_dir, "extract_snps.txt")
  allele_file <- file.path(out_dir, "counted_A1_alleles.txt")

  writeLines(a$rsid, extract_file)

  data.table::fwrite(
    a[, c("rsid", "ref_A1"), drop = FALSE],
    allele_file,
    sep = "\t",
    col.names = FALSE
  )

  geno_prefix <- file.path(out_dir, "reference_A1_genotypes")
  console_file <- file.path(out_dir, "plink_console.txt")

  cat("\n正在导出参考基因型……\n")

  status <- system2(
    plink_bin,
    args = c(
      "--bfile", shQuote(ref_prefix),
      "--extract", shQuote(extract_file),
      "--keep-allele-order",
      "--recode", "A",
      "--recode-allele", shQuote(allele_file),
      "--memory", "2048",
      "--threads", "1",
      "--out", shQuote(geno_prefix)
    ),
    stdout = console_file,
    stderr = console_file
  )

  if (length(status) != 1L ||
      is.na(status) || status != 0L) {
    cat(
      tail(readLines(console_file, warn = FALSE), 40L),
      sep = "\n"
    )
    stop("PLINK导出失败。")
  }

  raw_file <- paste0(geno_prefix, ".raw")
  if (!file.exists(raw_file)) stop("未生成.raw文件。")

  raw <- data.table::fread(
    raw_file,
    na.strings = c("NA", "nan"),
    check.names = FALSE,
    showProgress = FALSE
  )

  expected_cols <- paste(a$rsid, a$ref_A1, sep = "_")

  if (anyDuplicated(names(raw)) ||
      !all(expected_cols %in% names(raw))) {
    stop("导出基因型列与预期A1计数不一致。")
  }

  G <- as.matrix(raw[, expected_cols, with = FALSE])
  storage.mode(G) <- "double"
  colnames(G) <- a$rsid

  if (any(!is.na(G) & !(G %in% c(0, 1, 2)))) {
    stop("基因型中存在0/1/2/NA之外的值。")
  }

  if (nrow(G) != obj$reference_n) {
    stop("参考样本数与上一步不一致。")
  }

  cat("参考个体数：", nrow(G), "\n")
  cat("导出候选数：", ncol(G), "\n")

  # 2. 基因型质量检查
  missing_rate <- colMeans(is.na(G))
  ref_AF_A1 <- colMeans(G, na.rm = TRUE) / 2
  genotype_variance <- apply(G, 2L, stats::var, na.rm = TRUE)

  qc <- a[, c(
    "key", "rsid", "pos", "REF", "ALT",
    "ref_A1", "ref_A2", "sign_ALT_to_A1"
  ), drop = FALSE]

  qc$missing_rate <- missing_rate
  qc$ref_AF_A1 <- ref_AF_A1
  qc$genotype_variance <- genotype_variance

  # 本步预设缺失率上限2%；不按GWAS P值筛选
  qc$keep_genotype <- (
    is.finite(qc$missing_rate) &
    qc$missing_rate <= 0.02 &
    is.finite(qc$ref_AF_A1) &
    qc$ref_AF_A1 > 0 & qc$ref_AF_A1 < 1 &
    is.finite(qc$genotype_variance) &
    qc$genotype_variance > 0
  )

  qc$keep_genotype[is.na(qc$keep_genotype)] <- FALSE

  # 3. 频率按相同A1方向比较
  for (trait in traits) {

    d <- as.data.frame(obj$datasets[[trait]])
    idx <- match(qc$key, d$key)

    if (anyNA(idx)) stop("GWAS变异匹配失败：", trait)

    af_A1 <- ifelse(
      qc$sign_ALT_to_A1 == 1,
      d$AF[idx],
      1 - d$AF[idx]
    )

    qc[[paste0(trait, "_AF_A1")]] <- af_A1

    qc[[paste0(trait, "_AF_difference")]] <-
      abs(af_A1 - qc$ref_AF_A1)
  }

  delta_cols <- paste0(traits, "_AF_difference")

  qc$max_AF_difference <- apply(
    as.matrix(qc[, delta_cols, drop = FALSE]),
    1L,
    max
  )

  # 0.10仅作核查标记，不自动删除，不宣称是通用阈值
  qc$frequency_review_flag <- (
    !is.finite(qc$max_AF_difference) |
    qc$max_AF_difference > 0.10
  )

  data.table::fwrite(
    qc,
    file.path(out_dir, "LD_variant_QC.tsv"),
    sep = "\t"
  )

  key_qc <- qc[
    match(key_snps, qc$rsid),
    , drop = FALSE
  ]

  data.table::fwrite(
    key_qc,
    file.path(out_dir, "five_candidate_QC.tsv"),
    sep = "\t"
  )

  cat("\n===== 基因型质量汇总 =====\n")
  cat("输入候选数：", nrow(qc), "\n")
  cat("缺失率超过2%的记录：", sum(qc$missing_rate > 0.02), "\n")
  cat(
    "无有效方差或单态记录：",
    sum(!is.finite(genotype_variance) |
          genotype_variance <= 0),
    "\n"
  )
  cat("基因型检查保留数：", sum(qc$keep_genotype), "\n")
  cat(
    "保留记录中频率差需复核数：",
    sum(qc$keep_genotype & qc$frequency_review_flag),
    "\n"
  )

  if (!all(key_qc$keep_genotype)) {
    print(key_qc, row.names = FALSE)
    stop("关键候选未通过基因型检查，请先返回QC表。")
  }

  if (sum(qc$keep_genotype) < 2L) {
    stop("有效变异不足以构建LD。")
  }

  # 4. 保留有效变异，按位点均值填补少量缺失
  keep <- qc$keep_genotype
  G <- G[, keep, drop = FALSE]
  kept <- qc[keep, , drop = FALSE]

  for (j in seq_len(ncol(G))) {
    miss <- is.na(G[, j])
    if (any(miss)) {
      G[miss, j] <- mean(G[, j], na.rm = TRUE)
    }
  }

  # 标准化后计算相关矩阵，避免逐对删除缺失造成非一致矩阵
  Z <- scale(G, center = TRUE, scale = TRUE)

  if (any(!is.finite(Z))) {
    stop("标准化基因型出现非有限值。")
  }

  LD <- crossprod(Z) / (nrow(Z) - 1)
  dimnames(LD) <- list(kept$rsid, kept$rsid)

  symmetry_error <- max(abs(LD - t(LD)))
  diagonal_error <- max(abs(diag(LD) - 1))

  if (any(!is.finite(LD)) ||
      symmetry_error > 1e-8 ||
      diagonal_error > 1e-8 ||
      max(abs(LD)) > 1 + 1e-8) {
    stop("LD矩阵基本数值检查失败。")
  }

  # 5. 三种表型统一到参考A1方向和同一变异顺序
  aligned <- list()

  for (trait in traits) {

    d <- as.data.frame(obj$datasets[[trait]])
    idx <- match(kept$key, d$key)

    if (anyNA(idx)) stop("对齐输入缺少变异：", trait)

    beta_A1 <- d$ES[idx] * kept$sign_ALT_to_A1

    af_A1 <- ifelse(
      kept$sign_ALT_to_A1 == 1,
      d$AF[idx],
      1 - d$AF[idx]
    )

    if (any(!is.finite(beta_A1)) ||
        any(!is.finite(d$SE[idx])) ||
        any(d$SE[idx] <= 0)) {
      stop("对齐后的效应值异常：", trait)
    }

    aligned[[trait]] <- list(
      snp = kept$rsid,
      position = kept$pos,
      beta = beta_A1,
      varbeta = d$SE[idx]^2,
      MAF = pmin(af_A1, 1 - af_A1),
      N = obj$metadata[[trait]]$N,
      type = "quant"
    )

    stopifnot(
      identical(aligned[[trait]]$snp, colnames(LD))
    )
  }

  # 6. 五个候选之间的LD
  key_r <- LD[key_snps, key_snps, drop = FALSE]
  key_r2 <- key_r^2

  save_matrix_table <- function(m, filename) {
    tab <- data.frame(
      SNP = rownames(m),
      m,
      check.names = FALSE,
      row.names = NULL
    )
    data.table::fwrite(
      tab, file.path(out_dir, filename), sep = "\t"
    )
  }

  save_matrix_table(key_r, "five_candidate_LD_r.tsv")
  save_matrix_table(key_r2, "five_candidate_LD_r2.tsv")

  # 不重复保存三份大型LD矩阵
  saveRDS(
    list(
      LD = LD,
      aligned = aligned,
      variant_QC = qc,
      retained_variants = kept,
      reference_n = nrow(G),
      input_audit_file = audit_file,
      counted_allele = "reference A1",
      missing_method = "per-variant mean imputation",
      missing_rate_limit = 0.02,
      frequency_flag_threshold = 0.10,
      ready_for_susie = FALSE
    ),
    file.path(out_dir, "LD_and_aligned_GWAS.rds")
  )

  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, "sessionInfo.txt")
  )

  cat("\n===== 频率差异汇总：保留变异 =====\n")
  print(summary(kept$max_AF_difference))

  cat("\n===== 五个关键候选QC =====\n")
  print(
    key_qc[, c(
      "rsid", "missing_rate", "ref_AF_A1",
      "max_AF_difference", "frequency_review_flag",
      "keep_genotype"
    )],
    row.names = FALSE
  )

  cat("\n===== 五个关键候选：有符号r =====\n")
  print(round(key_r, 4))

  cat("\n===== 五个关键候选：r平方 =====\n")
  print(round(key_r2, 4))

  cat("\n===== LD矩阵检查 =====\n")
  cat("矩阵维度：", paste(dim(LD), collapse = " × "), "\n")
  cat("最大对称误差：", symmetry_error, "\n")
  cat("最大对角线误差：", diagonal_error, "\n")

  cat("\nLD已保存，但尚未完成GWAS-LD一致性诊断。\n")
  cat("频率差标记未自动用于删除变异。\n")
  cat("尚未运行SuSiE，也未判断独立信号数量。\n")
  cat("结果目录：", out_dir, "\n")

})()