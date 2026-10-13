(function() {

  # ============================================================
  # 09：chr8区域LD参考数据覆盖与等位基因匹配检查
  # 本步不运行SuSiE，不进行LD剪枝，不按P值筛选
  # ============================================================

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("缺少data.table包。")
  }

  # 1. 已确认的路径
  data_dir <- path.expand("~/胰腺")
  plink_bin <- path.expand("~/mr_tools/plink")

  ref_prefix <- file.path(
    data_dir, "00_ld_reference", "EUR"
  )

  input_file <- file.path(
    data_dir,
    "08_pancreas_volume_coloc",
    "chr8_volume_c35bc3fd85136",
    "three_trait_region_inputs.rds"
  )

  required_files <- c(
    input_file,
    plink_bin,
    paste0(ref_prefix, c(".bed", ".bim", ".fam"))
  )

  if (!all(file.exists(required_files))) {
    stop(
      "以下文件缺失：\n",
      paste(
        required_files[!file.exists(required_files)],
        collapse = "\n"
      )
    )
  }

  if (file.access(plink_bin, mode = 1) != 0L) {
    stop("PLINK没有执行权限：", plink_bin)
  }

  parent_dir <- file.path(
    data_dir, "09_LD_reference_check"
  )

  dir.create(
    parent_dir, recursive = TRUE, showWarnings = FALSE
  )

  out_dir <- tempfile(
    pattern = "chr8_check_",
    tmpdir = parent_dir
  )

  if (!dir.create(out_dir)) {
    stop("无法创建输出目录。")
  }

  sink(file.path(out_dir, "run_log.txt"), split = TRUE)
  on.exit(sink(), add = TRUE)

  cat("结果目录：", out_dir, "\n")

  # 2. 读取三表型区域数据
  obj <- readRDS(input_file)

  traits <- c(
    "pancreas_fat",
    "pancreas_iron",
    "pancreas_volume"
  )

  if (!all(traits %in% names(obj$datasets))) {
    stop("RDS中缺少预期表型。")
  }

  if (!identical(obj$build, "GRCh37") ||
      as.character(obj$chromosome) != "8") {
    stop("输入数据的基因组版本或染色体不符。")
  }

  region_start <- as.integer(obj$start)
  region_end <- as.integer(obj$end)

  stopifnot(
    region_start == 129047378L,
    region_end == 130047378L
  )

  datasets <- lapply(
    obj$datasets[traits],
    as.data.frame
  )

  required_cols <- c(
    "key", "rsid", "chr", "pos",
    "REF", "ALT", "ES", "SE", "AF", "LP"
  )

  for (trait in traits) {
    d <- datasets[[trait]]

    if (!all(required_cols %in% names(d))) {
      stop("缺少必要列：", trait)
    }

    if (anyNA(d$key) || anyDuplicated(d$key)) {
      stop("变异键缺失或重复：", trait)
    }
  }

  common_keys <- Reduce(
    intersect,
    lapply(datasets, function(d) d$key)
  )

  if (!length(common_keys)) {
    stop("三种表型没有共同变异。")
  }

  # 按脂肪数据的顺序组织三者共同变异
  anchor <- datasets$pancreas_fat[
    match(common_keys, datasets$pancreas_fat$key),
    , drop = FALSE
  ]

  anchor <- anchor[
    order(anchor$pos, anchor$key),
    , drop = FALSE
  ]

  common_keys <- anchor$key

  cat(
    "\n三表型共同变异数：",
    length(common_keys), "\n"
  )

  # 3. 调用本地PLINK，仅提取目标区域
  # --keep-allele-order保留参考数据中的A1/A2顺序
  region_prefix <- file.path(out_dir, "EUR_chr8_region")

  run_plink <- function(args, log_name) {

    console_file <- file.path(out_dir, log_name)

    status <- system2(
      command = plink_bin,
      args = args,
      stdout = console_file,
      stderr = console_file
    )

    if (length(status) != 1L ||
        is.na(status) || status != 0L) {

      if (file.exists(console_file)) {
        cat(
          tail(readLines(console_file, warn = FALSE), 40L),
          sep = "\n"
        )
      }

      stop("PLINK执行失败，请保留并返回日志。")
    }
  }

  cat("\n正在提取欧洲参考数据中的chr8目标区域……\n")

  run_plink(
    args = c(
      "--bfile", shQuote(ref_prefix),
      "--chr", "8",
      "--from-bp", as.character(region_start),
      "--to-bp", as.character(region_end),
      "--keep-allele-order",
      "--make-bed",
      "--memory", "2048",
      "--threads", "1",
      "--out", shQuote(region_prefix)
    ),
    log_name = "plink_extract_console.txt"
  )

  regional_files <- paste0(
    region_prefix, c(".bed", ".bim", ".fam")
  )

  if (!all(file.exists(regional_files))) {
    stop("PLINK未生成完整区域参考文件。")
  }

  # 4. 读取区域BIM及参考样本数
  bim <- data.table::fread(
    paste0(region_prefix, ".bim"),
    header = FALSE,
    colClasses = "character",
    showProgress = FALSE
  )

  if (ncol(bim) != 6L) {
    stop("区域BIM不是预期的6列。")
  }

  data.table::setnames(
    bim, c("chr", "rsid", "cm", "pos", "A1", "A2")
  )

  bim <- as.data.frame(bim)
  bim$pos <- suppressWarnings(as.integer(bim$pos))

  if (!nrow(bim) || anyNA(bim$pos)) {
    stop("区域BIM为空或位置异常。")
  }

  n_ref <- length(readLines(
    paste0(region_prefix, ".fam"),
    warn = FALSE
  ))

  cat("区域参考变异数：", nrow(bim), "\n")
  cat("参考个体数：", n_ref, "\n")

  # 5. 构建逐变异核查表
  # 本阶段采取保守策略：
  # 仅将rsID唯一、位置一致、非回文SNV且等位基因直接匹配
  # 的变异列为下一步候选；其他记录保留原因，不强制翻转链。
  audit <- anchor[, c(
    "key", "rsid", "chr", "pos", "REF", "ALT"
  ), drop = FALSE]

  audit$is_snv <- (
    audit$REF %in% c("A", "C", "G", "T") &
    audit$ALT %in% c("A", "C", "G", "T") &
    audit$REF != audit$ALT
  )

  audit$palindromic <- (
    paste0(audit$REF, audit$ALT) %in%
      c("AT", "TA", "CG", "GC")
  )

  missing_id <- (
    is.na(audit$rsid) |
    audit$rsid == "" |
    audit$rsid == "."
  )

  duplicate_gwas_id <- (
    duplicated(audit$rsid) |
    duplicated(audit$rsid, fromLast = TRUE)
  )

  duplicate_ref_ids <- unique(
    bim$rsid[
      duplicated(bim$rsid) |
      duplicated(bim$rsid, fromLast = TRUE)
    ]
  )

  j <- match(audit$rsid, bim$rsid)

  audit$ref_pos <- bim$pos[j]
  audit$ref_A1 <- bim$A1[j]
  audit$ref_A2 <- bim$A2[j]

  audit$status <- "pending"

  assign_status <- function(mask, value) {
    mask[is.na(mask)] <- FALSE
    ix <- which(audit$status == "pending" & mask)
    audit$status[ix] <<- value
  }

  assign_status(missing_id, "missing_rsid")
  assign_status(duplicate_gwas_id, "duplicate_gwas_rsid")
  assign_status(is.na(j), "not_in_reference")
  assign_status(
    audit$rsid %in% duplicate_ref_ids,
    "duplicate_reference_rsid"
  )
  assign_status(
    audit$pos != audit$ref_pos,
    "position_mismatch"
  )
  assign_status(!audit$is_snv, "indel_or_non_snv")
  assign_status(audit$palindromic, "palindromic_pending")

  direct <- (
    audit$ALT == audit$ref_A1 &
    audit$REF == audit$ref_A2
  )

  swapped <- (
    audit$ALT == audit$ref_A2 &
    audit$REF == audit$ref_A1
  )

  assign_status(direct, "matched_ALT_is_A1")
  assign_status(swapped, "matched_ALT_is_A2")

  audit$status[audit$status == "pending"] <-
    "allele_mismatch_or_strand_check"

  # 方向因子：未来若LD以A1计数，可将ALT效应转到A1
  # 本步只记录，不修改任何原始beta。
  audit$sign_ALT_to_A1 <- NA_integer_

  audit$sign_ALT_to_A1[
    audit$status == "matched_ALT_is_A1"
  ] <- 1L

  audit$sign_ALT_to_A1[
    audit$status == "matched_ALT_is_A2"
  ] <- -1L

  # 6. 核对三种表型对同一key是否使用一致rsID
  for (trait in traits) {

    d <- datasets[[trait]]
    k <- match(audit$key, d$key)

    stopifnot(
      identical(as.character(d$REF[k]), as.character(audit$REF)),
      identical(as.character(d$ALT[k]), as.character(audit$ALT))
    )

    audit[[paste0(trait, "_beta_ALT")]] <- d$ES[k]
    audit[[paste0(trait, "_se")]] <- d$SE[k]
    audit[[paste0(trait, "_AF_ALT")]] <- d$AF[k]

    same_id <- (
      !is.na(d$rsid[k]) &
      !is.na(audit$rsid) &
      d$rsid[k] == audit$rsid
    )

    same_id[is.na(same_id)] <- FALSE
    audit[[paste0(trait, "_rsid_consistent")]] <- same_id
  }

  id_columns <- paste0(traits, "_rsid_consistent")

  all_ids_consistent <- Reduce(
    `&`, audit[, id_columns, drop = FALSE]
  )

  audit$eligible_for_next_check <- (
    audit$status %in%
      c("matched_ALT_is_A1", "matched_ALT_is_A2") &
    all_ids_consistent
  )

  # 7. 保存覆盖结果和五个关键候选的核查结果
  key_snps <- c(
    "rs4733612", "rs1561927", "rs7824269",
    "rs7819154", "rs12676584"
  )

  key_check <- data.frame(
    rsid = key_snps,
    stringsAsFactors = FALSE
  )

  key_check <- merge(
    key_check,
    audit[, c(
      "rsid", "pos", "REF", "ALT",
      "ref_pos", "ref_A1", "ref_A2",
      "status", "sign_ALT_to_A1",
      "eligible_for_next_check"
    ), drop = FALSE],
    by = "rsid",
    all.x = TRUE,
    sort = FALSE
  )

  key_check <- key_check[
    match(key_snps, key_check$rsid),
    , drop = FALSE
  ]

  key_check$status[is.na(key_check$status)] <-
    "not_in_common_input"

  summary_table <- as.data.frame(
    table(audit$status),
    stringsAsFactors = FALSE
  )

  names(summary_table) <- c("status", "n_variants")

  data.table::fwrite(
    audit,
    file.path(out_dir, "variant_reference_audit.tsv"),
    sep = "\t"
  )

  data.table::fwrite(
    summary_table,
    file.path(out_dir, "reference_coverage_summary.tsv"),
    sep = "\t"
  )

  data.table::fwrite(
    key_check,
    file.path(out_dir, "five_candidate_reference_check.tsv"),
    sep = "\t"
  )

  writeLines(
    audit$rsid[audit$eligible_for_next_check],
    file.path(out_dir, "matched_nonpalindromic_snps.txt")
  )

  saveRDS(
    list(
      input_file = input_file,
      region_prefix = region_prefix,
      reference_prefix = ref_prefix,
      reference_n = n_ref,
      region_start = region_start,
      region_end = region_end,
      audit = audit,
      key_check = key_check,
      datasets = datasets,
      metadata = obj$metadata
    ),
    file.path(out_dir, "LD_reference_audit.rds")
  )

  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, "sessionInfo.txt")
  )

  # 8. 输出汇总
  cat("\n===== 参考数据覆盖情况 =====\n")
  print(summary_table, row.names = FALSE)

  cat(
    "\n可进入下一轮LD质量检查的非回文SNV：",
    sum(audit$eligible_for_next_check),
    "/", nrow(audit), "\n"
  )

  cat("\n===== 五个关键候选检查 =====\n")
  print(key_check, row.names = FALSE)

  cat("\n===== 软件可用性 =====\n")
  cat(
    "susieR可用：",
    requireNamespace("susieR", quietly = TRUE), "\n"
  )
  cat(
    "coloc可用：",
    requireNamespace("coloc", quietly = TRUE), "\n"
  )

  cat("\n本步尚未构建LD矩阵，尚未判断独立信号数量。\n")
  cat("回文位点和indel暂未列入下一步候选，已保留核查记录。\n")
  cat("参考数据中的单态、缺失率和频率差异将在下一步检查。\n")
  cat("结果目录：", out_dir, "\n")

})()