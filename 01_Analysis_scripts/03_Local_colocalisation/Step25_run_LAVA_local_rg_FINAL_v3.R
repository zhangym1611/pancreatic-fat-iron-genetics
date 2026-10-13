# 第25步：LAVA全基因组局部遗传相关（pancreas fat vs pancreas iron）
# 第三次修正版：将LAVA负局部方差视为预期的不可分析结果；捕获其print输出；移除不合理的90%可分析率门槛。
# 主目的：
# 1) 将第22-23步的全基因组rg分解到预定义GRCh37局部遗传区块；
# 2) 找出chr8之外可能贡献于pancreas fat–pancreas iron全局共享遗传背景的区域；
# 3) 为下一步cross-trait colocalization筛选优先区域。
#
# 方法学固定：
# - LAVA v0.1.5；
# - 官方推荐UK Biobank欧洲人群LD reference v1.1（N=100,000），而不是503人的1000G参考；
# - 官方2495个GRCh37/hg19 LAVA partition blocks；
# - 输入沿用第20步严格QC后的HapMap3共同框架（SNP/A1/A2/Z/N）；
# - sample-overlap相关矩阵由第21步单表型LDSC intercept和第22步cross-trait gencov intercept构建；
# - 先做每个locus的local h2；只有两个表型均通过Bonferroni local-h2阈值(0.05/2495)
#   才进行bivariate local rg，避免在局部遗传信号不足时解释不稳定的相关系数；
# - 对实际完成的bivariate tests同时给出BH-FDR和Bonferroni校正；
# - 不把local rg解释为因果关系，也不把物理区域信号直接解释为共定位。
#
# 第一次运行将下载并解压约15 GB的UKB LAVA LD reference；下载为断点续传，解压后zip会删除。
# 后续重复运行会直接复用已经完整的reference。

(function() {
  base <- path.expand("~/胰腺")
  step20 <- file.path(base, "20_LDSC_inputs", "GWAS_inputs_16c90c74e8c0e1")
  step21 <- file.path(base, "21_LDSC_heritability", "h2_16c90c661ce88d")

  if (!requireNamespace("data.table", quietly = TRUE))
    stop("缺少R包 data.table。")

  check_status <- function(path, expected) {
    if (!file.exists(path)) stop("缺少完成标记：", path)
    x <- readLines(path, warn = FALSE)
    if (!length(x) || !identical(x[1], expected))
      stop("完成标记不匹配：", path, "\n实际：", paste(x, collapse = " | "))
  }

  check_status(
    file.path(step20, "reports", "completion_status.txt"),
    "COMPLETED_INPUT_PREPARATION_REVIEW_QC_BEFORE_LDSC"
  )
  check_status(
    file.path(step21, "reports", "completion_status.txt"),
    "COMPLETED_8_TRAIT_H2_REVIEW_ESTIMATES"
  )

  # -------- 1. 自动定位唯一成功的第22步目录 --------
  step22_parent <- file.path(base, "22_LDSC_genetic_correlation")
  if (!dir.exists(step22_parent)) stop("未找到第22步目录：", step22_parent)
  dirs22 <- list.dirs(step22_parent, recursive = FALSE, full.names = TRUE)
  ok22 <- vapply(dirs22, function(d) {
    f <- file.path(d, "reports", "completion_status.txt")
    file.exists(f) &&
      identical(readLines(f, warn = FALSE)[1], "COMPLETED_28_PAIR_FITS_REVIEW_ESTIMATES")
  }, logical(1))
  good22 <- dirs22[ok22]
  if (length(good22) != 1L) {
    cat("检测到的成功第22步目录：\n", paste(good22, collapse = "\n"), "\n")
    stop("要求恰好1个已通过核查的第22步目录；请保留目标批次或在脚本顶部手动固定路径。")
  }
  step22 <- good22[1]

  # -------- 2. 输出及持久化软件/reference目录 --------
  parent <- file.path(base, "25_LAVA_local_rg")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile("fat_iron_local_rg_", tmpdir = parent)
  if (!dir.create(out)) stop("无法建立第25步输出目录：", out)

  reports <- file.path(out, "reports")
  logs <- file.path(out, "logs")
  checkpoints <- file.path(out, "checkpoints")
  for (d in c(reports, logs, checkpoints)) dir.create(d)

  writeLines("INCOMPLETE_CHECK_LOGS", file.path(reports, "completion_status.txt"))
  writeLines(capture.output(sessionInfo()), file.path(reports, "R_sessionInfo_before.txt"))

  script_source <- file.path(base, "25_run_LAVA_local_rg.R")
  if (file.exists(script_source))
    file.copy(script_source, file.path(reports, basename(script_source)), overwrite = FALSE)

  # 独立用户库，避免改变此前R环境。
  lava_lib <- path.expand("~/R_pancreas_LAVA_015")
  dir.create(lava_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lava_lib, .libPaths()))

  # 持久reference，不放在临时结果目录。
  lava_root <- file.path(base, "lava_reference")
  refroot <- file.path(lava_root, "UKB_v1.1")
  dir.create(refroot, recursive = TRUE, showWarnings = FALSE)

  # -------- 3. 安装并固定LAVA v0.1.5 --------
  install_log <- file.path(logs, "LAVA_install.txt")
  if (!requireNamespace("LAVA", quietly = TRUE) ||
      as.character(utils::packageVersion("LAVA")) != "0.1.5") {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", lib = lava_lib, repos = "https://cloud.r-project.org")
    }
    zz <- file(install_log, open = "wt")
    sink(zz, type = "output")
    sink(zz, type = "message")
    install_ok <- TRUE
    tryCatch(
      remotes::install_github(
        "josefin-werme/LAVA@v0.1.5",
        lib = lava_lib,
        dependencies = TRUE,
        upgrade = "never",
        build_vignettes = FALSE,
        quiet = FALSE
      ),
      error = function(e) {
        install_ok <<- FALSE
        message(conditionMessage(e))
      }
    )
    sink(type = "message")
    sink(type = "output")
    close(zz)
    if (!install_ok) stop("LAVA v0.1.5安装失败，请查看：", install_log)
  }
  if (!requireNamespace("LAVA", quietly = TRUE))
    stop("LAVA安装后仍不可用。")
  if (as.character(utils::packageVersion("LAVA")) != "0.1.5")
    stop("要求LAVA==0.1.5；当前版本：", as.character(utils::packageVersion("LAVA")))

  suppressPackageStartupMessages(library(LAVA))
  writeLines(
    paste0("LAVA_version\t", as.character(utils::packageVersion("LAVA"))),
    file.path(reports, "software_version.tsv")
  )

  # -------- 4. 下载官方UKB v1.1 LAVA binary LD reference --------
  # 官方REFERENCE.md列出的7个归档链接。
  archives <- data.frame(
    label = c("chr1-2","chr3-4","chr5-6","chr7-9","chr10-12","chr13-16","chr17-23"),
    url = c(
      "https://vu.data.surfsara.nl/index.php/s/7NBVIvtPRdu7Qhz/download",
      "https://vu.data.surfsara.nl/index.php/s/fy6ITboMojHrQXr/download",
      "https://vu.data.surfsara.nl/index.php/s/mRz31q0lq7KMcuI/download",
      "https://vu.data.surfsara.nl/index.php/s/89RXxPN2BlOqxLs/download",
      "https://vu.data.surfsara.nl/index.php/s/3dU1L2Hap43xuCs/download",
      "https://vu.data.surfsara.nl/index.php/s/DQqJ2Sqr49RP4xe/download",
      "https://vu.data.surfsara.nl/index.php/s/U8eg5XfTr8qPzPp/download"
    ),
    stringsAsFactors = FALSE
  )
  archives$zip <- file.path(refroot, paste0("lava-ukb-v1.1_", archives$label, ".zip"))
  data.table::fwrite(archives, file.path(reports, "UKB_reference_download_manifest.tsv"), sep = "\t")

  expected_ref <- unlist(lapply(1:22, function(chr)
    file.path(refroot, paste0("lava-ukb-v1.1_chr", chr, c(".info", ".bcor")))))
  reference_complete <- all(file.exists(expected_ref) & file.info(expected_ref)$size > 0)

  if (!reference_complete) {
    wget <- Sys.which("wget")
    curl <- Sys.which("curl")
    unzip_bin <- Sys.which("unzip")
    if (!nzchar(unzip_bin))
      stop("系统缺少unzip。请让服务器提供unzip后再运行；不建议用不明第三方解压工具处理15GB reference。")

    # 只给出磁盘空间提示，不自动删除用户数据。
    df <- tryCatch(system2("df", c("-Pk", shQuote(refroot)), stdout = TRUE, stderr = TRUE),
                   error = function(e) character())
    if (length(df) >= 2) {
      writeLines(df, file.path(reports, "disk_space_before_reference_download.txt"))
      cat("参考文件尚不完整。UKB LAVA reference解压后约15GB，请确认服务器剩余空间充足。\n")
    }

    for (i in seq_len(nrow(archives))) {
      # 如果该归档负责的染色体已经全部存在，则跳过下载。
      rng <- switch(
        archives$label[i],
        "chr1-2" = 1:2,
        "chr3-4" = 3:4,
        "chr5-6" = 5:6,
        "chr7-9" = 7:9,
        "chr10-12" = 10:12,
        "chr13-16" = 13:16,
        "chr17-23" = 17:22
      )
      needed <- unlist(lapply(rng, function(chr)
        file.path(refroot, paste0("lava-ukb-v1.1_chr", chr, c(".info", ".bcor")))))
      if (all(file.exists(needed) & file.info(needed)$size > 0)) {
        cat("已存在：", archives$label[i], "\n")
        next
      }

      zf <- archives$zip[i]
      url <- archives$url[i]
      cat("\n下载UKB LAVA reference：", archives$label[i], "\n")
      flush.console()

      if (nzchar(wget)) {
        rc <- system2(
          wget,
          c("-c", "--tries=5", "--timeout=60", "-O", shQuote(zf), shQuote(url)),
          stdout = file.path(logs, paste0("wget_", archives$label[i], ".txt")),
          stderr = file.path(logs, paste0("wget_", archives$label[i], ".txt"))
        )
      } else if (nzchar(curl)) {
        rc <- system2(
          curl,
          c("-L", "--retry", "5", "--connect-timeout", "60", "-C", "-", "-o", shQuote(zf), shQuote(url)),
          stdout = file.path(logs, paste0("curl_", archives$label[i], ".txt")),
          stderr = file.path(logs, paste0("curl_", archives$label[i], ".txt"))
        )
      } else {
        stop("系统没有wget或curl，无法可靠断点续传大型reference。")
      }
      if (!identical(as.integer(rc), 0L) || !file.exists(zf) || file.info(zf)$size <= 0)
        stop("reference下载失败：", archives$label[i], "；查看logs。")

      # -j：扁平化目录，使LAVA前缀固定为refroot/lava-ukb-v1.1
      ulog <- file.path(logs, paste0("unzip_", archives$label[i], ".txt"))
      urc <- system2(
        unzip_bin,
        c("-j", "-o", shQuote(zf), "-d", shQuote(refroot)),
        stdout = ulog, stderr = ulog
      )
      if (!identical(as.integer(urc), 0L))
        stop("reference解压失败：", archives$label[i], "；查看：", ulog)

      # 确认负责的常染色体都已产生非空.info/.bcor后才删除zip。
      if (!all(file.exists(needed) & file.info(needed)$size > 0))
        stop("解压后预期reference文件不完整：", archives$label[i], "。请保留zip并核查日志。")
      unlink(zf)
    }
  }

  if (!all(file.exists(expected_ref) & file.info(expected_ref)$size > 0))
    stop("UKB LAVA v1.1 reference仍不完整。")

  ref_manifest <- data.frame(
    file = expected_ref,
    bytes = file.info(expected_ref)$size,
    mtime = as.character(file.info(expected_ref)$mtime),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(ref_manifest, file.path(reports, "UKB_reference_file_manifest.tsv"), sep = "\t")
  ref_prefix <- file.path(refroot, "lava-ukb-v1.1")

  # -------- 5. 固定官方2495个GRCh37/hg19 LAVA blocks --------
  # 先验证已有文件；若上一次下载中断留下残缺文件，则删除并可靠重下。
  locus_file <- file.path(lava_root, "blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
  locus_tmp <- paste0(locus_file, ".part")
  locus_url <- paste0(
    "https://raw.githubusercontent.com/josefin-werme/LAVA/v0.1.5/",
    "support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
  )

  validate_locus_file <- function(f) {
    if (!file.exists(f) || file.info(f)$size <= 0) return(FALSE)
    ok <- tryCatch({
      x <- LAVA::read.loci(f)
      nrow(x) == 2495L &&
        all(c("LOC","CHR","START","STOP") %in% names(x)) &&
        all(x$CHR %in% 1:22)
    }, error=function(e) FALSE)
    isTRUE(ok)
  }

  if (file.exists(locus_file) && !validate_locus_file(locus_file)) {
    cat("检测到不完整或异常的LAVA locus文件，删除并重新下载。\n")
    unlink(locus_file)
  }
  if (file.exists(locus_tmp)) unlink(locus_tmp)

  if (!validate_locus_file(locus_file)) {
    wget <- Sys.which("wget")
    curl <- Sys.which("curl")
    cat("下载LAVA官方2495-locus文件……\n")
    flush.console()

    success <- FALSE
    if (nzchar(wget)) {
      rc <- system2(
        wget,
        c("--tries=10","--timeout=120","--retry-connrefused",
          "-O", shQuote(locus_tmp), shQuote(locus_url)),
        stdout=file.path(logs,"wget_locus_file.txt"),
        stderr=file.path(logs,"wget_locus_file.txt")
      )
      success <- identical(as.integer(rc),0L) &&
        file.exists(locus_tmp) && file.info(locus_tmp)$size > 60000
    } else if (nzchar(curl)) {
      rc <- system2(
        curl,
        c("-L","--retry","10","--retry-delay","3",
          "--connect-timeout","120","--max-time","600",
          "-o", shQuote(locus_tmp), shQuote(locus_url)),
        stdout=file.path(logs,"curl_locus_file.txt"),
        stderr=file.path(logs,"curl_locus_file.txt")
      )
      success <- identical(as.integer(rc),0L) &&
        file.exists(locus_tmp) && file.info(locus_tmp)$size > 60000
    } else {
      old_timeout <- getOption("timeout")
      options(timeout=600)
      rc <- tryCatch(
        download.file(locus_url, locus_tmp, mode="wb", method="libcurl", quiet=FALSE),
        error=function(e) 1L
      )
      options(timeout=old_timeout)
      success <- identical(as.integer(rc),0L) &&
        file.exists(locus_tmp) && file.info(locus_tmp)$size > 60000
    }

    if (!success)
      stop("LAVA locus文件未完整下载。请检查网络或logs目录。")
    if (!validate_locus_file(locus_tmp))
      stop("下载的LAVA locus文件未通过2495-locus完整性验证。")
    if (!file.rename(locus_tmp, locus_file))
      stop("LAVA locus文件验证通过，但无法移动到正式路径。")
  }

  loci <- LAVA::read.loci(locus_file)
  if (nrow(loci) != 2495L) stop("预期2495个LAVA loci，实际：", nrow(loci))
  if (!all(c("LOC","CHR","START","STOP") %in% names(loci)))
    stop("locus文件字段异常。")
  if (!all(loci$CHR %in% 1:22)) stop("locus文件含非预期染色体。")
  data.table::fwrite(loci, file.path(reports, "LAVA_2495_loci.tsv"), sep = "\t")

  # -------- 6. 构建LAVA输入（复用第20步高质量SNP/A1/A2/Z/N文件） --------
  traits <- c("pancreas_fat", "pancreas_iron")
  ids <- c(
    pancreas_fat  = "ebi-a-GCST90016675",
    pancreas_iron = "ebi-a-GCST90016676"
  )
  sumstats <- file.path(step20, "sumstats", paste0(ids, ".sumstats.gz"))
  names(sumstats) <- traits
  if (!all(file.exists(sumstats))) stop("缺少第20步sumstats：", paste(sumstats[!file.exists(sumstats)], collapse = "\n"))

  # 再次验证第20步两个目标sumstats的表头、行数和基本数值。
  # 不再依赖shell的 zcat | wc -l，避免不同服务器环境下返回非数值结果。
  expected_n <- c(
    pancreas_fat  = 1171383L,
    pancreas_iron = 1171383L
  )
  expected_N <- c(
    pancreas_fat  = 25617L,
    pancreas_iron = 25617L
  )

  sumstat_qc <- lapply(traits, function(t) {
    f <- unname(sumstats[t])
    if (!file.exists(f)) stop("第20步sumstats文件不存在：", f)

    cat("核查LAVA输入：", t, "\n")
    flush.console()

    d <- data.table::fread(f, showProgress = FALSE)
    req <- c("SNP","A1","A2","Z","N")

    if (!identical(names(d), req)) {
      stop(
        "第20步sumstats表头与预期不一致：", t,
        "\n预期：", paste(req, collapse=","),
        "\n实际：", paste(names(d), collapse=",")
      )
    }

    n <- nrow(d)
    if (as.integer(n) != as.integer(expected_n[t])) {
      stop(
        "第20步sumstats SNP数量与既往QC结果不一致：", t,
        "\n预期：", expected_n[t],
        "\n实际：", n
      )
    }

    if (anyNA(d$SNP) || anyDuplicated(d$SNP))
      stop("SNP存在缺失或重复：", t)

    if (anyNA(d$A1) || anyNA(d$A2) ||
        any(!d$A1 %in% c("A","C","G","T")) ||
        any(!d$A2 %in% c("A","C","G","T")))
      stop("等位基因字段异常：", t)

    if (any(!is.finite(d$Z)))
      stop("存在非有限Z值：", t)

    if (any(!is.finite(d$N)) || any(d$N <= 0))
      stop("样本量N异常：", t)

    if (!all(d$N == expected_N[t])) {
      stop(
        "逐SNP N与第20步固定样本量不一致：", t,
        "\n预期N：", expected_N[t]
      )
    }

    ans <- data.frame(
      trait = t,
      id = unname(ids[t]),
      file = f,
      n_rows = n,
      expected_n_rows = unname(expected_n[t]),
      N = unname(expected_N[t]),
      SNP_unique = !anyDuplicated(d$SNP),
      basic_QC = "PASS",
      stringsAsFactors = FALSE
    )

    rm(d)
    gc()
    ans
  })

  sumstat_qc <- do.call(rbind, sumstat_qc)

  if (nrow(sumstat_qc) != 2L || any(sumstat_qc$basic_QC != "PASS"))
    stop("第20步两个目标sumstats复核未通过。")

  data.table::fwrite(
    sumstat_qc,
    file.path(reports, "LAVA_sumstats_audit.tsv"),
    sep = "\t"
  )

  cat("\n第20步两个LAVA输入文件复核通过：\n")
  print(sumstat_qc)

  info <- data.frame(
    phenotype = traits,
    cases = NA_integer_,
    controls = NA_integer_,
    filename = unname(sumstats[traits]),
    stringsAsFactors = FALSE
  )
  info_file <- file.path(out, "LAVA_input_info.tsv")
  data.table::fwrite(info, info_file, sep = "\t", na = "NA")

  # -------- 7. 从第21/22步构建sampling-correlation矩阵 --------
  h2_file <- file.path(step21, "reports", "LDSC_h2_summary.tsv")
  rg_file <- file.path(step22, "reports", "LDSC_rg_raw_summary.tsv")
  if (!all(file.exists(c(h2_file, rg_file))))
    stop("缺少第21/22步LDSC结果。")

  h2 <- as.data.frame(data.table::fread(h2_file))
  rg <- as.data.frame(data.table::fread(rg_file))
  if (!all(c("trait","intercept","run_status") %in% names(h2)))
    stop("第21步h2表缺少intercept字段。")
  if (!all(c("trait1","trait2","gencov_intercept","run_status") %in% names(rg)))
    stop("第22步rg表缺少gencov_intercept字段。")

  hh <- h2[h2$trait %in% traits & h2$run_status == "COMPLETED", , drop=FALSE]
  if (nrow(hh) != 2L) stop("第21步无法唯一提取两个目标表型intercept。")
  int_fat <- hh$intercept[match("pancreas_fat", hh$trait)]
  int_iron <- hh$intercept[match("pancreas_iron", hh$trait)]

  rr <- rg[
    ((rg$trait1=="pancreas_fat" & rg$trait2=="pancreas_iron") |
     (rg$trait1=="pancreas_iron" & rg$trait2=="pancreas_fat")) &
      rg$run_status=="COMPLETED",
    , drop=FALSE
  ]
  if (nrow(rr) != 1L) stop("第22步无法唯一提取pancreas_fat–pancreas_iron gencov intercept。")
  gcov_int <- as.numeric(rr$gencov_intercept[1])

  if (any(!is.finite(c(int_fat,int_iron,gcov_int))) || int_fat <= 0 || int_iron <= 0)
    stop("LDSC intercept无法构建有效sampling-correlation矩阵。")

  int_cov <- matrix(
    c(int_fat, gcov_int,
      gcov_int, int_iron),
    nrow=2, byrow=TRUE,
    dimnames=list(traits,traits)
  )
  eig <- eigen(int_cov, symmetric=TRUE, only.values=TRUE)$values
  if (min(eig) <= 0)
    stop("由LDSC intercept构成的sampling covariance matrix非正定；不应强行修正。")

  overlap_cor <- cov2cor(int_cov)
  if (any(!is.finite(overlap_cor)) || any(abs(overlap_cor) > 1 + 1e-10))
    stop("sampling correlation矩阵异常。")
  overlap_file <- file.path(out, "LAVA_sample_overlap.txt")
  write.table(overlap_cor, overlap_file, quote=FALSE, sep="\t", col.names=NA)

  overlap_audit <- data.frame(
    quantity = c("pancreas_fat_h2_intercept","pancreas_iron_h2_intercept",
                 "fat_iron_gencov_intercept","derived_sampling_correlation",
                 "sampling_cov_min_eigenvalue"),
    value = c(int_fat,int_iron,gcov_int,overlap_cor[1,2],min(eig))
  )
  data.table::fwrite(overlap_audit, file.path(reports, "sample_overlap_audit.tsv"), sep="\t")

  # -------- 8. 记录预先固定的分析规则 --------
  univ_thresh <- 0.05 / nrow(loci)
  writeLines(c(
    "Step25: genome-wide local genetic correlation using LAVA v0.1.5.",
    "Primary pair: pancreas_fat vs pancreas_iron.",
    "Reference: official LAVA UK Biobank European LD v1.1, N=100,000.",
    "Locus partition: official 2495 GRCh37/hg19 LAVA blocks.",
    "GWAS input: Step20 harmonized HapMap3 SNP/A1/A2/Z/N files; this is a conservative common-SNP local analysis.",
    paste0("Univariate local-h2 screening threshold: 0.05/2495 = ", format(univ_thresh, scientific=TRUE, digits=12), "."),
    "LAVA process.locus may drop a phenotype when its local genetic variance estimate is negative; this is treated as an expected non-analysable phenotype-locus result, not as a data error.",
    "If all target phenotypes have negative local variance and process.locus returns NULL, the locus is recorded as EXPECTED_ALL_NEGATIVE_VARIANCE rather than a pipeline failure.",
    "The analysis is considered complete when all 2495 loci are attempted and no unclassified/R-exception processing events remain; there is no arbitrary minimum percentage of analyzable loci.",
    "Univariate analysis is run only for phenotypes retained in locus$phenos after process.locus filtering.",
    "Bivariate local rg is only estimated when BOTH target phenotypes remain available and BOTH pass the prespecified local-h2 threshold.",
    "Sampling correlation is derived from the Step21 h2 LDSC intercepts and Step22 cross-trait genetic-covariance intercept, standardized with cov2cor, following LAVA guidance.",
    "For bivariate tests, raw P, BH-FDR across all performed local-rg tests, Bonferroni across performed tests, and an additional conservative 2495-locus Bonferroni value are all retained.",
    "No locus is selected for testing based on whether its bivariate rg appears favorable.",
    "LAVA local rg identifies localized genetic sharing; it does not establish causality or colocalization.",
    "The previously studied chr8 interval GRCh37 129047378-130047378 is flagged for context only."
  ), file.path(reports, "methods_and_interpretation.txt"))

  # -------- 9. process.input --------
  process_log <- file.path(logs, "process_input_console.txt")
  zz <- file(process_log, open="wt")
  sink(zz, type="output"); sink(zz, type="message")
  input <- tryCatch(
    LAVA::process.input(
      input.info.file = info_file,
      sample.overlap.file = overlap_file,
      ref.prefix = ref_prefix,
      phenos = traits
    ),
    error=function(e) e
  )
  sink(type="message"); sink(type="output"); close(zz)
  if (inherits(input, "error"))
    stop("LAVA process.input失败：", conditionMessage(input), "；查看：", process_log)

  # 输入审计
  input_audit <- data.frame(
    trait = traits,
    n_aligned_sumstats = vapply(traits, function(t) nrow(input$sum.stats[[t]]), integer(1)),
    stringsAsFactors=FALSE
  )
  if ("unalignable.snps" %in% ls(input)) {
    input_audit$n_unalignable <- vapply(traits, function(t) {
      x <- input$unalignable.snps[[t]]
      if (is.null(x)) 0L else length(x)
    }, integer(1))
  } else input_audit$n_unalignable <- NA_integer_
  data.table::fwrite(input_audit, file.path(reports, "LAVA_process_input_audit.tsv"), sep="\t")
  saveRDS(input$info, file.path(reports, "LAVA_processed_input_info.rds"))

  # -------- 10. 顺序分析2495 loci --------
  chr8_start <- 129047378L
  chr8_end <- 130047378L

  univ_rows <- list()
  bivar_rows <- list()
  locus_qc <- list()
  warnings_seen <- list()

  cat("\n开始LAVA局部分析，共", nrow(loci), "个loci。\n")
  cat("local h2 Bonferroni阈值：", signif(univ_thresh, 5), "\n")
  flush.console()

  for (i in seq_len(nrow(loci))) {
    if (i %% 50L == 0L || i == 1L)
      cat("LAVA locus", i, "/", nrow(loci), " (chr", loci$CHR[i], ":", loci$START[i], "-", loci$STOP[i], ")\n", sep="")

    locus <- NULL
    warn_i <- character()
    process_error <- NULL
    process_print <- character()

    # LAVA v0.1.5 对“负方差”使用 print() 输出 "Error:" / "Warning:"，
    # 这并不是R异常。这里将其捕获到日志，避免控制台被数千条预期信息淹没。
    process_print <- capture.output({
      locus <- tryCatch(
        withCallingHandlers(
          LAVA::process.locus(loci[i, , drop=FALSE], input),
          warning=function(w) {
            warn_i <<- c(warn_i, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        ),
        error=function(e) {
          process_error <<- conditionMessage(e)
          NULL
        }
      )
    }, type="output")

    if (length(warn_i)) {
      warnings_seen[[length(warnings_seen)+1L]] <- data.frame(
        LOC=loci$LOC[i],
        warning=paste(unique(warn_i), collapse=" | "),
        stringsAsFactors=FALSE
      )
    }

    # 保存LAVA自身print的信息用于审计，但不再把预期的negative-variance当作程序错误。
    if (length(process_print)) {
      warnings_seen[[length(warnings_seen)+1L]] <- data.frame(
        LOC=loci$LOC[i],
        warning=paste(process_print, collapse=" | "),
        stringsAsFactors=FALSE
      )
    }

    if (!is.null(process_error)) {
      locus_qc[[length(locus_qc)+1L]] <- data.frame(
        LOC=loci$LOC[i], CHR=loci$CHR[i], START=loci$START[i], STOP=loci$STOP[i],
        process_status="R_EXCEPTION",
        n_snps=NA_integer_, n_pcs=NA_integer_,
        n_available_phenos=NA_integer_, available_phenos="",
        both_target_phenos_available=FALSE,
        message=process_error, stringsAsFactors=FALSE
      )
      next
    }

    if (is.null(locus)) {
      msg <- paste(process_print, collapse=" | ")

      expected_all_negative <- grepl(
        "Negative variance estimate for all phenotypes",
        msg, fixed=TRUE
      )
      expected_all_failed <- grepl(
        "Processing of all phenotypes",
        msg, fixed=TRUE
      )
      too_few_snps <- grepl("Fewer than", msg, fixed=TRUE) &&
        grepl("SNPs in locus", msg, fixed=TRUE)
      too_few_pcs <- grepl("Fewer than", msg, fixed=TRUE) &&
        grepl("PCs in locus", msg, fixed=TRUE)

      status <- if (expected_all_negative) {
        "EXPECTED_ALL_NEGATIVE_VARIANCE"
      } else if (expected_all_failed) {
        "EXPECTED_ALL_PHENOTYPES_FAILED"
      } else if (too_few_snps) {
        "EXPECTED_TOO_FEW_SNPS"
      } else if (too_few_pcs) {
        "EXPECTED_TOO_FEW_PCS"
      } else {
        "NULL_UNCLASSIFIED_REVIEW"
      }

      locus_qc[[length(locus_qc)+1L]] <- data.frame(
        LOC=loci$LOC[i], CHR=loci$CHR[i], START=loci$START[i], STOP=loci$STOP[i],
        process_status=status,
        n_snps=NA_integer_, n_pcs=NA_integer_,
        n_available_phenos=0L, available_phenos="",
        both_target_phenos_available=FALSE,
        message=if (nzchar(msg)) msg else "process.locus returned NULL without printed message",
        stringsAsFactors=FALSE
      )
      next
    }

    # process.locus()会在局部遗传方差估计为负时自动丢弃对应表型。
    # 因此后续只能分析仍保留在locus$phenos中的目标表型，不能再强制传入全部traits。
    available_phenos <- intersect(traits, as.character(locus$phenos))
    n_available <- length(available_phenos)

    returned_status <- if (
      n_available == 2L && all(traits %in% available_phenos)
    ) {
      "ANALYSABLE_BOTH_PHENOTYPES"
    } else if (n_available == 1L) {
      "ANALYSABLE_ONE_PHENOTYPE"
    } else {
      "RETURNED_NO_TARGET_PHENOTYPE_REVIEW"
    }

    locus_qc[[length(locus_qc)+1L]] <- data.frame(
      LOC=loci$LOC[i], CHR=loci$CHR[i], START=loci$START[i], STOP=loci$STOP[i],
      process_status=returned_status, n_snps=locus$n.snps, n_pcs=locus$K,
      n_available_phenos=n_available,
      available_phenos=if (n_available) paste(available_phenos, collapse=";") else "",
      both_target_phenos_available=(n_available==2L && all(traits %in% available_phenos)),
      message=if (length(process_print)) paste(process_print, collapse=" | ") else "",
      stringsAsFactors=FALSE
    )

    # 极少数情况下locus对象存在但两个目标表型均已被LAVA丢弃：只记录QC，不继续分析。
    if (n_available == 0L) {
      warnings_seen[[length(warnings_seen)+1L]] <- data.frame(
        LOC=loci$LOC[i],
        warning="No target phenotype remained in locus after process.locus variance filtering",
        stringsAsFactors=FALSE
      )
      next
    }

    # 只对LAVA在当前locus中实际保留的表型运行univariate local h2。
    u <- tryCatch(
      LAVA::run.univ(
        locus,
        phenos=available_phenos,
        var=TRUE,
        cap.estimates=FALSE
      ),
      error=function(e) e
    )

    if (inherits(u, "error")) {
      warnings_seen[[length(warnings_seen)+1L]] <- data.frame(
        LOC=loci$LOC[i], warning=paste0("run.univ ERROR: ", conditionMessage(u)),
        stringsAsFactors=FALSE
      )
      next
    }

    if (is.null(u) ||
        nrow(u) != n_available ||
        !all(available_phenos %in% u$phen)) {
      warnings_seen[[length(warnings_seen)+1L]] <- data.frame(
        LOC=loci$LOC[i],
        warning=paste0(
          "run.univ returned unexpected structure; available phenotypes=",
          paste(available_phenos, collapse=";")
        ),
        stringsAsFactors=FALSE
      )
      next
    }

    u$LOC <- loci$LOC[i]
    u$CHR <- loci$CHR[i]
    u$START <- loci$START[i]
    u$STOP <- loci$STOP[i]
    u$n_snps <- locus$n.snps
    u$n_pcs <- locus$K
    u$n_available_phenos <- n_available
    u$both_target_phenos_available <- (
      n_available == 2L && all(traits %in% available_phenos)
    )
    u$passes_univ_bonf <- is.finite(u$p) & u$p < univ_thresh
    u$overlaps_known_chr8 <- loci$CHR[i] == 8L &&
      loci$START[i] <= chr8_end && loci$STOP[i] >= chr8_start

    univ_rows[[length(univ_rows)+1L]] <- u

    # bivariate local rg只有在两个目标表型：
    # (1) 均被process.locus保留；
    # (2) 均完成univariate分析；
    # (3) 均通过预设local-h2 Bonferroni阈值
    # 时才允许计算。
    can_bivar <- FALSE
    if (n_available == 2L &&
        all(traits %in% available_phenos) &&
        all(traits %in% u$phen)) {
      pass <- setNames(as.logical(u$passes_univ_bonf), as.character(u$phen))
      can_bivar <- all(!is.na(pass[traits])) && all(pass[traits])
    }

    if (isTRUE(can_bivar)) {
      # LAVA bivariate P/CI包含模拟步骤；按locus固定seed保证重复运行可复现。
      set.seed(250000L + as.integer(loci$LOC[i]))

      b <- tryCatch(
        LAVA::run.bivar(
          locus,
          phenos=traits
        ),
        error=function(e) e
      )

      if (inherits(b, "error")) {
        warnings_seen[[length(warnings_seen)+1L]] <- data.frame(
          LOC=loci$LOC[i], warning=paste0("run.bivar ERROR: ", conditionMessage(b)),
          stringsAsFactors=FALSE
        )
      } else if (!is.null(b) &&
                 nrow(b) == 1L &&
                 all(c("phen1","phen2","rho","rho.lower","rho.upper","p") %in% names(b))) {
        b$LOC <- loci$LOC[i]
        b$CHR <- loci$CHR[i]
        b$START <- loci$START[i]
        b$STOP <- loci$STOP[i]
        b$n_snps <- locus$n.snps
        b$n_pcs <- locus$K
        b$overlaps_known_chr8 <- loci$CHR[i] == 8L &&
          loci$START[i] <= chr8_end && loci$STOP[i] >= chr8_start
        bivar_rows[[length(bivar_rows)+1L]] <- b
      } else {
        warnings_seen[[length(warnings_seen)+1L]] <- data.frame(
          LOC=loci$LOC[i], warning="run.bivar returned unexpected structure",
          stringsAsFactors=FALSE
        )
      }
    }

    # 每50个locus保存一次检查点。
    if (i %% 50L == 0L) {
      saveRDS(
        list(last_i=i, univ=univ_rows, bivar=bivar_rows,
             locus_qc=locus_qc, warnings=warnings_seen),
        file.path(checkpoints, "analysis_checkpoint.rds")
      )
    }
  }

  # -------- 11. 汇总及多重校正 --------
  univ <- if (length(univ_rows))
    data.table::rbindlist(univ_rows, fill=TRUE) else data.table::data.table()
  bivar <- if (length(bivar_rows))
    data.table::rbindlist(bivar_rows, fill=TRUE) else data.table::data.table()
  qc <- if (length(locus_qc))
    data.table::rbindlist(locus_qc, fill=TRUE) else data.table::data.table()
  warns <- if (length(warnings_seen))
    data.table::rbindlist(warnings_seen, fill=TRUE) else
      data.table::data.table(LOC=character(), warning=character())

  data.table::fwrite(univ, file.path(reports, "LAVA_local_h2_all.tsv"), sep="\t", na="NA")
  data.table::fwrite(qc, file.path(reports, "LAVA_locus_QC.tsv"), sep="\t", na="NA")
  data.table::fwrite(warns, file.path(reports, "LAVA_warnings.tsv"), sep="\t", na="NA")

  n_biv <- nrow(bivar)
  if (n_biv > 0L) {
    if (!all(c("rho","rho.lower","rho.upper","p") %in% names(bivar)))
      stop("LAVA bivariate输出缺少预期字段。")
    bivar[, FDR_BH_tested := p.adjust(p, method="BH", n=n_biv)]
    bivar[, Bonferroni_tested := p.adjust(p, method="bonferroni", n=n_biv)]
    bivar[, Bonferroni_2495 := pmin(1, p * nrow(loci))]
    bivar[, significant_FDR := FDR_BH_tested < 0.05]
    bivar[, significant_Bonf_tested := Bonferroni_tested < 0.05]
    bivar[, significant_Bonf_2495 := Bonferroni_2495 < 0.05]
    data.table::setorder(bivar, p)
    data.table::fwrite(
      bivar, file.path(reports, "LAVA_local_rg_tested.tsv"),
      sep="\t", na="NA"
    )
    data.table::fwrite(
      bivar[significant_Bonf_tested == TRUE],
      file.path(reports, "LAVA_local_rg_Bonferroni_significant.tsv"),
      sep="\t", na="NA"
    )
    data.table::fwrite(
      bivar[CHR != 8],
      file.path(reports, "LAVA_local_rg_outside_chr8.tsv"),
      sep="\t", na="NA"
    )
    data.table::fwrite(
      bivar[overlaps_known_chr8 == TRUE],
      file.path(reports, "LAVA_chr8_local_rg.tsv"),
      sep="\t", na="NA"
    )
  } else {
    data.table::fwrite(
      data.table::data.table(note="No loci passed the prespecified dual local-h2 Bonferroni screen; no bivariate rg was estimated."),
      file.path(reports, "LAVA_local_rg_tested.tsv"), sep="\t"
    )
  }

  # 每个locus两个表型local-h2是否共同过阈值的审计。
  # 某一表型若被process.locus因负方差删除，则其pass明确记为FALSE，而不是NA。
  dual_pass_count <- 0L
  if (nrow(univ)) {
    wide_pass <- data.table::dcast(
      univ[, .(LOC, CHR, START, STOP, phen, p, h2.obs, passes_univ_bonf)],
      LOC + CHR + START + STOP ~ phen,
      value.var=c("p","h2.obs","passes_univ_bonf")
    )

    required_pass_cols <- paste0("passes_univ_bonf_", traits)
    for (cc in required_pass_cols) {
      if (!cc %in% names(wide_pass)) wide_pass[, (cc) := FALSE]
      wide_pass[is.na(get(cc)), (cc) := FALSE]
    }

    wide_pass[, dual_pass :=
      get("passes_univ_bonf_pancreas_fat") &
      get("passes_univ_bonf_pancreas_iron")
    ]
    dual_pass_count <- sum(wide_pass$dual_pass, na.rm=TRUE)

    data.table::fwrite(
      wide_pass,
      file.path(reports, "LAVA_dual_local_h2_screen.tsv"),
      sep="\t", na="NA"
    )
  }

  # -------- 12. QC摘要与解释 --------
  summary <- data.frame(
    metric = c(
      "n_partition_loci",
      "n_loci_analysable_both_phenotypes",
      "n_loci_analysable_one_phenotype",
      "n_loci_expected_all_negative_variance",
      "n_loci_other_expected_nonanalysable",
      "n_loci_unclassified_or_R_exception",
      "univ_threshold",
      "n_loci_with_univ_results",
      "n_process_locus_with_both_target_phenotypes",
      "n_process_locus_with_one_target_phenotype",
      "n_loci_both_traits_pass_local_h2",
      "n_bivariate_local_rg_tests",
      "n_bivariate_FDR_significant",
      "n_bivariate_Bonferroni_tested_significant",
      "n_bivariate_Bonferroni_2495_significant",
      "n_tested_local_rg_outside_chr8",
      "n_Bonferroni_significant_outside_chr8"
    ),
    value = c(
      nrow(loci),
      sum(qc$process_status=="ANALYSABLE_BOTH_PHENOTYPES", na.rm=TRUE),
      sum(qc$process_status=="ANALYSABLE_ONE_PHENOTYPE", na.rm=TRUE),
      sum(qc$process_status=="EXPECTED_ALL_NEGATIVE_VARIANCE", na.rm=TRUE),
      sum(qc$process_status %in% c(
        "EXPECTED_ALL_PHENOTYPES_FAILED","EXPECTED_TOO_FEW_SNPS","EXPECTED_TOO_FEW_PCS"
      ), na.rm=TRUE),
      sum(qc$process_status %in% c(
        "NULL_UNCLASSIFIED_REVIEW","R_EXCEPTION","RETURNED_NO_TARGET_PHENOTYPE_REVIEW"
      ), na.rm=TRUE),
      univ_thresh,
      length(unique(univ$LOC)),
      if ("n_available_phenos" %in% names(qc))
        sum(qc$n_available_phenos==2L, na.rm=TRUE) else NA_integer_,
      if ("n_available_phenos" %in% names(qc))
        sum(qc$n_available_phenos==1L, na.rm=TRUE) else NA_integer_,
      dual_pass_count,
      n_biv,
      if (n_biv) sum(bivar$significant_FDR, na.rm=TRUE) else 0L,
      if (n_biv) sum(bivar$significant_Bonf_tested, na.rm=TRUE) else 0L,
      if (n_biv) sum(bivar$significant_Bonf_2495, na.rm=TRUE) else 0L,
      if (n_biv) sum(bivar$CHR != 8) else 0L,
      if (n_biv) sum(bivar$CHR != 8 & bivar$significant_Bonf_tested, na.rm=TRUE) else 0L
    ),
    stringsAsFactors=FALSE
  )
  data.table::fwrite(summary, file.path(reports, "Step25_summary.tsv"), sep="\t")

  # chr8是否被LAVA partition覆盖并能测试。
  chr8_loci <- loci[loci$CHR==8 & loci$START<=chr8_end & loci$STOP>=chr8_start, , drop=FALSE]
  data.table::fwrite(chr8_loci, file.path(reports, "chr8_partition_context.tsv"), sep="\t")
  if (nrow(chr8_loci) < 1L)
    writeLines("REVIEW: no LAVA partition overlaps known chr8 interval",
               file.path(reports, "chr8_context_status.txt"))
  else
    writeLines("PASS: at least one LAVA partition overlaps known chr8 interval",
               file.path(reports, "chr8_context_status.txt"))

  # 轻量图：仅对实际测试的local rg画全基因组位置图。
  if (n_biv > 0L) {
    tryCatch({
      pdf(file.path(reports, "LAVA_local_rg_manhattan_style.pdf"), width=12, height=5)
      oldpar <- par(no.readonly=TRUE); on.exit(par(oldpar), add=TRUE)
      chr_lengths <- tapply(loci$STOP, loci$CHR, max)
      offsets <- c(0, cumsum(as.numeric(chr_lengths[-length(chr_lengths)])))
      names(offsets) <- names(chr_lengths)
      mid <- (bivar$START + bivar$STOP)/2
      x <- mid + offsets[as.character(bivar$CHR)]
      y <- -log10(pmax(bivar$p, .Machine$double.xmin))
      plot(x, y, pch=19, cex=.7, xaxt="n",
           xlab="Genome position", ylab="-log10(local rg P)",
           main="LAVA: pancreatic fat vs pancreatic iron")
      chr_mid <- offsets + as.numeric(chr_lengths)/2
      axis(1, at=chr_mid, labels=names(chr_lengths), cex.axis=.7)
      if (n_biv > 0) {
        abline(h=-log10(0.05/n_biv), lty=2)
      }
      dev.off()
    }, error=function(e) {
      writeLines(conditionMessage(e), file.path(reports, "plot_error.txt"))
    })
  }

  # -------- 13. 完成状态及复核包 --------
  writeLines(capture.output(sessionInfo()), file.path(reports, "R_sessionInfo_after.txt"))

  # 完整性条件：
  # 1) 2495个partition必须全部被尝试并各有一条QC记录；
  # 2) LAVA预期的负方差/局部不可分析结果不构成失败；
  # 3) 只有真实R异常或无法分类的NULL才要求人工复核。
  if (nrow(qc) != nrow(loci))
    stop("locus QC行数不等于2495，分析可能不完整。")

  if (anyDuplicated(qc$LOC))
    stop("locus QC出现重复LOC，分析记录异常。")

  unexpected_n <- sum(qc$process_status %in% c(
    "NULL_UNCLASSIFIED_REVIEW",
    "R_EXCEPTION",
    "RETURNED_NO_TARGET_PHENOTYPE_REVIEW"
  ), na.rm=TRUE)

  invalid_pheno_msg <- any(grepl(
    "Invalid phenotype ID",
    warns$warning,
    fixed=TRUE
  ))

  if (invalid_pheno_msg)
    stop("仍检测到Invalid phenotype ID，说明分析流程存在逻辑异常。")

  if (unexpected_n > 0L) {
    writeLines(
      paste0("REVIEW_UNEXPECTED_LOCUS_PROCESSING_EVENTS=", unexpected_n),
      file.path(reports, "unexpected_processing_review.txt")
    )
    stop("存在无法归类的process.locus/R异常，请上传当前结果目录或日志复核。")
  }

  writeLines(
    "COMPLETED_LAVA_LOCAL_RG_REVIEW_RESULTS",
    file.path(reports, "completion_status.txt")
  )

  # 仅打包小型报告，不打包15GB reference和原始sumstats。
  zipfile <- file.path(out, "25_LAVA_local_rg_review.zip")
  oldwd <- getwd()
  on.exit(setwd(oldwd), add=TRUE)
  setwd(out)
  files_to_zip <- c(
    list.files("reports", full.names=TRUE),
    list.files("logs", full.names=TRUE)
  )
  files_to_zip <- files_to_zip[file.exists(files_to_zip)]
  utils::zip(zipfile, files=files_to_zip, flags="-q")

  cat("\n===== 第25步完成 =====\n")
  print(summary)
  if (n_biv > 0L) {
    cat("\n按P值排序的前20个local rg结果：\n")
    print(head(bivar[, .(LOC,CHR,START,STOP,rho,rho.lower,rho.upper,p,
                         FDR_BH_tested,Bonferroni_tested,overlaps_known_chr8)], 20))
  } else {
    cat("\n没有locus同时通过两个表型的严格local-h2 Bonferroni筛选，因此未进行bivariate local rg。\n")
  }
  cat("\n结果目录：", out, "\n")
  cat("复核包：", zipfile, "\n")
  cat("请保留 ~/胰腺/lava_reference/UKB_v1.1；后续local analyses会直接复用。\n")
})()
