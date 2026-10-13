# 第25b步：LAVA局部遗传相关敏感性分析
# pancreas fat vs pancreas iron
#
# 目的：
# 在第25步严格local-h2 Bonferroni筛选（0.05/2495）下无locus进入bivariate local rg后，
# 采用预先定义的探索性/敏感性阈值 P_local-h2 < 0.05，
# 识别两种表型均具有名义显著局部遗传力的区域，并对这些区域计算LAVA bivariate local rg。
#
# 解释边界：
# - 本步为 sensitivity / exploratory analysis，不替代第25步严格主分析；
# - 只有两个表型在同一locus均满足 local h2 P<0.05 才进入bivariate local rg；
# - 对实际完成的local rg检验统一计算BH-FDR和Bonferroni；
# - 同时保留更保守的2495-locus Bonferroni值；
# - local rg ≠ causality；local rg ≠ colocalization；
# - chr8仅作为既有重点区域标记，重点关注chr8之外是否存在探索性局部共享信号。

(function() {

  if (!requireNamespace("data.table", quietly = TRUE))
    stop("缺少R包：data.table")

  base <- path.expand("~/胰腺")

  # 1. 固定LAVA环境
  lava_lib <- path.expand("~/R_pancreas_LAVA_015")
  .libPaths(c(lava_lib, .libPaths()))

  if (!requireNamespace("LAVA", quietly = TRUE))
    stop("未检测到LAVA，请确认第25步已经成功安装LAVA v0.1.5。")

  if (as.character(utils::packageVersion("LAVA")) != "0.1.5")
    stop("要求LAVA v0.1.5；当前版本：",
         as.character(utils::packageVersion("LAVA")))

  suppressPackageStartupMessages(library(LAVA))
  suppressPackageStartupMessages(library(data.table))

  cat("Running LAVA version ",
      as.character(utils::packageVersion("LAVA")), "\n", sep="")

  # 2. 自动定位最新完成的第25步主分析目录
  step25_parent <- file.path(base, "25_LAVA_local_rg")

  if (!dir.exists(step25_parent))
    stop("未找到第25步目录：", step25_parent)

  dirs25 <- list.dirs(
    step25_parent,
    recursive = FALSE,
    full.names = TRUE
  )

  if (!length(dirs25))
    stop("第25步目录下没有结果批次。")

  completed25 <- dirs25[vapply(dirs25, function(d) {
    f <- file.path(d, "reports", "completion_status.txt")
    file.exists(f) &&
      identical(
        readLines(f, warn=FALSE)[1],
        "COMPLETED_LAVA_LOCAL_RG_REVIEW_RESULTS"
      )
  }, logical(1))]

  if (!length(completed25)) {
    stop(
      "没有找到已完成的第25步主分析目录。\n",
      "请先完成第25步最终QC并写入 COMPLETED_LAVA_LOCAL_RG_REVIEW_RESULTS。"
    )
  }

  latest25 <- completed25[
    which.max(file.info(completed25)$mtime)
  ]

  cat("使用第25步主分析目录：\n", latest25, "\n")

  # 3. 核查第25步核心文件
  report25 <- file.path(latest25, "reports")

  f_univ <- file.path(report25, "LAVA_local_h2_all.tsv")
  f_qc   <- file.path(report25, "LAVA_locus_QC.tsv")
  f_loci <- file.path(report25, "LAVA_2495_loci.tsv")
  f_input_info <- file.path(latest25, "LAVA_input_info.tsv")
  f_overlap <- file.path(latest25, "LAVA_sample_overlap.txt")

  required <- c(
    f_univ, f_qc, f_loci, f_input_info, f_overlap
  )

  if (!all(file.exists(required))) {
    stop(
      "第25步缺少必要文件：\n",
      paste(required[!file.exists(required)], collapse="\n")
    )
  }

  univ25 <- fread(f_univ)
  qc25 <- fread(f_qc)
  loci <- fread(f_loci)

  if (nrow(loci) != 2495L)
    stop("第25步LAVA partition数量不是2495：", nrow(loci))

  traits <- c("pancreas_fat", "pancreas_iron")

  if (!all(c("LOC","CHR","START","STOP") %in% names(loci)))
    stop("LAVA loci字段异常。")

  if (!all(c("LOC","phen","p","h2.obs") %in% names(univ25)))
    stop("第25步local h2文件缺少预期字段。实际字段：",
         paste(names(univ25), collapse=", "))

  # 4. 输出目录
  parent <- file.path(base, "25b_LAVA_local_rg_sensitivity")
  dir.create(parent, recursive=TRUE, showWarnings=FALSE)

  out <- tempfile(
    "fat_iron_nominal_h2_",
    tmpdir=parent
  )

  if (!dir.create(out))
    stop("无法创建25b输出目录：", out)

  reports <- file.path(out, "reports")
  logs <- file.path(out, "logs")
  dir.create(reports)
  dir.create(logs)

  writeLines(
    "INCOMPLETE_CHECK_LOGS",
    file.path(reports, "completion_status.txt")
  )

  writeLines(
    c(
      paste0("Step25_source_dir\t", latest25),
      paste0("LAVA_version\t", as.character(utils::packageVersion("LAVA"))),
      "Sensitivity_univariate_threshold\t0.05",
      "Primary_step25_univariate_threshold\t0.05/2495",
      "Interpretation\tExploratory/sensitivity analysis; does not replace strict Step25 primary analysis."
    ),
    file.path(reports, "analysis_context.txt")
  )

  # 5. 验证LAVA reference
  lava_root <- file.path(base, "lava_reference")
  refroot <- file.path(lava_root, "UKB_v1.1")
  ref_prefix <- file.path(refroot, "lava-ukb-v1.1")

  expected_ref <- unlist(lapply(1:22, function(chr) {
    file.path(
      refroot,
      paste0(
        "lava-ukb-v1.1_chr",
        chr,
        c(".info", ".bcor")
      )
    )
  }))

  if (!all(file.exists(expected_ref) & file.info(expected_ref)$size > 0)) {
    stop(
      "LAVA UKB v1.1 reference不完整。\n",
      "本步不会重新下载reference，请先确认第25步reference完整。"
    )
  }

  # 6. 从第25步local-h2结果定义25b候选loci
  uu <- univ25[
    phen %in% traits &
      is.finite(p)
  ]

  if (!nrow(uu))
    stop("第25步local h2结果中没有目标表型数据。")

  nominal_wide <- dcast(
    uu[, .(
      LOC, CHR, START, STOP,
      phen,
      p,
      h2.obs
    )],
    LOC + CHR + START + STOP ~ phen,
    value.var = c("p","h2.obs")
  )

  p_fat_col <- "p_pancreas_fat"
  p_iron_col <- "p_pancreas_iron"

  if (!all(c(p_fat_col, p_iron_col) %in% names(nominal_wide))) {
    stop(
      "无法构建两个表型的wide local-h2结果。\n字段：",
      paste(names(nominal_wide), collapse=", ")
    )
  }

  nominal_wide[, fat_nominal :=
    is.finite(get(p_fat_col)) &
      get(p_fat_col) < 0.05
  ]

  nominal_wide[, iron_nominal :=
    is.finite(get(p_iron_col)) &
      get(p_iron_col) < 0.05
  ]

  nominal_wide[, dual_nominal :=
    fat_nominal & iron_nominal
  ]

  candidate <- nominal_wide[
    dual_nominal == TRUE
  ]

  setorder(candidate, CHR, START, STOP)

  fwrite(
    nominal_wide,
    file.path(
      reports,
      "Step25b_nominal_local_h2_screen_all.tsv"
    ),
    sep="\t",
    na="NA"
  )

  fwrite(
    candidate,
    file.path(
      reports,
      "Step25b_candidate_loci_dual_nominal_h2.tsv"
    ),
    sep="\t",
    na="NA"
  )

  cat(
    "\n两个表型均满足local h2 P<0.05的候选loci数量：",
    nrow(candidate),
    "\n"
  )

  # 没有候选也属于合法结果
  if (nrow(candidate) == 0L) {

    summary <- data.table(
      metric=c(
        "n_partition_loci",
        "n_candidate_loci_dual_nominal_h2",
        "n_bivariate_tests",
        "n_FDR_significant",
        "n_Bonferroni_significant"
      ),
      value=c(
        2495L,0L,0L,0L,0L
      )
    )

    fwrite(
      summary,
      file.path(reports, "Step25b_summary.tsv"),
      sep="\t"
    )

    writeLines(
      "COMPLETED_NO_DUAL_NOMINAL_H2_LOCI",
      file.path(reports, "completion_status.txt")
    )

    cat(
      "\n没有任何locus同时满足两个表型local h2 P<0.05；",
      "第25b步无需进行bivariate local rg。\n"
    )

    return(invisible(NULL))
  }

  # 7. 重新处理LAVA输入
  process_log <- file.path(
    logs,
    "process_input_console.txt"
  )

  zz <- file(process_log, open="wt")
  sink(zz, type="output")
  sink(zz, type="message")

  input <- tryCatch(
    LAVA::process.input(
      input.info.file=f_input_info,
      sample.overlap.file=f_overlap,
      ref.prefix=ref_prefix,
      phenos=traits
    ),
    error=function(e) e
  )

  sink(type="message")
  sink(type="output")
  close(zz)

  if (inherits(input, "error")) {
    stop(
      "LAVA process.input失败：",
      conditionMessage(input),
      "\n查看：",
      process_log
    )
  }

  # 8. 仅对25b候选loci重新process并run.bivar
  chr8_start <- 129047378L
  chr8_end   <- 130047378L

  bivar_rows <- list()
  qc_rows <- list()
  print_rows <- list()

  cat(
    "\n开始第25b步候选loci bivariate local rg分析，共",
    nrow(candidate),
    "个候选loci。\n"
  )

  for (j in seq_len(nrow(candidate))) {

    loc_id <- candidate$LOC[j]

    loc_row <- loci[
      LOC == loc_id
    ]

    if (nrow(loc_row) != 1L) {
      qc_rows[[length(qc_rows)+1L]] <- data.table(
        LOC=loc_id,
        status="LOCUS_LOOKUP_ERROR",
        message="LOC not uniquely found in 2495 partition file"
      )
      next
    }

    if (j %% 10L == 0L || j == 1L) {
      cat(
        "25b locus ",
        j, "/", nrow(candidate),
        " (LOC=", loc_id,
        ", chr", loc_row$CHR,
        ":", loc_row$START,
        "-", loc_row$STOP,
        ")\n",
        sep=""
      )
    }

    locus <- NULL
    err <- NULL
    lava_print <- character()

    lava_print <- capture.output({
      locus <- tryCatch(
        LAVA::process.locus(
          as.data.frame(loc_row),
          input
        ),
        error=function(e) {
          err <<- conditionMessage(e)
          NULL
        }
      )
    })

    if (length(lava_print)) {
      print_rows[[length(print_rows)+1L]] <- data.table(
        LOC=loc_id,
        lava_message=paste(lava_print, collapse=" | ")
      )
    }

    if (!is.null(err)) {
      qc_rows[[length(qc_rows)+1L]] <- data.table(
        LOC=loc_id,
        status="PROCESS_LOCUS_R_EXCEPTION",
        message=err
      )
      next
    }

    if (is.null(locus)) {
      qc_rows[[length(qc_rows)+1L]] <- data.table(
        LOC=loc_id,
        status="PROCESS_LOCUS_NULL_REVIEW",
        message=if (length(lava_print))
          paste(lava_print, collapse=" | ")
        else
          "process.locus returned NULL"
      )
      next
    }

    available <- intersect(
      traits,
      as.character(locus$phenos)
    )

    if (length(available) != 2L ||
        !all(traits %in% available)) {

      qc_rows[[length(qc_rows)+1L]] <- data.table(
        LOC=loc_id,
        status="TARGET_PHENOTYPE_DROPPED_ON_REPROCESS",
        message=paste(
          "available=",
          paste(available, collapse=";")
        )
      )
      next
    }

    # 重新运行univ确认两个表型仍满足P<0.05
    u <- tryCatch(
      LAVA::run.univ(
        locus,
        phenos=traits,
        var=TRUE,
        cap.estimates=FALSE
      ),
      error=function(e) e
    )

    if (inherits(u, "error")) {
      qc_rows[[length(qc_rows)+1L]] <- data.table(
        LOC=loc_id,
        status="RUN_UNIV_EXCEPTION",
        message=conditionMessage(u)
      )
      next
    }

    if (is.null(u) ||
        nrow(u) != 2L ||
        !all(traits %in% u$phen)) {

      qc_rows[[length(qc_rows)+1L]] <- data.table(
        LOC=loc_id,
        status="RUN_UNIV_UNEXPECTED",
        message="Unexpected univariate result structure"
      )
      next
    }

    p_now <- setNames(
      u$p,
      u$phen
    )

    if (!all(is.finite(p_now[traits])) ||
        !all(p_now[traits] < 0.05)) {

      qc_rows[[length(qc_rows)+1L]] <- data.table(
        LOC=loc_id,
        status="NO_LONGER_DUAL_NOMINAL_ON_REPROCESS",
        message=paste0(
          "P_fat=",
          signif(p_now["pancreas_fat"],6),
          ";P_iron=",
          signif(p_now["pancreas_iron"],6)
        )
      )
      next
    }

    set.seed(
      250000L + as.integer(loc_id)
    )

    b <- tryCatch(
      LAVA::run.bivar(
        locus,
        phenos=traits
      ),
      error=function(e) e
    )

    if (inherits(b, "error")) {
      qc_rows[[length(qc_rows)+1L]] <- data.table(
        LOC=loc_id,
        status="RUN_BIVAR_EXCEPTION",
        message=conditionMessage(b)
      )
      next
    }

    if (is.null(b) || nrow(b) != 1L) {
      qc_rows[[length(qc_rows)+1L]] <- data.table(
        LOC=loc_id,
        status="RUN_BIVAR_UNEXPECTED",
        message="Unexpected bivariate result structure"
      )
      next
    }

    b <- as.data.table(b)

    b[, LOC := loc_id]
    b[, CHR := loc_row$CHR]
    b[, START := loc_row$START]
    b[, STOP := loc_row$STOP]
    b[, local_h2_p_pancreas_fat :=
        as.numeric(p_now["pancreas_fat"])]
    b[, local_h2_p_pancreas_iron :=
        as.numeric(p_now["pancreas_iron"])]
    b[, n_snps := locus$n.snps]
    b[, n_pcs := locus$K]

    b[, overlaps_known_chr8 :=
      CHR == 8L &
      START <= chr8_end &
      STOP >= chr8_start
    ]

    bivar_rows[[length(bivar_rows)+1L]] <- b

    qc_rows[[length(qc_rows)+1L]] <- data.table(
      LOC=loc_id,
      status="BIVARIATE_COMPLETED",
      message=""
    )
  }

  # 9. 汇总QC
  qc <- if (length(qc_rows))
    rbindlist(qc_rows, fill=TRUE)
  else
    data.table(
      LOC=integer(),
      status=character(),
      message=character()
    )

  fwrite(
    qc,
    file.path(reports, "Step25b_locus_QC.tsv"),
    sep="\t",
    na="NA"
  )

  if (length(print_rows)) {
    pr <- rbindlist(print_rows, fill=TRUE)
  } else {
    pr <- data.table(
      LOC=integer(),
      lava_message=character()
    )
  }

  fwrite(
    pr,
    file.path(reports, "Step25b_LAVA_messages.tsv"),
    sep="\t",
    na="NA"
  )

  # 10. 汇总bivariate结果与多重校正
  bivar <- if (length(bivar_rows))
    rbindlist(bivar_rows, fill=TRUE)
  else
    data.table()

  n_biv <- nrow(bivar)

  if (n_biv > 0L) {

    if (!all(c(
      "rho",
      "rho.lower",
      "rho.upper",
      "p"
    ) %in% names(bivar))) {
      stop(
        "LAVA bivariate输出缺少预期字段：",
        paste(names(bivar), collapse=", ")
      )
    }

    bivar[, FDR_BH_tested :=
      p.adjust(
        p,
        method="BH",
        n=n_biv
      )
    ]

    bivar[, Bonferroni_tested :=
      p.adjust(
        p,
        method="bonferroni",
        n=n_biv
      )
    ]

    bivar[, Bonferroni_2495 :=
      pmin(
        1,
        p * 2495
      )
    ]

    bivar[, significant_FDR :=
      FDR_BH_tested < 0.05
    ]

    bivar[, significant_Bonferroni_tested :=
      Bonferroni_tested < 0.05
    ]

    bivar[, significant_Bonferroni_2495 :=
      Bonferroni_2495 < 0.05
    ]

    setorder(
      bivar,
      p
    )

    fwrite(
      bivar,
      file.path(
        reports,
        "Step25b_LAVA_local_rg_all_tested.tsv"
      ),
      sep="\t",
      na="NA"
    )

    fwrite(
      bivar[
        significant_FDR == TRUE
      ],
      file.path(
        reports,
        "Step25b_LAVA_local_rg_FDR_significant.tsv"
      ),
      sep="\t",
      na="NA"
    )

    fwrite(
      bivar[
        significant_Bonferroni_tested == TRUE
      ],
      file.path(
        reports,
        "Step25b_LAVA_local_rg_Bonferroni_significant.tsv"
      ),
      sep="\t",
      na="NA"
    )

    fwrite(
      bivar[
        CHR != 8
      ],
      file.path(
        reports,
        "Step25b_LAVA_local_rg_outside_chr8.tsv"
      ),
      sep="\t",
      na="NA"
    )

    fwrite(
      bivar[
        overlaps_known_chr8 == TRUE
      ],
      file.path(
        reports,
        "Step25b_LAVA_chr8_local_rg.tsv"
      ),
      sep="\t",
      na="NA"
    )

  } else {

    fwrite(
      data.table(
        note=
          "No bivariate local rg test completed after reprocessing nominal dual-local-h2 loci."
      ),
      file.path(
        reports,
        "Step25b_LAVA_local_rg_all_tested.tsv"
      ),
      sep="\t"
    )
  }

  # 11. 汇总统计
  summary <- data.table(
    metric=c(
      "n_partition_loci",
      "n_loci_with_local_h2_results_step25",
      "n_candidate_loci_dual_nominal_h2_P_lt_0.05",
      "n_candidate_loci_reprocessed",
      "n_bivariate_local_rg_completed",
      "n_FDR_significant",
      "n_Bonferroni_tested_significant",
      "n_Bonferroni_2495_significant",
      "n_tested_outside_chr8",
      "n_FDR_significant_outside_chr8",
      "n_Bonferroni_tested_significant_outside_chr8",
      "n_tested_overlapping_known_chr8"
    ),
    value=c(
      2495L,
      length(unique(uu$LOC)),
      nrow(candidate),
      nrow(candidate),
      n_biv,
      if (n_biv) sum(bivar$significant_FDR, na.rm=TRUE) else 0L,
      if (n_biv) sum(bivar$significant_Bonferroni_tested, na.rm=TRUE) else 0L,
      if (n_biv) sum(bivar$significant_Bonferroni_2495, na.rm=TRUE) else 0L,
      if (n_biv) sum(bivar$CHR != 8) else 0L,
      if (n_biv) sum(
        bivar$CHR != 8 &
        bivar$significant_FDR,
        na.rm=TRUE
      ) else 0L,
      if (n_biv) sum(
        bivar$CHR != 8 &
        bivar$significant_Bonferroni_tested,
        na.rm=TRUE
      ) else 0L,
      if (n_biv) sum(
        bivar$overlaps_known_chr8,
        na.rm=TRUE
      ) else 0L
    )
  )

  fwrite(
    summary,
    file.path(
      reports,
      "Step25b_summary.tsv"
    ),
    sep="\t"
  )

  # 12. 解释边界
  writeLines(
    c(
      "Step25b is a sensitivity/exploratory LAVA analysis.",
      "The strict Step25 primary analysis used local-h2 P < 0.05/2495 and yielded zero loci eligible for bivariate local rg.",
      "Step25b relaxes the local-h2 inclusion threshold to nominal P < 0.05 for BOTH pancreatic fat and pancreatic iron in the same locus.",
      "Only loci satisfying this dual nominal local-h2 criterion are tested for local rg.",
      "BH-FDR and Bonferroni corrections are applied across all actually completed Step25b local-rg tests.",
      "Bonferroni_2495 is additionally retained as a deliberately conservative genome-wide reference correction.",
      "Step25b findings must be labelled exploratory/sensitivity results and must not replace the strict Step25 primary analysis.",
      "Local genetic correlation does not establish causality or colocalization.",
      "Any locus prioritized for mechanistic interpretation should undergo separate cross-trait colocalization/fine-mapping."
    ),
    file.path(
      reports,
      "Step25b_methods_and_interpretation.txt"
    )
  )

  # 13. 简单可视化
  if (n_biv > 0L) {

    tryCatch({

      pdf(
        file.path(
          reports,
          "Step25b_local_rg_genome_plot.pdf"
        ),
        width=12,
        height=5
      )

      chr_lengths <- tapply(
        loci$STOP,
        loci$CHR,
        max
      )

      offsets <- c(
        0,
        cumsum(
          as.numeric(
            chr_lengths[
              -length(chr_lengths)
            ]
          )
        )
      )

      names(offsets) <- names(chr_lengths)

      mid <- (
        bivar$START +
          bivar$STOP
      ) / 2

      x <- mid +
        offsets[
          as.character(
            bivar$CHR
          )
        ]

      y <- -log10(
        pmax(
          bivar$p,
          .Machine$double.xmin
        )
      )

      plot(
        x,
        y,
        pch=19,
        cex=.8,
        xaxt="n",
        xlab="Genome position",
        ylab="-log10(local rg P)",
        main=
          "Step25b LAVA local rg: pancreatic fat vs pancreatic iron"
      )

      chr_mid <- offsets +
        as.numeric(chr_lengths)/2

      axis(
        1,
        at=chr_mid,
        labels=names(chr_lengths),
        cex.axis=.7
      )

      abline(
        h=-log10(
          0.05/n_biv
        ),
        lty=2
      )

      dev.off()

    }, error=function(e) {
      writeLines(
        conditionMessage(e),
        file.path(
          reports,
          "Step25b_plot_error.txt"
        )
      )
    })
  }

  # 14. 完整性检查
  serious_status <- c(
    "LOCUS_LOOKUP_ERROR",
    "PROCESS_LOCUS_R_EXCEPTION",
    "PROCESS_LOCUS_NULL_REVIEW",
    "TARGET_PHENOTYPE_DROPPED_ON_REPROCESS",
    "RUN_UNIV_EXCEPTION",
    "RUN_UNIV_UNEXPECTED",
    "NO_LONGER_DUAL_NOMINAL_ON_REPROCESS",
    "RUN_BIVAR_EXCEPTION",
    "RUN_BIVAR_UNEXPECTED"
  )

  n_serious <- qc[
    status %in% serious_status,
    .N
  ]

  if (nrow(qc) != nrow(candidate))
    stop(
      "25b candidate QC行数与候选locus数不一致：",
      nrow(qc), " vs ", nrow(candidate)
    )

  writeLines(
    "COMPLETED_STEP25B_LAVA_NOMINAL_H2_SENSITIVITY",
    file.path(
      reports,
      "completion_status.txt"
    )
  )

  # 15. 打包复核文件
  zipfile <- file.path(
    out,
    "25b_LAVA_local_rg_review.zip"
  )

  oldwd <- getwd()
  on.exit(
    setwd(oldwd),
    add=TRUE
  )

  setwd(out)

  files_to_zip <- c(
    list.files(
      "reports",
      full.names=TRUE
    ),
    list.files(
      "logs",
      full.names=TRUE
    )
  )

  files_to_zip <- files_to_zip[
    file.exists(files_to_zip)
  ]

  utils::zip(
    zipfile,
    files=files_to_zip,
    flags="-q"
  )

  # 16. 控制台摘要
  cat("\n===== 第25b步完成 =====\n")
  print(summary)

  cat(
    "\n候选loci：",
    nrow(candidate),
    "\n"
  )

  cat(
    "完成bivariate local rg：",
    n_biv,
    "\n"
  )

  cat(
    "候选重处理过程中出现需复核状态：",
    n_serious,
    "\n"
  )

  if (n_biv > 0L) {

    cat(
      "\n按P值排序的前20个local rg结果：\n"
    )

    keep_cols <- intersect(
      c(
        "LOC","CHR","START","STOP",
        "rho","rho.lower","rho.upper","p",
        "FDR_BH_tested",
        "Bonferroni_tested",
        "Bonferroni_2495",
        "overlaps_known_chr8"
      ),
      names(bivar)
    )

    print(
      head(
        bivar[
          ,
          ..keep_cols
        ],
        20
      )
    )
  }

  cat(
    "\n结果目录：",
    out,
    "\n"
  )

  cat(
    "复核包：",
    zipfile,
    "\n"
  )

})()
