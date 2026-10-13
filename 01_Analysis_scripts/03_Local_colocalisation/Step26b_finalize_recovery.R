# 第26b步结果汇总恢复脚本
# 适用场景：chr19和chr8的SuSiE-RSS + coloc.susie均已成功，
# 但最终汇总因 duplicated column name [region_id] 报错。
# 本脚本不会重跑SuSiE，只做结果汇总、完成标记和重新打包。

(function() {
  if (!requireNamespace("data.table", quietly=TRUE)) stop("缺少data.table")
  suppressPackageStartupMessages(library(data.table))

  base <- path.expand("~/胰腺")
  parent <- file.path(base, "26b_coloc_susie_fix")
  dirs <- list.dirs(parent, recursive=FALSE, full.names=TRUE)
  if (!length(dirs)) stop("未找到第26b结果目录。")
  latest <- dirs[which.max(file.info(dirs)$mtime)]

  cat("使用第26b目录：\n", latest, "\n")

  reports <- file.path(latest, "reports")
  susiedir <- file.path(latest, "susie")
  dir.create(reports, recursive=TRUE, showWarnings=FALSE)

  regions <- data.table(
    region_id=c("chr19_LOC2340_primary","chr8_positive_control"),
    role=c("STEP25B_EXPLORATORY_CANDIDATE","POSITIVE_CONTROL")
  )

  status_rows <- list()

  for (rid in regions$region_id) {
    rdir <- file.path(susiedir, rid)
    qcf <- file.path(rdir, "region_QC.tsv")
    sumf <- file.path(rdir, "coloc_susie_summary.tsv")
    objf <- file.path(rdir, "Step26b_susie_objects.rds")

    if (!all(file.exists(c(qcf,sumf,objf)))) {
      status_rows[[rid]] <- data.table(
        region_id=rid, status="MISSING_OUTPUT_FILES",
        message=paste(c(qcf,sumf,objf)[!file.exists(c(qcf,sumf,objf))], collapse=";")
      )
      next
    }

    qc <- fread(qcf)
    if (nrow(qc)!=1L) {
      status_rows[[rid]] <- data.table(
        region_id=rid, status="REGION_QC_ROWCOUNT_UNEXPECTED", message=""
      )
      next
    }

    # 关键修复：去掉region_QC里已有的region_id，避免重复列名
    if ("region_id" %in% names(qc)) qc[, region_id := NULL]

    ss <- fread(sumf)
    h4cols <- grep("^PP\\.H4", names(ss), value=TRUE)
    maxh4 <- NA_real_
    if (length(h4cols) && nrow(ss)>0) {
      v <- suppressWarnings(as.numeric(unlist(ss[, ..h4cols])))
      v <- v[is.finite(v)]
      if (length(v)) maxh4 <- max(v)
    }

    row <- data.table(region_id=rid, status="OK", message="")
    for (nm in names(qc)) {
      if (!nm %in% names(row)) row[[nm]] <- qc[[nm]][1]
    }
    row[, n_coloc_susie_signal_pairs := nrow(ss)]
    row[, max_PP_H4_susie := maxh4]

    status_rows[[rid]] <- row
  }

  status_tab <- rbindlist(status_rows, fill=TRUE)
  status_tab[, role := regions$role[match(region_id, regions$region_id)]]
  setcolorder(status_tab, c("region_id","role","status","message",
                            setdiff(names(status_tab),c("region_id","role","status","message"))))

  fwrite(status_tab, file.path(reports,"Step26b_susie_status.tsv"),
         sep="\t", na="NA")

  # 定位第26步ABF结果
  context_file <- file.path(reports,"analysis_context.txt")
  step26 <- NA_character_
  if (file.exists(context_file)) {
    x <- readLines(context_file, warn=FALSE)
    hit <- grep("^Step26_source_dir\\t", x, value=TRUE)
    if (length(hit)) step26 <- sub("^Step26_source_dir\\t","",hit[1])
  }

  if (is.na(step26) || !dir.exists(step26)) {
    p26 <- file.path(base,"26_cross_trait_colocalization")
    d26 <- list.dirs(p26, recursive=FALSE, full.names=TRUE)
    good <- d26[vapply(d26,function(d) all(file.exists(c(
      file.path(d,"coloc_abf","Step26_chr19_primary_coloc.tsv"),
      file.path(d,"coloc_abf","Step26_chr8_positive_control_coloc.tsv")
    ))), logical(1))]
    if (!length(good)) stop("无法定位第26步ABF结果。")
    step26 <- good[which.max(file.info(good)$mtime)]
  }

  a19 <- fread(file.path(step26,"coloc_abf","Step26_chr19_primary_coloc.tsv"))
  a8  <- fread(file.path(step26,"coloc_abf","Step26_chr8_positive_control_coloc.tsv"))

  gv <- function(dt,nm) if (nm %in% names(dt)) as.numeric(dt[[nm]][1]) else NA_real_

  abf_tab <- data.table(
    region_id=c("chr19_LOC2340_primary","chr8_positive_control"),
    ABF_PP_H3=c(gv(a19,"PP.H3"),gv(a8,"PP.H3")),
    ABF_PP_H4=c(gv(a19,"PP.H4"),gv(a8,"PP.H4"))
  )

  final_summary <- merge(abf_tab, status_tab, by="region_id", all.x=TRUE, sort=FALSE)
  fwrite(final_summary,
         file.path(reports,"Step26b_ABF_plus_SuSiE_summary.tsv"),
         sep="\t",na="NA")

  all_ok <- nrow(status_tab)==2L && all(status_tab$status=="OK")
  completion <- if (all_ok) {
    "COMPLETED_STEP26B_SUSIE_COLOCALIZATION"
  } else {
    "COMPLETED_STEP26B_WITH_REVIEW_FLAGS"
  }

  writeLines(completion, file.path(reports,"completion_status.txt"))

  writeLines(c(
    "Step26b recovery after final-table construction error.",
    "Both regional SuSiE-RSS and coloc.susie computations had already completed.",
    "The original failure was due only to a duplicated region_id column in post-processing.",
    "No SuSiE or coloc.susie model was rerun in this recovery step."
  ), file.path(reports,"Step26b_recovery_note.txt"))

  zipfile <- file.path(latest,"26b_coloc_susie_review.zip")
  if (file.exists(zipfile)) unlink(zipfile)

  oldwd <- getwd()
  on.exit(setwd(oldwd), add=TRUE)
  setwd(latest)

  zfiles <- c(
    list.files("reports",recursive=TRUE,full.names=TRUE),
    list.files("susie",recursive=TRUE,full.names=TRUE),
    list.files("logs",recursive=TRUE,full.names=TRUE)
  )
  zfiles <- zfiles[file.exists(zfiles) & !dir.exists(zfiles)]

  utils::zip(zipfile, files=zfiles, flags="-q")

  cat("\n===== 第26b恢复完成 =====\n")
  print(final_summary)
  cat("\n完成状态：", completion, "\n")
  cat("复核包：", zipfile, "\n")
})()
