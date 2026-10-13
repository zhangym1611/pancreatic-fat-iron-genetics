# ============================================================
# 第27b步：PEPD / rs33824 / rs430752 胰腺eQTL验证与GWAS-eQTL共定位
# 修正版v6：保留v5修复；修正run_coloc中的window_id作用域错误，并静默读取首行时预期的curl broken-pipe提示。
#
# 研究问题：
# Step27将chr19 LOC2340的共享credible-set信号收敛到PEPD内的
# rs33824和rs430752。本步进一步回答：
# 1) PEPD在现有pancreas / pancreatic-islet eQTL数据中是否可评估？
# 2) rs33824 / rs430752在这些eQTL数据中的覆盖和效应如何？
# 3) PEPD expression QTL是否与pancreas fat / pancreas iron GWAS共享局部信号？
#
# 预设eQTL Catalogue数据集（与此前chr8 Step16/17保持一致）：
#   QTD000296 : GTEx pancreas
#   QTD000554 : pancreatic islet (van de Bunt)
#   QTD000574 : pancreatic/islet dataset used in prior Step16/17
# 实际study/sample/tissue/sample_size全部从官方tabix metadata动态读取，
# 不依赖手工写死标签。
#
# 数据源：
# eQTL Catalogue官方tabix metadata：
# https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv
#
# 重要原则：
# - eQTL Catalogue summary statistics为GRCh38；
# - MRI GWAS为GRCh37；
# - 本步不通过坐标直接合并，而通过唯一rsID + allele compatibility进行跨build匹配；
# - eQTL ALT为effect allele；
# - coloc.abf回答共享因果信号模型，不回答因果方向；
# - chr19是Step25b探索性候选，因此本步属于机制跟进而非独立复制。
# ============================================================

(function() {

  options(stringsAsFactors = FALSE)
  base <- path.expand("~/胰腺")

  # ------------------------------------------------------------
  # 0. R环境
  # ------------------------------------------------------------
  lib <- path.expand("~/R_pancreas_step27b")
  dir.create(lib, recursive=TRUE, showWarnings=FALSE)
  .libPaths(c(lib, .libPaths()))

  need <- function(pkg) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
      install.packages(pkg, lib=lib, repos="https://cloud.r-project.org")
    }
    if (!requireNamespace(pkg, quietly=TRUE)) stop("无法加载R包：", pkg)
  }

  for (p in c("data.table","coloc","httr","jsonlite")) need(p)

  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(httr))
  suppressPackageStartupMessages(library(jsonlite))

  cat("coloc version: ", as.character(packageVersion("coloc")), "\n", sep="")

  # ------------------------------------------------------------
  # 1. 外部命令：tabix
  #    v3：若服务器没有tabix/conda/mamba，则自动bootstrap micromamba，
  #    再在用户目录安装htslib 1.24。无需sudo，不污染既有LDSC/LAVA环境。
  # ------------------------------------------------------------
  find_tabix <- function() {
    cand <- unique(c(
      Sys.which("tabix"),
      path.expand("~/mr_tools/step27b_htslib/bin/tabix"),
      path.expand("~/miniconda3/bin/tabix"),
      path.expand("~/anaconda3/bin/tabix"),
      "/usr/bin/tabix",
      "/usr/local/bin/tabix"
    ))
    cand <- cand[nzchar(cand) & file.exists(cand)]
    if (length(cand)) cand[1] else ""
  }

  tabix_bin <- find_tabix()

  if (!nzchar(tabix_bin)) {

    cat(
      "\n未检测到tabix；开始在用户目录自举micromamba并安装htslib 1.24（无需sudo）...\n"
    )

    tools_root <- path.expand("~/mr_tools")
    bootstrap_dir <- file.path(tools_root, "micromamba_bootstrap")
    mamba_root <- file.path(tools_root, "micromamba_root")
    envdir <- file.path(tools_root, "step27b_htslib")

    dir.create(tools_root, recursive=TRUE, showWarnings=FALSE)
    dir.create(bootstrap_dir, recursive=TRUE, showWarnings=FALSE)
    dir.create(mamba_root, recursive=TRUE, showWarnings=FALSE)

    micromamba_bin <- file.path(bootstrap_dir, "bin", "micromamba")

    install_log <- file.path(
      base,
      "step27b_tabix_install.log"
    )

    # ----------------------------------------------------------
    # 1A. 如果没有micromamba，直接从官方micromamba endpoint下载
    # ----------------------------------------------------------
    if (!file.exists(micromamba_bin)) {

      archive <- file.path(
        bootstrap_dir,
        "micromamba-linux-64.tar.bz2"
      )

      curl_bin <- Sys.which("curl")
      wget_bin <- Sys.which("wget")

      download_ok <- FALSE

      if (nzchar(curl_bin)) {

        cat("使用curl下载micromamba...\n")

        rc <- system2(
          curl_bin,
          c(
            "-L",
            "--fail",
            "--retry","3",
            "--connect-timeout","30",
            "-o", shQuote(archive),
            shQuote("https://micro.mamba.pm/api/micromamba/linux-64/latest")
          ),
          stdout=install_log,
          stderr=install_log
        )

        download_ok <- identical(as.integer(rc),0L) &&
          file.exists(archive) &&
          file.info(archive)$size > 100000

      } else if (nzchar(wget_bin)) {

        cat("使用wget下载micromamba...\n")

        rc <- system2(
          wget_bin,
          c(
            "-O", shQuote(archive),
            shQuote("https://micro.mamba.pm/api/micromamba/linux-64/latest")
          ),
          stdout=install_log,
          stderr=install_log
        )

        download_ok <- identical(as.integer(rc),0L) &&
          file.exists(archive) &&
          file.info(archive)$size > 100000
      }

      if (!download_ok) {
        stop(
          "无法下载micromamba。\n",
          "服务器已有curl/wget，但下载失败。\n",
          "请查看日志：", install_log
        )
      }

      cat("解压micromamba...\n")

      untar_status <- tryCatch(
        {
          utils::untar(
            archive,
            exdir=bootstrap_dir
          )
          0L
        },
        error=function(e) {
          write(
            paste0("untar error: ",conditionMessage(e)),
            file=install_log,
            append=TRUE
          )
          1L
        }
      )

      if (!identical(untar_status,0L) ||
          !file.exists(micromamba_bin)) {

        # 如果R的untar没有正确处理bz2，尝试系统tar
        tar_bin <- Sys.which("tar")

        if (nzchar(tar_bin)) {

          rc_tar <- system2(
            tar_bin,
            c(
              "-xjf",
              shQuote(archive),
              "-C",
              shQuote(bootstrap_dir)
            ),
            stdout=install_log,
            stderr=install_log
          )

        } else {
          rc_tar <- 1L
        }

        if (!identical(as.integer(rc_tar),0L) ||
            !file.exists(micromamba_bin)) {
          stop(
            "micromamba下载成功，但解压失败。\n",
            "请查看日志：", install_log
          )
        }
      }

      Sys.chmod(
        micromamba_bin,
        mode="0755"
      )
    }

    if (!file.exists(micromamba_bin)) {
      stop("micromamba bootstrap失败：", micromamba_bin)
    }

    cat("micromamba：", micromamba_bin, "\n")

    # ----------------------------------------------------------
    # 1B. 安装htslib/tabix
    # ----------------------------------------------------------
    target_tabix <- file.path(envdir, "bin", "tabix")

    if (!file.exists(target_tabix)) {

      cat(
        "使用micromamba创建独立htslib环境：",
        envdir,
        "\n"
      )

      rc_env <- system2(
        micromamba_bin,
        c(
          "create",
          "-y",
          "-p", shQuote(envdir),
          "-c", "conda-forge",
          "-c", "bioconda",
          "htslib=1.24"
        ),
        stdout=install_log,
        stderr=install_log,
        env=c(
          paste0(
            "MAMBA_ROOT_PREFIX=",
            mamba_root
          )
        )
      )

      if (!identical(as.integer(rc_env),0L) ||
          !file.exists(target_tabix)) {

        stop(
          "micromamba已安装，但htslib/tabix环境创建失败。\n",
          "请查看日志：", install_log
        )
      }
    }

    tabix_bin <- target_tabix

    # ----------------------------------------------------------
    # 1C. 最终功能检查
    # ----------------------------------------------------------
    version_out <- tryCatch(
      system2(
        tabix_bin,
        "--version",
        stdout=TRUE,
        stderr=TRUE
      ),
      error=function(e) character()
    )

    if (!length(version_out)) {
      stop("tabix文件存在，但无法执行：",tabix_bin)
    }

    cat(
      "tabix安装成功：",
      tabix_bin,
      "\n"
    )

    cat(
      "tabix版本：",
      paste(version_out,collapse=" "),
      "\n"
    )
  }

  cat("tabix: ", tabix_bin, "\n", sep="")

  # ------------------------------------------------------------
  # 2. 自动定位Step26区域GWAS数据
  # ------------------------------------------------------------
  p26 <- file.path(base, "26_cross_trait_colocalization")
  if (!dir.exists(p26)) stop("未找到Step26目录：", p26)

  d26 <- list.dirs(p26, recursive=FALSE, full.names=TRUE)
  good26 <- d26[vapply(d26, function(d) {
    all(file.exists(c(
      file.path(
        d, "regional_sumstats",
        "pancreas_fat__chr19_LOC2340_primary.tsv.gz"
      ),
      file.path(
        d, "regional_sumstats",
        "pancreas_iron__chr19_LOC2340_primary.tsv.gz"
      )
    )))
  }, logical(1))]

  if (!length(good26)) {
    stop("未找到包含chr19两个原始区域GWAS文件的Step26批次。")
  }

  step26 <- good26[which.max(file.info(good26)$mtime)]
  cat("使用Step26目录：\n", step26, "\n")

  # ------------------------------------------------------------
  # 3. 自动定位Step27，确认PEPD结论来源
  # ------------------------------------------------------------
  p27 <- file.path(base, "27_functional_annotation")
  if (!dir.exists(p27)) stop("未找到Step27目录：", p27)

  d27 <- list.dirs(p27, recursive=FALSE, full.names=TRUE)
  good27 <- d27[vapply(d27, function(d) {
    f <- file.path(d, "reports", "completion_status.txt")
    file.exists(f) &&
      identical(
        readLines(f, warn=FALSE)[1],
        "COMPLETED_STEP27_FUNCTIONAL_ANNOTATION"
      )
  }, logical(1))]

  if (!length(good27)) {
    stop("未找到已完成的Step27目录。")
  }

  step27 <- good27[which.max(file.info(good27)$mtime)]
  priority_file <- file.path(
    step27, "tables", "Step27_priority_shared_candidates.tsv"
  )
  if (!file.exists(priority_file)) {
    stop("Step27缺少优先候选表：", priority_file)
  }

  priority27 <- fread(priority_file)

  # 必须能核对rs33824和rs430752
  expected_chr19 <- c("rs33824","rs430752")
  missing_expected <- setdiff(
    expected_chr19,
    priority27[region_id=="chr19_LOC2340_primary", unique(SNP)]
  )
  if (length(missing_expected)) {
    stop(
      "Step27优先候选表中缺少预期chr19共享SNP：",
      paste(missing_expected, collapse=",")
    )
  }

  # ------------------------------------------------------------
  # 4. 输出目录
  # ------------------------------------------------------------
  parent <- file.path(base, "27b_PEPD_eQTL_coloc")
  dir.create(parent, recursive=TRUE, showWarnings=FALSE)

  out <- tempfile("PEPD_pancreas_eQTL_", tmpdir=parent)
  if (!dir.create(out)) stop("无法创建Step27b输出目录。")

  reports <- file.path(out, "reports")
  rawdir  <- file.path(out, "raw_eqtl")
  tables  <- file.path(out, "tables")
  colocdir <- file.path(out, "coloc")
  logs <- file.path(out, "logs")

  for (d in c(reports,rawdir,tables,colocdir,logs))
    dir.create(d, recursive=TRUE, showWarnings=FALSE)

  writeLines(
    "INCOMPLETE_CHECK_LOGS",
    file.path(reports, "completion_status.txt")
  )

  # ------------------------------------------------------------
  # 5. Ensembl GRCh38：PEPD基因和两个候选SNP映射
  # ------------------------------------------------------------
  ens_server <- "https://rest.ensembl.org"

  ens_get <- function(path, tries=4L) {
    for (i in seq_len(tries)) {
      r <- tryCatch(
        GET(
          paste0(ens_server,path),
          add_headers(
            `Content-Type`="application/json",
            Accept="application/json"
          ),
          timeout(60)
        ),
        error=function(e) e
      )
      if (!inherits(r,"error") &&
          status_code(r)>=200 &&
          status_code(r)<300) {
        return(
          content(
            r,
            as="parsed",
            type="application/json",
            simplifyVector=FALSE
          )
        )
      }
      Sys.sleep(i)
    }
    NULL
  }

  pepd <- ens_get("/lookup/symbol/homo_sapiens/PEPD?expand=0")

  if (is.null(pepd) || is.null(pepd$id)) {
    stop("无法从Ensembl GRCh38解析PEPD。")
  }

  pepd_info <- data.table(
    gene_symbol="PEPD",
    gene_id=as.character(pepd$id),
    assembly=as.character(pepd$assembly_name),
    CHR=as.character(pepd$seq_region_name),
    START=as.numeric(pepd$start),
    STOP=as.numeric(pepd$end),
    strand=as.numeric(pepd$strand),
    biotype=as.character(pepd$biotype)
  )

  fwrite(
    pepd_info,
    file.path(tables, "Step27b_PEPD_GRCh38_gene_info.tsv"),
    sep="\t"
  )

  get_grch38_mapping <- function(rsid) {

    v <- ens_get(paste0("/variation/human/",rsid))

    if (is.null(v) || is.null(v$mappings)) {
      return(data.table(
        SNP=rsid,
        assembly="GRCh38",
        CHR=NA_character_,
        BP=NA_real_,
        allele_string=NA_character_,
        strand=NA_real_
      ))
    }

    mm <- v$mappings
    outm <- list()

    for (m in mm) {
      if (!is.null(m$assembly_name) &&
          m$assembly_name=="GRCh38" &&
          !is.null(m$seq_region_name) &&
          m$seq_region_name=="19") {

        outm[[length(outm)+1L]] <- data.table(
          SNP=rsid,
          assembly="GRCh38",
          CHR=as.character(m$seq_region_name),
          BP=as.numeric(m$start),
          allele_string=if (!is.null(m$allele_string))
            as.character(m$allele_string) else NA_character_,
          strand=if (!is.null(m$strand))
            as.numeric(m$strand) else NA_real_
        )
      }
    }

    if (!length(outm)) {
      return(data.table(
        SNP=rsid,
        assembly="GRCh38",
        CHR=NA_character_,
        BP=NA_real_,
        allele_string=NA_character_,
        strand=NA_real_
      ))
    }

    unique(rbindlist(outm, fill=TRUE))
  }

  map_shared <- rbindlist(
    lapply(expected_chr19, get_grch38_mapping),
    fill=TRUE
  )

  fwrite(
    map_shared,
    file.path(tables, "Step27b_rs33824_rs430752_GRCh38_mapping.tsv"),
    sep="\t",
    na="NA"
  )

  if (!all(is.finite(map_shared$BP))) {
    stop("无法完整获得rs33824/rs430752的GRCh38坐标。")
  }

  # 用两个共享SNP的中点定义共定位窗口中心
  signal_center <- round(median(map_shared$BP))

  windows <- data.table(
    window_id=c("250kb","500kb","1000kb"),
    flank=c(250000L,500000L,1000000L)
  )
  windows[, START := pmax(1L, signal_center-flank)]
  windows[, STOP := signal_center+flank]
  windows[, CHR := 19L]

  fwrite(
    windows,
    file.path(tables, "Step27b_coloc_windows_GRCh38.tsv"),
    sep="\t"
  )

  # 主分析窗口：±500kb
  main_window <- windows[window_id=="500kb"]
  fetch_window <- windows[window_id=="1000kb"]

  # ------------------------------------------------------------
  # 6. eQTL Catalogue官方metadata
  #    v4修复：
  #    1) 不再让fread直接访问GitHub URL，避免R默认60秒timeout；
  #    2) 使用curl/wget下载到本地，设置重试与更长timeout；
  #    3) 同时读取 uniformly processed 与 imported 两张tabix metadata，
  #       因为GTEx等数据可能位于 imported_tabix_paths.tsv。
  # ------------------------------------------------------------
  metadata_urls <- c(
    processed=paste0(
      "https://raw.githubusercontent.com/",
      "eQTL-Catalogue/eQTL-Catalogue-resources/master/",
      "tabix/tabix_ftp_paths.tsv"
    ),
    imported=paste0(
      "https://raw.githubusercontent.com/",
      "eQTL-Catalogue/eQTL-Catalogue-resources/master/",
      "tabix/tabix_ftp_paths_imported.tsv"
    )
  )

  curl_bin <- Sys.which("curl")
  wget_bin <- Sys.which("wget")

  download_metadata <- function(url, dest, label) {

    ok <- FALSE
    msg <- ""

    if (nzchar(curl_bin)) {

      rc <- system2(
        curl_bin,
        c(
          "-L",
          "--fail",
          "--retry","5",
          "--retry-delay","3",
          "--connect-timeout","30",
          "--max-time","300",
          "-A","Mozilla/5.0",
          "-o", shQuote(dest),
          shQuote(url)
        ),
        stdout=file.path(logs,paste0("metadata_",label,"_curl.log")),
        stderr=file.path(logs,paste0("metadata_",label,"_curl.log"))
      )

      ok <- identical(as.integer(rc),0L) &&
        file.exists(dest) &&
        file.info(dest)$size > 100

      if (!ok) msg <- paste0("curl rc=",rc)
    }

    if (!ok && nzchar(wget_bin)) {

      rc <- system2(
        wget_bin,
        c(
          "--tries=5",
          "--timeout=60",
          "--read-timeout=120",
          "-O", shQuote(dest),
          shQuote(url)
        ),
        stdout=file.path(logs,paste0("metadata_",label,"_wget.log")),
        stderr=file.path(logs,paste0("metadata_",label,"_wget.log"))
      )

      ok <- identical(as.integer(rc),0L) &&
        file.exists(dest) &&
        file.info(dest)$size > 100

      if (!ok) msg <- paste0(msg,"; wget rc=",rc)
    }

    list(ok=ok,message=msg,path=dest)
  }

  meta_parts <- list()
  meta_download_qc <- list()

  for (nm in names(metadata_urls)) {

    dest <- file.path(
      rawdir,
      paste0("eQTL_Catalogue_",nm,"_tabix_metadata.tsv")
    )

    cat("下载eQTL Catalogue metadata：",nm,"\n")

    dl <- download_metadata(
      metadata_urls[[nm]],
      dest,
      nm
    )

    if (dl$ok) {

      mm <- tryCatch(
        fread(dest),
        error=function(e) e
      )

      if (!inherits(mm,"error") && nrow(mm)) {

        mm[, metadata_source := nm]
        meta_parts[[nm]] <- mm

        meta_download_qc[[nm]] <- data.table(
          source=nm,
          status="OK",
          url=metadata_urls[[nm]],
          local_file=dest,
          file_size=file.info(dest)$size,
          n_rows=nrow(mm),
          message=""
        )

      } else {

        meta_download_qc[[nm]] <- data.table(
          source=nm,
          status="PARSE_FAILED",
          url=metadata_urls[[nm]],
          local_file=dest,
          file_size=if (file.exists(dest)) file.info(dest)$size else NA_real_,
          n_rows=0L,
          message=if (inherits(mm,"error")) conditionMessage(mm) else "zero rows"
        )
      }

    } else {

      meta_download_qc[[nm]] <- data.table(
        source=nm,
        status="DOWNLOAD_FAILED",
        url=metadata_urls[[nm]],
        local_file=dest,
        file_size=if (file.exists(dest)) file.info(dest)$size else NA_real_,
        n_rows=0L,
        message=dl$message
      )
    }
  }

  meta_qc <- rbindlist(meta_download_qc, fill=TRUE)

  fwrite(
    meta_qc,
    file.path(reports,"Step27b_metadata_download_QC.tsv"),
    sep="\t",
    na="NA"
  )

  if (!length(meta_parts)) {
    stop(
      "无法下载或解析eQTL Catalogue tabix metadata。\n",
      "请查看：",
      file.path(reports,"Step27b_metadata_download_QC.tsv"),
      "\n以及logs目录中的metadata下载日志。"
    )
  }

  meta <- rbindlist(
    meta_parts,
    fill=TRUE,
    use.names=TRUE
  )

  # 不同metadata版本中ID字段原则上应为dataset_id；若命名发生改变，明确报错。
  if (!"dataset_id" %in% names(meta)) {
    stop(
      "eQTL Catalogue metadata中未找到dataset_id字段。\n实际字段：",
      paste(names(meta),collapse=",")
    )
  }

  if (!"quant_method" %in% names(meta)) {
    stop(
      "eQTL Catalogue metadata中未找到quant_method字段。\n实际字段：",
      paste(names(meta),collapse=",")
    )
  }

  target_ids <- c(
    "QTD000296",
    "QTD000554",
    "QTD000574"
  )

  md <- meta[
    dataset_id %in% target_ids &
    quant_method=="ge"
  ]

  # 若同一dataset在两张metadata中同时出现，按dataset_id+ftp_path去重；
  # 不凭空选择不同路径。
  if ("ftp_path" %in% names(md)) {
    md <- unique(md, by=c("dataset_id","ftp_path"))
  } else {
    md <- unique(md)
  }

  if (!nrow(md)) {
    fwrite(
      meta[
        dataset_id %in% target_ids
      ],
      file.path(
        reports,
        "Step27b_target_dataset_metadata_any_quant_method.tsv"
      ),
      sep="\t",
      na="NA"
    )

    stop(
      "当前官方metadata中没有找到目标数据集的ge记录。\n",
      "已输出目标dataset的其他quant_method记录供复核。"
    )
  }

  fwrite(
    md,
    file.path(tables, "Step27b_eQTL_dataset_metadata.tsv"),
    sep="\t",
    na="NA"
  )

  # 对metadata缺失的数据集记录，而不是静默丢失
  missing_datasets <- setdiff(target_ids, unique(md$dataset_id))

  # ------------------------------------------------------------
  # 7. tabix远程读取函数
  # ------------------------------------------------------------
  make_https <- function(x) {
    sub("^ftp://ftp\\.ebi\\.ac\\.uk", "https://ftp.ebi.ac.uk", x)
  }

  # eQTL Catalogue的all.tsv.gz第一行是普通TSV列名，
  # 并不一定以#开头，因此tabix -H可能正常返回空结果。
  # 官方教程本身也是通过“curl | zcat | head -n 1”读取列名。
  # 这里先尝试tabix -H；若为空，则流式读取压缩文件第一行。
  stream_first_line_header <- function(url) {

    curl_bin2 <- Sys.which("curl")
    wget_bin2 <- Sys.which("wget")

    cmd <- NULL

    if (nzchar(curl_bin2)) {

      # head -n 1读取首行后会主动关闭pipe，curl可能产生23号broken-pipe提示；
      # 这是预期现象，不代表远程文件失败，因此静默curl stderr。
      cmd <- paste(
        shQuote(curl_bin2),
        "-L -s --fail --retry 3 --connect-timeout 30 --max-time 180",
        shQuote(url),
        "2>/dev/null | gzip -dc 2>/dev/null | head -n 1"
      )

    } else if (nzchar(wget_bin2)) {

      cmd <- paste(
        shQuote(wget_bin2),
        "-q -O -",
        shQuote(url),
        "| gzip -dc 2>/dev/null | head -n 1"
      )
    }

    if (is.null(cmd)) return(NULL)

    con <- pipe(cmd, open="r")
    on.exit(try(close(con), silent=TRUE), add=TRUE)

    h <- tryCatch(
      readLines(con, n=1, warn=FALSE),
      error=function(e) character()
    )

    if (!length(h) || !nzchar(h[1])) return(NULL)

    h <- sub("^#+","",h[1])

    if (!grepl("\t",h,fixed=TRUE)) return(NULL)

    strsplit(h,"\t",fixed=TRUE)[[1]]
  }

  tabix_header <- function(url) {

    tf <- tempfile(fileext=".header.txt")
    ef <- tempfile(fileext=".header.err")

    rc <- system2(
      tabix_bin,
      c("-H", shQuote(url)),
      stdout=tf,
      stderr=ef
    )

    header_names <- NULL

    if (identical(as.integer(rc),0L) && file.exists(tf)) {

      x <- readLines(tf, warn=FALSE)
      x <- x[nzchar(x)]

      cand <- x[grepl("\t",x,fixed=TRUE)]

      if (length(cand)) {
        h <- tail(cand,1)
        h <- sub("^#+","",h)
        header_names <- strsplit(h,"\t",fixed=TRUE)[[1]]
      }
    }

    # tabix -H空并不代表文件失败：eQTL Catalogue列名通常是普通首行。
    if (is.null(header_names) || !length(header_names)) {

      header_names <- stream_first_line_header(url)
    }

    # 保存诊断信息
    diag_file <- file.path(
      logs,
      "tabix_header_diagnostics.log"
    )

    cat(
      paste0(
        "URL=",url,
        "\ttabix_rc=",rc,
        "\tn_header_cols=",
        if (is.null(header_names)) 0L else length(header_names),
        "\n"
      ),
      file=diag_file,
      append=TRUE
    )

    if (file.exists(ef) && file.info(ef)$size>0) {
      errtxt <- readLines(ef,warn=FALSE)
      if (length(errtxt)) {
        cat(
          paste0("stderr=",paste(errtxt,collapse=" | "),"\n"),
          file=diag_file,
          append=TRUE
        )
      }
    }

    header_names
  }

  tabix_region <- function(url, region_string, header_names) {

    tf <- tempfile(fileext=".tsv")
    ef <- tempfile(fileext=".err")

    rc <- system2(
      tabix_bin,
      c(shQuote(url), shQuote(region_string)),
      stdout=tf,
      stderr=ef
    )

    if (!identical(as.integer(rc),0L)) {
      return(list(
        ok=FALSE,
        data=NULL,
        error=paste(readLines(ef,warn=FALSE),collapse=" | ")
      ))
    }

    if (!file.exists(tf) || file.info(tf)$size==0) {
      return(list(
        ok=TRUE,
        data=data.table(),
        error=""
      ))
    }

    x <- tryCatch(
      fread(tf, header=FALSE, showProgress=FALSE),
      error=function(e) e
    )

    if (inherits(x,"error")) {
      return(list(
        ok=FALSE,
        data=NULL,
        error=conditionMessage(x)
      ))
    }

    if (ncol(x)!=length(header_names)) {
      return(list(
        ok=FALSE,
        data=NULL,
        error=paste0(
          "column_count_mismatch: data=",
          ncol(x),
          "; header=",
          length(header_names)
        )
      ))
    }

    setnames(x, header_names)

    list(ok=TRUE,data=x,error="")
  }

  # ------------------------------------------------------------
  # 8. 提取PEPD eQTL
  # ------------------------------------------------------------
  pepd_gene_id <- sub("\\..*$","",pepd_info$gene_id)

  eqtl_list <- list()
  extract_qc <- list()

  fetch_region_string <- paste0(
    "19:",
    fetch_window$START,
    "-",
    fetch_window$STOP
  )

  for (i in seq_len(nrow(md))) {

    row <- md[i]
    ds <- row$dataset_id

    cat(
      "\n提取eQTL：", ds,
      " / ", row$study_label,
      " / ", row$tissue_label,
      "\n",
      sep=""
    )

    urls <- unique(c(
      as.character(row$ftp_path),
      make_https(as.character(row$ftp_path))
    ))

    success <- FALSE
    last_error <- ""
    qdat <- NULL
    used_url <- NA_character_

    for (u in urls) {

      hdr <- tabix_header(u)

      if (is.null(hdr)) {
        last_error <- paste0("无法读取header: ",u)
        next
      }

      rr <- tabix_region(
        u,
        fetch_region_string,
        hdr
      )

      if (!rr$ok) {
        last_error <- rr$error
        next
      }

      qdat <- rr$data
      used_url <- u
      success <- TRUE
      break
    }

    if (!success) {

      extract_qc[[ds]] <- data.table(
        dataset_id=ds,
        status="TABIX_FAILED",
        study_label=row$study_label,
        tissue_label=row$tissue_label,
        sample_size=row$sample_size,
        n_region_rows=0L,
        n_PEPD_rows=0L,
        used_url=NA_character_,
        message=last_error
      )

      next
    }

    if (!nrow(qdat)) {

      extract_qc[[ds]] <- data.table(
        dataset_id=ds,
        status="NO_REGION_ROWS",
        study_label=row$study_label,
        tissue_label=row$tissue_label,
        sample_size=row$sample_size,
        n_region_rows=0L,
        n_PEPD_rows=0L,
        used_url=used_url,
        message=""
      )

      next
    }

    # 字段要求
    reqq <- c(
      "chromosome","position","ref","alt",
      "rsid","maf","pvalue","beta","se","gene_id"
    )

    miss <- setdiff(reqq,names(qdat))

    if (length(miss)) {

      extract_qc[[ds]] <- data.table(
        dataset_id=ds,
        status="MISSING_COLUMNS",
        study_label=row$study_label,
        tissue_label=row$tissue_label,
        sample_size=row$sample_size,
        n_region_rows=nrow(qdat),
        n_PEPD_rows=0L,
        used_url=used_url,
        message=paste(miss,collapse=",")
      )

      next
    }

    qdat[, gene_id_clean := sub("\\..*$","",gene_id)]

    qpep <- qdat[
      gene_id_clean==pepd_gene_id
    ]

    if (!nrow(qpep)) {

      extract_qc[[ds]] <- data.table(
        dataset_id=ds,
        status="PEPD_NOT_PRESENT_IN_REGION",
        study_label=row$study_label,
        tissue_label=row$tissue_label,
        sample_size=row$sample_size,
        n_region_rows=nrow(qdat),
        n_PEPD_rows=0L,
        used_url=used_url,
        message=""
      )

      next
    }

    # 保留有效SNV/eQTL记录
    qpep <- qpep[
      nchar(ref)==1 &
      nchar(alt)==1 &
      ref %chin% c("A","C","G","T") &
      alt %chin% c("A","C","G","T") &
      is.finite(beta) &
      is.finite(se) &
      se>0 &
      is.finite(maf) &
      maf>0 &
      maf<=0.5
    ]

    # eQTL Catalogue同一variant可因multiple rsids重复。
    # 先构建variant key，再对完全重复variant保留一条；
    # 后续与GWAS merge时要求rsID唯一。
    qpep[, variant_key := paste(
      chromosome,position,ref,alt,sep=":"
    )]

    # 对完全相同variant+rsid去重
    qpep <- unique(
      qpep,
      by=c("variant_key","rsid")
    )

    qpep[, dataset_id := ds]
    qpep[, study_label := row$study_label]
    qpep[, tissue_label := row$tissue_label]
    qpep[, sample_size_metadata := as.numeric(row$sample_size)]

    # raw an/2作为样本量交叉检查
    if ("an" %in% names(qpep)) {
      qpep[, N_from_AN := as.numeric(an)/2]
    } else {
      qpep[, N_from_AN := NA_real_]
    }

    eqtl_list[[ds]] <- qpep

    fwrite(
      qpep,
      file.path(
        rawdir,
        paste0(ds,"_PEPD_chr19_GRCh38.tsv.gz")
      ),
      sep="\t",
      na="NA"
    )

    extract_qc[[ds]] <- data.table(
      dataset_id=ds,
      status="OK",
      study_label=row$study_label,
      tissue_label=row$tissue_label,
      sample_size=row$sample_size,
      n_region_rows=nrow(qdat),
      n_PEPD_rows=nrow(qpep),
      used_url=used_url,
      message=""
    )
  }

  if (length(missing_datasets)) {
    for (ds in missing_datasets) {
      extract_qc[[paste0(ds,"_missing")]] <- data.table(
        dataset_id=ds,
        status="NOT_FOUND_IN_CURRENT_METADATA",
        study_label=NA_character_,
        tissue_label=NA_character_,
        sample_size=NA_real_,
        n_region_rows=0L,
        n_PEPD_rows=0L,
        used_url=NA_character_,
        message="Target dataset ID not present among current ge metadata rows."
      )
    }
  }

  qctl <- rbindlist(extract_qc, fill=TRUE)

  fwrite(
    qctl,
    file.path(reports, "Step27b_eQTL_extraction_QC.tsv"),
    sep="\t",
    na="NA"
  )

  # 至少一个PEPD eQTL数据集必须成功
  ok_ds <- qctl[status=="OK",dataset_id]

  if (!length(ok_ds)) {
    writeLines(
      "COMPLETED_STEP27B_WITH_REVIEW_FLAGS",
      file.path(reports,"completion_status.txt")
    )
    stop("三个目标数据集中没有任何一个成功提取到PEPD eQTL。")
  }

  # ------------------------------------------------------------
  # 9. rs33824 / rs430752 eQTL覆盖
  # ------------------------------------------------------------
  qtl_priority_rows <- list()

  for (ds in names(eqtl_list)) {

    q <- eqtl_list[[ds]]

    qq <- q[
      rsid %in% expected_chr19
    ]

    if (nrow(qq)) {
      qtl_priority_rows[[ds]] <- qq
    }
  }

  priority_qtl <- if (length(qtl_priority_rows))
    rbindlist(qtl_priority_rows, fill=TRUE) else data.table()

  if (nrow(priority_qtl)) {
    fwrite(
      priority_qtl,
      file.path(
        tables,
        "Step27b_rs33824_rs430752_PEPD_eQTL_stats.tsv"
      ),
      sep="\t",
      na="NA"
    )
  } else {
    fwrite(
      data.table(
        note="Neither rs33824 nor rs430752 was directly present among extracted PEPD eQTL rows."
      ),
      file.path(
        tables,
        "Step27b_rs33824_rs430752_PEPD_eQTL_stats.tsv"
      ),
      sep="\t"
    )
  }

  coverage <- CJ(
    dataset_id=ok_ds,
    SNP=expected_chr19,
    unique=TRUE
  )

  coverage[, present := FALSE]

  for (i in seq_len(nrow(coverage))) {
    ds <- coverage$dataset_id[i]
    coverage$present[i] <- coverage$SNP[i] %in% eqtl_list[[ds]]$rsid
  }

  fwrite(
    coverage,
    file.path(
      tables,
      "Step27b_priority_SNP_eQTL_coverage.tsv"
    ),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 10. 读取两个GWAS区域summary statistics
  # ------------------------------------------------------------
  gwas_files <- c(
    pancreas_fat=file.path(
      step26,"regional_sumstats",
      "pancreas_fat__chr19_LOC2340_primary.tsv.gz"
    ),
    pancreas_iron=file.path(
      step26,"regional_sumstats",
      "pancreas_iron__chr19_LOC2340_primary.tsv.gz"
    )
  )

  gwas <- lapply(gwas_files, fread)

  gN <- c(
    pancreas_fat=25617,
    pancreas_iron=25617
  )

  for (tr in names(gwas)) {

    x <- gwas[[tr]]

    req <- c(
      "SNP","REF","ALT","BETA","SE","VARBETA","AF","MAF"
    )

    if (!all(req %in% names(x))) {
      stop(
        "GWAS区域文件字段异常：", tr,
        "\n缺少：",
        paste(setdiff(req,names(x)),collapse=",")
      )
    }

    x <- x[
      grepl("^rs[0-9]+$",SNP) &
      is.finite(BETA) &
      is.finite(SE) &
      SE>0 &
      is.finite(MAF) &
      MAF>0 &
      MAF<=0.5
    ]

    # rsID必须唯一，重复rsID全部剔除
    dup <- x[, .N, by=SNP][N>1,SNP]
    if (length(dup)) x <- x[!SNP %in% dup]

    gwas[[tr]] <- x
  }

  # ------------------------------------------------------------
  # 11. harmonization：GWAS effect allele=ALT；eQTL effect allele=ALT
  # ------------------------------------------------------------
  comp <- function(a) chartr("ACGT","TGCA",toupper(a))

  harmonize_pair <- function(g, q, dataset_id, trait) {

    # eQTL rsID唯一化：同一个rsid若对应多个variant记录，全部剔除
    q2 <- copy(q)

    dupq <- q2[, .N, by=rsid][N>1,rsid]
    if (length(dupq)) q2 <- q2[!rsid %in% dupq]

    m <- merge(
      g,
      q2,
      by.x="SNP",
      by.y="rsid",
      suffixes=c("_gwas","_eqtl"),
      all=FALSE
    )

    if (!nrow(m)) return(m)

    # GWAS REF/ALT与eQTL ref/alt
    same <- (
      m$REF == m$ref &
      m$ALT == m$alt
    )

    swap <- (
      m$REF == m$alt &
      m$ALT == m$ref
    )

    comp_same <- (
      comp(m$REF) == m$ref &
      comp(m$ALT) == m$alt
    )

    comp_swap <- (
      comp(m$REF) == m$alt &
      comp(m$ALT) == m$ref
    )

    ok <- same | swap | comp_same | comp_swap
    m <- m[ok]

    same <- same[ok]
    swap <- swap[ok]
    comp_same <- comp_same[ok]
    comp_swap <- comp_swap[ok]

    if (!nrow(m)) return(m)

    flip_eqtl <- swap | comp_swap

    # 将eQTL beta对齐到GWAS ALT方向，便于方向性描述
    m[, beta_eqtl_aligned :=
      ifelse(flip_eqtl, -beta, beta)
    ]

    # eQTL maf本身不依赖effect allele方向
    m[, harmonization :=
      fifelse(
        same,"same",
        fifelse(
          swap,"swap",
          fifelse(
            comp_same,"complement_same",
            "complement_swap"
          )
        )
      )
    ]

    m[, dataset_id := dataset_id]
    m[, gwas_trait := trait]

    # 主分析要求双方MAF>=0.01
    m[, common_MAF_ge_001 :=
      MAF>=0.01 &
      maf>=0.01
    ]

    setorder(m, position)

    m
  }

  harm_list <- list()
  harm_qc <- list()

  for (ds in names(eqtl_list)) {

    q <- eqtl_list[[ds]]

    for (tr in names(gwas)) {

      h <- harmonize_pair(
        gwas[[tr]],
        q,
        ds,
        tr
      )

      key <- paste(ds,tr,sep="__")
      harm_list[[key]] <- h

      harm_qc[[key]] <- data.table(
        dataset_id=ds,
        gwas_trait=tr,
        n_gwas_rsids=nrow(gwas[[tr]]),
        n_eqtl_rows=nrow(q),
        n_common_harmonized=nrow(h),
        n_common_MAF_ge_001=if (nrow(h))
          sum(h$common_MAF_ge_001) else 0L,
        rs33824_common=if (nrow(h))
          "rs33824" %in% h$SNP else FALSE,
        rs430752_common=if (nrow(h))
          "rs430752" %in% h$SNP else FALSE
      )

      if (nrow(h)) {
        fwrite(
          h,
          file.path(
            tables,
            paste0(
              "harmonized_",
              ds,"_",tr,".tsv.gz"
            )
          ),
          sep="\t",
          na="NA"
        )
      }
    }
  }

  harm_qc_tab <- rbindlist(harm_qc, fill=TRUE)

  fwrite(
    harm_qc_tab,
    file.path(reports, "Step27b_harmonization_QC.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 12. coloc.abf
  # ------------------------------------------------------------
  run_coloc <- function(h, dataset_row, trait, window_id, p12) {

    if (is.null(h) || !nrow(h)) {
      return(list(status="NO_COMMON_VARIANTS"))
    }

    win_id <- as.character(window_id)
    ww <- windows[window_id == win_id]

    x <- h[
      position>=ww$START &
      position<=ww$STOP &
      common_MAF_ge_001==TRUE
    ]

    # coloc至少需要足够区域覆盖
    if (nrow(x)<100) {
      return(list(
        status="TOO_FEW_COMMON_VARIANTS",
        n=nrow(x)
      ))
    }

    # qtl N：优先metadata sample_size；同时保留AN/2 QC
    qN <- suppressWarnings(
      as.numeric(dataset_row$sample_size)
    )

    if (!is.finite(qN) || qN<=0) {

      if ("N_from_AN" %in% names(x)) {
        qN <- median(
          x$N_from_AN[
            is.finite(x$N_from_AN) &
            x$N_from_AN>0
          ],
          na.rm=TRUE
        )
      }
    }

    if (!is.finite(qN) || qN<=0) {
      return(list(
        status="EQTL_SAMPLE_SIZE_UNAVAILABLE",
        n=nrow(x)
      ))
    }

    # eQTL Catalogue官方tutorial采用beta/se/maf/N/type=quant，
    # 不强制指定sdY，让coloc根据summary statistics估计。
    D_eqtl <- list(
      beta=x$beta,
      varbeta=x$se^2,
      MAF=x$maf,
      N=qN,
      type="quant",
      snp=x$SNP
    )

    D_gwas <- list(
      beta=x$BETA,
      varbeta=x$VARBETA,
      MAF=x$MAF,
      N=as.numeric(gN[trait]),
      type="quant",
      snp=x$SNP
    )

    res <- tryCatch(
      coloc::coloc.abf(
        dataset1=D_eqtl,
        dataset2=D_gwas,
        p1=1e-4,
        p2=1e-4,
        p12=p12
      ),
      error=function(e) e
    )

    if (inherits(res,"error")) {
      return(list(
        status=paste0(
          "COLOC_ERROR:",
          conditionMessage(res)
        ),
        n=nrow(x),
        qN=qN
      ))
    }

    list(
      status="OK",
      n=nrow(x),
      qN=qN,
      summary=res$summary,
      results=res$results
    )
  }

  grid <- CJ(
    dataset_id=names(eqtl_list),
    gwas_trait=names(gwas),
    window_id=windows$window_id,
    p12=c(1e-4,1e-5,1e-6),
    unique=TRUE
  )

  coloc_rows <- list()
  coloc_objs <- list()

  for (i in seq_len(nrow(grid))) {

    g <- grid[i]
    key <- paste(g$dataset_id,g$gwas_trait,sep="__")

    mdr <- md[dataset_id==g$dataset_id][1]

    rr <- run_coloc(
      harm_list[[key]],
      mdr,
      g$gwas_trait,
      g$window_id,
      g$p12
    )

    objkey <- paste(
      g$dataset_id,
      g$gwas_trait,
      g$window_id,
      format(g$p12,scientific=TRUE),
      sep="__"
    )
    coloc_objs[[objkey]] <- rr

    row <- data.table(
      dataset_id=g$dataset_id,
      study_label=mdr$study_label,
      tissue_label=mdr$tissue_label,
      sample_size=mdr$sample_size,
      gwas_trait=g$gwas_trait,
      window_id=g$window_id,
      p1=1e-4,
      p2=1e-4,
      p12=g$p12,
      status=rr$status,
      n_common_variants=if (!is.null(rr$n)) rr$n else NA_integer_,
      eqtl_N_used=if (!is.null(rr$qN)) rr$qN else NA_real_,
      PP.H0=NA_real_,
      PP.H1=NA_real_,
      PP.H2=NA_real_,
      PP.H3=NA_real_,
      PP.H4=NA_real_
    )

    if (identical(rr$status,"OK")) {
      s <- rr$summary
      row[, `:=`(
        PP.H0=as.numeric(s["PP.H0.abf"]),
        PP.H1=as.numeric(s["PP.H1.abf"]),
        PP.H2=as.numeric(s["PP.H2.abf"]),
        PP.H3=as.numeric(s["PP.H3.abf"]),
        PP.H4=as.numeric(s["PP.H4.abf"])
      )]
    }

    coloc_rows[[i]] <- row
  }

  coloc_all <- rbindlist(coloc_rows, fill=TRUE)

  fwrite(
    coloc_all,
    file.path(
      colocdir,
      "Step27b_PEPD_GWAS_eQTL_coloc_all_sensitivity.tsv"
    ),
    sep="\t",
    na="NA"
  )

  # ------------------------------------------------------------
  # 13. 主分析：±500kb, p12=1e-5
  # ------------------------------------------------------------
  main <- coloc_all[
    window_id=="500kb" &
    abs(p12-1e-5)<1e-12
  ]

  fwrite(
    main,
    file.path(
      colocdir,
      "Step27b_PEPD_GWAS_eQTL_coloc_main.tsv"
    ),
    sep="\t",
    na="NA"
  )

  # 每个dataset×trait的窗口/先验稳健性
  robust <- coloc_all[
    status=="OK",
    .(
      n_runs=.N,
      min_PP_H4=min(PP.H4,na.rm=TRUE),
      median_PP_H4=median(PP.H4,na.rm=TRUE),
      max_PP_H4=max(PP.H4,na.rm=TRUE),
      n_PP_H4_ge_080=sum(PP.H4>=0.80,na.rm=TRUE),
      n_PP_H4_ge_090=sum(PP.H4>=0.90,na.rm=TRUE),
      max_PP_H3=max(PP.H3,na.rm=TRUE)
    ),
    by=.(dataset_id,study_label,tissue_label,gwas_trait)
  ]

  fwrite(
    robust,
    file.path(
      colocdir,
      "Step27b_PEPD_coloc_robustness_summary.tsv"
    ),
    sep="\t",
    na="NA"
  )

  # ------------------------------------------------------------
  # 14. 主分析SNP级H4 posterior
  # ------------------------------------------------------------
  snp_post_rows <- list()

  for (i in seq_len(nrow(main))) {

    if (main$status[i]!="OK") next

    ds <- main$dataset_id[i]
    tr <- main$gwas_trait[i]

    # 重新用固定主设定，避免依赖objkey字符串格式
    mdr <- md[dataset_id==ds][1]
    key <- paste(ds,tr,sep="__")

    rr <- run_coloc(
      harm_list[[key]],
      mdr,
      tr,
      "500kb",
      1e-5
    )

    if (!identical(rr$status,"OK")) next

    sr <- as.data.table(rr$results)

    if (!"SNP.PP.H4" %in% names(sr)) next

    sr[, dataset_id := ds]
    sr[, study_label := mdr$study_label]
    sr[, tissue_label := mdr$tissue_label]
    sr[, gwas_trait := tr]

    setorder(sr,-SNP.PP.H4)

    snp_post_rows[[paste(ds,tr,sep="__")]] <- sr
  }

  snp_post <- if (length(snp_post_rows))
    rbindlist(snp_post_rows, fill=TRUE) else data.table()

  if (nrow(snp_post)) {

    fwrite(
      snp_post,
      file.path(
        colocdir,
        "Step27b_main_coloc_SNP_H4_posterior.tsv.gz"
      ),
      sep="\t",
      na="NA"
    )

    top20 <- snp_post[
      order(-SNP.PP.H4),
      head(.SD,20),
      by=.(dataset_id,study_label,tissue_label,gwas_trait)
    ]

    fwrite(
      top20,
      file.path(
        colocdir,
        "Step27b_main_coloc_top20_H4_variants.tsv"
      ),
      sep="\t",
      na="NA"
    )
  }

  # ------------------------------------------------------------
  # 15. rs33824 / rs430752在harmonized数据中的方向性表
  # ------------------------------------------------------------
  priority_harm <- list()

  for (key in names(harm_list)) {
    h <- harm_list[[key]]
    if (!is.null(h) && nrow(h)) {
      z <- h[SNP %in% expected_chr19]
      if (nrow(z)) priority_harm[[key]] <- z
    }
  }

  priority_harm_tab <- if (length(priority_harm))
    rbindlist(priority_harm, fill=TRUE) else data.table()

  if (nrow(priority_harm_tab)) {

    keep <- intersect(
      c(
        "dataset_id","study_label","tissue_label",
        "gwas_trait","SNP",
        "BETA","SE","MAF",
        "beta","se","maf",
        "beta_eqtl_aligned",
        "pvalue","harmonization"
      ),
      names(priority_harm_tab)
    )

    fwrite(
      priority_harm_tab[, ..keep],
      file.path(
        tables,
        "Step27b_rs33824_rs430752_GWAS_eQTL_direction.tsv"
      ),
      sep="\t",
      na="NA"
    )
  }

  # ------------------------------------------------------------
  # 16. 自动证据分类（仅用于结果整理，不等于因果结论）
  # ------------------------------------------------------------
  classify <- function(h4,h3) {
    if (!is.finite(h4)) return("NOT_ESTIMABLE")
    if (h4>=0.90) return("STRONG_SHARED_SIGNAL_SUPPORT")
    if (h4>=0.80) return("MODERATE_TO_STRONG_SHARED_SIGNAL_SUPPORT")
    if (h4>=0.50) return("INTERMEDIATE_SHARED_SIGNAL_SUPPORT")
    if (is.finite(h3) && h3>=0.80) return("DISTINCT_SIGNAL_MODEL_DOMINANT")
    return("LIMITED_SHARED_SIGNAL_SUPPORT")
  }

  main[, evidence_label :=
    mapply(classify,PP.H4,PP.H3)
  ]

  fwrite(
    main,
    file.path(
      reports,
      "Step27b_main_evidence_labels.tsv"
    ),
    sep="\t",
    na="NA"
  )

  # ------------------------------------------------------------
  # 17. QC与总结
  # ------------------------------------------------------------
  summary <- data.table(
    metric=c(
      "n_target_datasets",
      "n_datasets_found_in_metadata",
      "n_datasets_with_PEPD_eQTL",
      "n_priority_SNP_dataset_combinations_covered",
      "n_harmonized_dataset_trait_pairs",
      "n_main_coloc_pairs",
      "n_main_coloc_OK",
      "n_main_PP_H4_ge_080",
      "n_main_PP_H4_ge_090"
    ),
    value=c(
      length(target_ids),
      uniqueN(md$dataset_id),
      length(ok_ds),
      coverage[present==TRUE,.N],
      harm_qc_tab[n_common_harmonized>0,.N],
      nrow(main),
      main[status=="OK",.N],
      main[status=="OK" & PP.H4>=0.80,.N],
      main[status=="OK" & PP.H4>=0.90,.N]
    )
  )

  fwrite(
    summary,
    file.path(reports, "Step27b_summary.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 18. 解释说明
  # ------------------------------------------------------------
  writeLines(
    c(
      "Step27b evaluates PEPD expression QTL evidence in pancreatic tissues/islets.",
      paste0(
        "Target eQTL Catalogue datasets: ",
        paste(target_ids,collapse=", "),
        "."
      ),
      "Dataset labels, tissue labels, sample sizes and FTP paths are resolved dynamically from the official eQTL Catalogue tabix metadata.",
      paste0(
        "PEPD Ensembl gene ID resolved from Ensembl GRCh38: ",
        pepd_info$gene_id,
        "."
      ),
      "MRI GWAS are GRCh37 and eQTL Catalogue data are GRCh38; cross-build harmonisation is performed using unique rsIDs plus allele compatibility rather than direct coordinate matching.",
      "The eQTL Catalogue ALT allele is treated as the eQTL effect allele.",
      "Primary colocalisation uses a +/-500 kb window around the GRCh38 midpoint of rs33824 and rs430752, MAF>=0.01 in both datasets, p1=p2=1e-4 and p12=1e-5.",
      "Sensitivity analyses use +/-250 kb, +/-500 kb and +/-1 Mb windows and p12=1e-4, 1e-5 and 1e-6.",
      "eQTL sample size is taken from official metadata, with AN/2 retained as a QC field where available.",
      "coloc PP.H4 represents support for a shared causal-signal model under the specified priors and assumptions; it does not establish that PEPD expression mediates pancreatic fat or pancreatic iron.",
      "A weak or absent PEPD eQTL colocalisation result does not exclude PEPD function, particularly when tissue sample size or eQTL power is limited.",
      "chr19 was selected in Step25b sensitivity analysis, so Step27b is mechanistic follow-up rather than independent replication."
    ),
    file.path(reports, "methods_and_interpretation.txt")
  )

  # ------------------------------------------------------------
  # 19. 完成判定
  # ------------------------------------------------------------
  review_flags <- character()

  if (length(ok_ds)<1)
    review_flags <- c(review_flags,"NO_PEPD_EQTL_DATASET")

  if (main[status=="OK",.N]<1)
    review_flags <- c(review_flags,"NO_SUCCESSFUL_MAIN_COLOC")

  if (length(missing_datasets))
    review_flags <- c(
      review_flags,
      paste0(
        "DATASET_NOT_IN_CURRENT_METADATA:",
        paste(missing_datasets,collapse=",")
      )
    )

  if (length(review_flags)) {

    writeLines(
      "COMPLETED_STEP27B_WITH_REVIEW_FLAGS",
      file.path(reports,"completion_status.txt")
    )

    writeLines(
      review_flags,
      file.path(reports,"review_flags.txt")
    )

  } else {

    writeLines(
      "COMPLETED_STEP27B_PEPD_EQTL_COLOCALIZATION",
      file.path(reports,"completion_status.txt")
    )
  }

  # ------------------------------------------------------------
  # 20. 精简review包
  # ------------------------------------------------------------
  zipfile <- file.path(
    out,
    "27b_PEPD_eQTL_coloc_review.zip"
  )

  oldwd <- getwd()
  on.exit(setwd(oldwd),add=TRUE)
  setwd(out)

  files_to_zip <- c(
    list.files("reports",recursive=TRUE,full.names=TRUE),
    list.files("tables",recursive=TRUE,full.names=TRUE),
    list.files("coloc",recursive=TRUE,full.names=TRUE)
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

  # ------------------------------------------------------------
  # 21. 控制台摘要
  # ------------------------------------------------------------
  cat("\n========================================\n")
  cat("第27b步完成\n")
  cat("========================================\n\n")

  cat("PEPD：",pepd_info$gene_id,"\n")
  cat(
    "rs33824 / rs430752 GRCh38中心：",
    signal_center,
    "\n"
  )

  cat("\neQTL提取QC：\n")
  print(qctl)

  cat("\n优先SNP覆盖：\n")
  print(coverage)

  cat("\n主coloc结果（±500kb, p12=1e-5）：\n")
  print(main)

  cat("\n稳健性摘要：\n")
  print(robust)

  cat("\nStep27b summary：\n")
  print(summary)

  cat("\n结果目录：",out,"\n")
  cat("复核包：",zipfile,"\n")

})()
