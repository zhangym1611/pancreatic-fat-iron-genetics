# 第26步：chr19候选区域的cross-trait colocalization
# 修正版v4：固定第19步Python环境；修复BIM字段冲突；修复SuSiE等位基因匹配中REF_fat/REF_A1未显式引用data.table列的问题。
#          + chr8既有强共享区域阳性参照
# 目标表型：pancreas fat vs pancreas iron
#
# 主问题：
# 25b在chr19 LOC2340 (GRCh37 chr19:33785836-34633274)发现探索性local rg。
# 本步检验该局部遗传相关是否与“两表型共享同一局部因果信号”的模型一致。
#
# 主分析：
# - chr19 LOC2340原始LAVA边界
# - 原始GWAS完整区域summary statistics（不做P值筛选）
# - 双等位SNV、PASS/.、有效ES/SE/AF
# - MAF >= 0.01
# - coloc.abf：p1=p2=1e-4，p12=1e-5
# - quantitative trait，sdY=1
#
# 预设敏感性：
# - chr19：原始LAVA边界、边界±250 kb、边界±500 kb
# - chr8阳性参照：GRCh37 chr8:129047378-130047378
# - MAF：全部有效变异 vs MAF>=0.01
# - sdY：1 vs coloc从summary statistics估计
# - p12：1e-4, 1e-5, 1e-6
#
# 补充多信号分析：
# - 若现有1000G EUR PLINK reference可用，则对chr19主区域及chr8参照区域
#   进行reference-matched SuSiE-RSS + coloc.susie。
# - SuSiE结果为补充证据；若LD覆盖或数值条件不足，不使ABF主分析失败。
#
# 解释边界：
# - PP.H4支持“区域内两表型共享因果变异/共享信号”的模型，不代表一个表型导致另一个。
# - local rg显著不等于colocalization；colocalization阴性也不否定广泛polygenic sharing。
# - chr19来自25b探索性筛选，因此chr19共定位结果属于候选区域后续验证，不能当作独立复制。

(function() {

  options(stringsAsFactors = FALSE)
  base <- path.expand("~/胰腺")

  # ------------------------------------------------------------
  # 1. 依赖包
  # ------------------------------------------------------------
  cran_lib <- path.expand("~/R_pancreas_step26")
  dir.create(cran_lib, recursive=TRUE, showWarnings=FALSE)
  .libPaths(c(cran_lib, .libPaths()))

  install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
      install.packages(pkg, lib=cran_lib, repos="https://cloud.r-project.org")
    }
    if (!requireNamespace(pkg, quietly=TRUE))
      stop("无法加载R包：", pkg)
  }

  install_if_missing("data.table")
  install_if_missing("coloc")
  install_if_missing("susieR")

  suppressPackageStartupMessages(library(data.table))

  # ------------------------------------------------------------
  # 2. 自动确认第25b候选
  # ------------------------------------------------------------
  step25b_parent <- file.path(base, "25b_LAVA_local_rg_sensitivity")
  if (!dir.exists(step25b_parent))
    stop("未找到第25b目录：", step25b_parent)

  dirs25b <- list.dirs(step25b_parent, recursive=FALSE, full.names=TRUE)
  completed25b <- dirs25b[vapply(dirs25b, function(d) {
    f <- file.path(d, "reports", "completion_status.txt")
    file.exists(f) &&
      identical(readLines(f, warn=FALSE)[1],
                "COMPLETED_STEP25B_LAVA_NOMINAL_H2_SENSITIVITY")
  }, logical(1))]

  if (!length(completed25b))
    stop("没有找到完成的第25b分析目录。")

  step25b <- completed25b[which.max(file.info(completed25b)$mtime)]
  f25b <- file.path(step25b, "reports", "Step25b_LAVA_local_rg_all_tested.tsv")
  if (!file.exists(f25b))
    stop("缺少第25b local rg结果：", f25b)

  lava25b <- fread(f25b)
  req25 <- c("LOC","CHR","START","STOP","p","FDR_BH_tested","Bonferroni_tested")
  if (!all(req25 %in% names(lava25b)))
    stop("第25b结果字段异常。")

  chr19_hit <- lava25b[
    LOC == 2340 & CHR == 19
  ]

  if (nrow(chr19_hit) != 1L)
    stop("无法唯一确认第25b LOC2340 chr19候选。")

  if (chr19_hit$START != 33785836L || chr19_hit$STOP != 34633274L)
    stop("LOC2340坐标与预期不一致，请先复核第25b。")

  # ------------------------------------------------------------
  # 3. 输出目录
  # ------------------------------------------------------------
  parent <- file.path(base, "26_cross_trait_colocalization")
  dir.create(parent, recursive=TRUE, showWarnings=FALSE)
  out <- tempfile("fat_iron_chr19_coloc_", tmpdir=parent)
  if (!dir.create(out)) stop("无法创建第26步输出目录：", out)

  reports <- file.path(out, "reports")
  extracted <- file.path(out, "regional_sumstats")
  abfdir <- file.path(out, "coloc_abf")
  susiedir <- file.path(out, "coloc_susie")
  logs <- file.path(out, "logs")
  for (d in c(reports, extracted, abfdir, susiedir, logs))
    dir.create(d)

  writeLines("INCOMPLETE_CHECK_LOGS",
             file.path(reports, "completion_status.txt"))

  # ------------------------------------------------------------
  # 4. 固定区域
  # ------------------------------------------------------------
  regions <- data.table(
    region_id = c(
      "chr19_LOC2340_primary",
      "chr19_LOC2340_flank250kb",
      "chr19_LOC2340_flank500kb",
      "chr8_positive_control"
    ),
    CHR = c(19L,19L,19L,8L),
    START = c(
      33785836L,
      max(1L,33785836L-250000L),
      max(1L,33785836L-500000L),
      129047378L
    ),
    STOP = c(
      34633274L,
      34633274L+250000L,
      34633274L+500000L,
      130047378L
    ),
    role = c(
      "PRIMARY_CANDIDATE",
      "BOUNDARY_SENSITIVITY",
      "BOUNDARY_SENSITIVITY",
      "POSITIVE_CONTROL"
    )
  )

  fwrite(regions, file.path(reports, "Step26_regions.tsv"), sep="\t")

  # ------------------------------------------------------------
  # 5. 原始GWAS VCF
  # ------------------------------------------------------------
  ids <- c(
    pancreas_fat  = "ebi-a-GCST90016675",
    pancreas_iron = "ebi-a-GCST90016676"
  )
  N <- c(
    pancreas_fat=25617L,
    pancreas_iron=25617L
  )

  vcfs <- file.path(base, paste0(ids, ".vcf.gz"))
  names(vcfs) <- names(ids)

  if (!all(file.exists(vcfs)))
    stop("缺少原始GWAS VCF：\n",
         paste(vcfs[!file.exists(vcfs)], collapse="\n"))

  # ------------------------------------------------------------
  # 6. Python一次扫描每份VCF并提取所有Step26区域
  # ------------------------------------------------------------
  pyfile <- file.path(out, "extract_step26_regions.py")

  py <- c(
'from pathlib import Path',
'import gzip, sys, math, traceback',
'import pandas as pd',
'',
'base=Path(sys.argv[1])',
'out=Path(sys.argv[2])',
'region_file=Path(sys.argv[3])',
'regions=pd.read_csv(region_file,sep="\\t")',
'traits={',
' "pancreas_fat":"ebi-a-GCST90016675",',
' "pancreas_iron":"ebi-a-GCST90016676"',
'}',
'',
'def num(x):',
'    try:',
'        v=float(x)',
'        return v if math.isfinite(v) else None',
'    except Exception:',
'        return None',
'',
'def parse_sample(fmt,val):',
'    k=fmt.split(":"); v=val.split(":")',
'    if len(k)!=len(v): return None',
'    return dict(zip(k,v))',
'',
'def scan(trait,sid):',
'    vcf=base/(sid+".vcf.gz")',
'    rows={r.region_id:[] for r in regions.itertuples(index=False)}',
'    counts={r.region_id:dict(raw_in_region=0,eligible=0) for r in regions.itertuples(index=False)}',
'    with gzip.open(vcf,"rt") as f:',
'        for line in f:',
'            if line.startswith("#"): continue',
'            p=line.rstrip("\\n").split("\\t")',
'            if len(p)<10: continue',
'            chrom=p[0].replace("chr","")',
'            if not chrom.isdigit(): continue',
'            ch=int(chrom); bp=int(p[1])',
'            hits=regions[(regions.CHR==ch)&(regions.START<=bp)&(regions.STOP>=bp)]',
'            if len(hits)==0: continue',
'            for rr in hits.itertuples(index=False): counts[rr.region_id]["raw_in_region"]+=1',
'            ref=p[3].upper(); alt=p[4].upper(); filt=p[6]',
'            if len(ref)!=1 or len(alt)!=1 or ref not in "ACGT" or alt not in "ACGT": continue',
'            if filt not in ("PASS","."): continue',
'            sm=parse_sample(p[8],p[9])',
'            if sm is None or not all(k in sm for k in ("ES","SE","AF")): continue',
'            beta=num(sm["ES"]); se=num(sm["SE"]); af=num(sm["AF"])',
'            if beta is None or se is None or af is None or se<=0 or not (0<af<1): continue',
'            rsid=p[2]',
'            maf=min(af,1-af)',
'            rec=dict(trait=trait,CHR=ch,BP=bp,SNP=rsid,REF=ref,ALT=alt,',
'                     BETA=beta,SE=se,VARBETA=se*se,AF=af,MAF=maf)',
'            for rr in hits.itertuples(index=False):',
'                rows[rr.region_id].append(dict(region_id=rr.region_id,**rec))',
'                counts[rr.region_id]["eligible"]+=1',
'    for rid,x in rows.items():',
'        d=pd.DataFrame(x)',
'        d.to_csv(out/(trait+"__"+rid+".tsv.gz"),sep="\\t",index=False,compression="gzip")',
'    q=[]',
'    for rid,c in counts.items(): q.append(dict(trait=trait,region_id=rid,**c))',
'    return q',
'',
'try:',
'    allq=[]',
'    for t,s in traits.items():',
'        allq.extend(scan(t,s))',
'    pd.DataFrame(allq).to_csv(out/"extraction_QC.tsv",sep="\\t",index=False)',
'except Exception:',
'    (out/"python_error.txt").write_text(traceback.format_exc())',
'    raise'
  )

  writeLines(py, pyfile, useBytes=TRUE)

  region_file <- file.path(out, "regions_for_python.tsv")
  fwrite(regions[,.(region_id,CHR,START,STOP)], region_file, sep="\t")

  # 固定使用第19步已验证的软件环境，避免系统python缺少pandas/numpy。
  step19 <- file.path(
    base,
    "19_LDSC_setup",
    "LDSC_setup_1566e430f1187"
  )

  cfg_file <- file.path(step19, "LDSC_paths.rds")

  if (!file.exists(cfg_file)) {
    stop("缺少第19步环境配置文件：", cfg_file)
  }

  cfg19 <- readRDS(cfg_file)

  if (is.null(cfg19$python) || !file.exists(cfg19$python)) {
    stop(
      "第19步Python环境不存在或路径无效：",
      cfg19$python
    )
  }

  python_bin <- cfg19$python

  cat(
    "第26步固定使用第19步Python：",
    python_bin,
    "\n"
  )

  pylog <- file.path(logs, "regional_extraction_console.txt")
  rc <- system2(
    python_bin,
    args=vapply(c(pyfile, base, extracted, region_file),
                shQuote, character(1)),
    stdout=pylog,
    stderr=pylog,
    env=c(
      "PYTHONNOUSERSITE=1",
      "PYTHONPATH=",
      "OPENBLAS_NUM_THREADS=1",
      "OMP_NUM_THREADS=1",
      "MKL_NUM_THREADS=1"
    )
  )
  if (!identical(as.integer(rc),0L))
    stop("区域VCF提取失败，请查看：", pylog)

  if (!file.exists(file.path(extracted,"extraction_QC.tsv")))
    stop("区域提取QC文件缺失。")

  # ------------------------------------------------------------
  # 7. harmonization函数
  # ------------------------------------------------------------
  comp <- function(a) chartr("ACGT","TGCA",a)

  harmonize_region <- function(region_id) {

    f1 <- file.path(extracted,
                    paste0("pancreas_fat__",region_id,".tsv.gz"))
    f2 <- file.path(extracted,
                    paste0("pancreas_iron__",region_id,".tsv.gz"))

    a <- fread(f1)
    b <- fread(f2)

    if (!nrow(a) || !nrow(b))
      stop("区域没有可用GWAS变异：", region_id)

    # 同一位置若存在多条记录，优先要求唯一；重复位置全部剔除，避免多等位/重复记录。
    dupA <- a[, .N, by=.(CHR,BP)][N>1,.(CHR,BP)]
    dupB <- b[, .N, by=.(CHR,BP)][N>1,.(CHR,BP)]
    if (nrow(dupA)) a <- a[!dupA, on=.(CHR,BP)]
    if (nrow(dupB)) b <- b[!dupB, on=.(CHR,BP)]

    m <- merge(
      a,b,
      by=c("CHR","BP"),
      suffixes=c("_fat","_iron"),
      all=FALSE
    )

    if (!nrow(m))
      stop("两表型在区域内没有共同坐标变异：", region_id)

    # 四种可接受方向
    same <- m$REF_fat==m$REF_iron & m$ALT_fat==m$ALT_iron
    swap <- m$REF_fat==m$ALT_iron & m$ALT_fat==m$REF_iron
    comp_same <- comp(m$REF_fat)==m$REF_iron & comp(m$ALT_fat)==m$ALT_iron
    comp_swap <- comp(m$REF_fat)==m$ALT_iron & comp(m$ALT_fat)==m$REF_iron

    ok <- same | swap | comp_same | comp_swap
    m <- m[ok]
    same <- same[ok]; swap <- swap[ok]
    comp_same <- comp_same[ok]; comp_swap <- comp_swap[ok]

    if (!nrow(m))
      stop("等位基因harmonization后无变异：",region_id)

    flip <- swap | comp_swap

    # 统一到fat的ALT等位基因方向
    m[, BETA_iron_aligned := ifelse(flip,-BETA_iron,BETA_iron)]
    m[, AF_iron_aligned := ifelse(flip,1-AF_iron,AF_iron)]
    m[, MAF_iron_aligned := pmin(AF_iron_aligned,1-AF_iron_aligned)]

    # SNP标识优先使用一致且唯一rsID，否则坐标等位基因
    rs <- ifelse(
      !is.na(m$SNP_fat) & m$SNP_fat!="." & m$SNP_fat!="" &
        m$SNP_fat==m$SNP_iron,
      m$SNP_fat,
      paste0("chr",m$CHR,":",m$BP,":",m$REF_fat,":",m$ALT_fat)
    )

    # 保证snp id唯一
    rs[duplicated(rs) | duplicated(rs,fromLast=TRUE)] <-
      paste0("chr",m$CHR[duplicated(rs) | duplicated(rs,fromLast=TRUE)],
             ":",m$BP[duplicated(rs) | duplicated(rs,fromLast=TRUE)],
             ":",m$REF_fat[duplicated(rs) | duplicated(rs,fromLast=TRUE)],
             ":",m$ALT_fat[duplicated(rs) | duplicated(rs,fromLast=TRUE)])

    m[, SNP_UNIQUE := rs]
    m[, harmonization := fifelse(
      same,"same",
      fifelse(swap,"swap",
              fifelse(comp_same,"complement_same","complement_swap"))
    )]

    # 两表型AF若严重不一致标记，不自动删；后续敏感性表可查看。
    m[, AF_abs_diff := abs(AF_fat-AF_iron_aligned)]

    setorder(m,BP)

    fwrite(
      m,
      file.path(extracted,paste0("harmonized__",region_id,".tsv.gz")),
      sep="\t", na="NA"
    )

    m
  }

  harm <- lapply(regions$region_id, harmonize_region)
  names(harm) <- regions$region_id

  harm_qc <- rbindlist(lapply(names(harm), function(rid) {
    d <- harm[[rid]]
    data.table(
      region_id=rid,
      n_harmonized=nrow(d),
      n_maf_ge_001=sum(pmin(d$MAF_fat,d$MAF_iron_aligned)>=0.01),
      n_AF_diff_gt_020=sum(d$AF_abs_diff>0.20),
      median_AF_diff=median(d$AF_abs_diff,na.rm=TRUE)
    )
  }))
  fwrite(harm_qc,file.path(reports,"harmonization_QC.tsv"),sep="\t")

  # ------------------------------------------------------------
  # 8. coloc.abf
  # ------------------------------------------------------------
  sdY_est_fun <- tryCatch(
    getFromNamespace("sdY.est","coloc"),
    error=function(e) NULL
  )

  get_sdY_est <- function(varbeta, maf, n) {
    if (is.null(sdY_est_fun)) return(NA_real_)
    tryCatch(
      as.numeric(sdY_est_fun(varbeta=varbeta, maf=maf, n=n)),
      error=function(e) NA_real_
    )
  }

  run_abf <- function(d, region_id, maf_mode, sdy_mode, p12) {

    x <- copy(d)

    if (maf_mode=="MAF_GE_001") {
      x <- x[
        pmin(MAF_fat,MAF_iron_aligned) >= 0.01
      ]
    }

    if (nrow(x)<100)
      return(list(status="TOO_FEW_VARIANTS",summary=NULL,results=NULL,
                  n=nrow(x),sdY1=NA_real_,sdY2=NA_real_))

    maf1 <- x$MAF_fat
    maf2 <- x$MAF_iron_aligned

    if (sdy_mode=="SDY_1") {
      sdy1 <- 1
      sdy2 <- 1
    } else {
      sdy1 <- get_sdY_est(x$VARBETA_fat,maf1,N["pancreas_fat"])
      sdy2 <- get_sdY_est(x$VARBETA_iron,maf2,N["pancreas_iron"])
      if (!is.finite(sdy1) || !is.finite(sdy2) || sdy1<=0 || sdy2<=0)
        return(list(status="SDY_EST_FAILED",summary=NULL,results=NULL,
                    n=nrow(x),sdY1=sdy1,sdY2=sdy2))
    }

    d1 <- list(
      beta=x$BETA_fat,
      varbeta=x$VARBETA_fat,
      snp=x$SNP_UNIQUE,
      position=x$BP,
      MAF=maf1,
      N=as.numeric(N["pancreas_fat"]),
      type="quant",
      sdY=sdy1
    )

    d2 <- list(
      beta=x$BETA_iron_aligned,
      varbeta=x$VARBETA_iron,
      snp=x$SNP_UNIQUE,
      position=x$BP,
      MAF=maf2,
      N=as.numeric(N["pancreas_iron"]),
      type="quant",
      sdY=sdy2
    )

    res <- tryCatch(
      coloc::coloc.abf(
        dataset1=d1,
        dataset2=d2,
        p1=1e-4,
        p2=1e-4,
        p12=p12
      ),
      error=function(e) e
    )

    if (inherits(res,"error"))
      return(list(status=paste0("ERROR: ",conditionMessage(res)),
                  summary=NULL,results=NULL,n=nrow(x),sdY1=sdy1,sdY2=sdy2))

    list(status="OK",summary=res$summary,results=res$results,
         n=nrow(x),sdY1=sdy1,sdY2=sdy2)
  }

  abf_grid <- CJ(
    region_id=regions$region_id,
    maf_mode=c("ALL_VALID","MAF_GE_001"),
    sdy_mode=c("SDY_1","SDY_EST"),
    p12=c(1e-4,1e-5,1e-6),
    unique=TRUE
  )

  abf_rows <- list()
  abf_objects <- list()

  for (i in seq_len(nrow(abf_grid))) {

    g <- abf_grid[i]
    res <- run_abf(
      harm[[g$region_id]],
      g$region_id,
      g$maf_mode,
      g$sdy_mode,
      g$p12
    )

    key <- paste(g$region_id,g$maf_mode,g$sdy_mode,
                 format(g$p12,scientific=TRUE),sep="__")
    abf_objects[[key]] <- res

    row <- data.table(
      region_id=g$region_id,
      maf_mode=g$maf_mode,
      sdy_mode=g$sdy_mode,
      p1=1e-4,
      p2=1e-4,
      p12=g$p12,
      status=res$status,
      n_variants=res$n,
      sdY_fat=res$sdY1,
      sdY_iron=res$sdY2,
      PP.H0=NA_real_,
      PP.H1=NA_real_,
      PP.H2=NA_real_,
      PP.H3=NA_real_,
      PP.H4=NA_real_
    )

    if (identical(res$status,"OK")) {
      s <- res$summary
      row[, `:=`(
        PP.H0=as.numeric(s["PP.H0.abf"]),
        PP.H1=as.numeric(s["PP.H1.abf"]),
        PP.H2=as.numeric(s["PP.H2.abf"]),
        PP.H3=as.numeric(s["PP.H3.abf"]),
        PP.H4=as.numeric(s["PP.H4.abf"])
      )]
    }

    abf_rows[[i]] <- row
  }

  abf <- rbindlist(abf_rows,fill=TRUE)
  fwrite(abf,file.path(abfdir,"Step26_coloc_abf_all_sensitivity.tsv"),
         sep="\t",na="NA")

  # ------------------------------------------------------------
  # 9. 主分析/阳性参照
  # ------------------------------------------------------------
  main_chr19 <- abf[
    region_id=="chr19_LOC2340_primary" &
      maf_mode=="MAF_GE_001" &
      sdy_mode=="SDY_1" &
      abs(p12-1e-5)<1e-12
  ]

  main_chr8 <- abf[
    region_id=="chr8_positive_control" &
      maf_mode=="MAF_GE_001" &
      sdy_mode=="SDY_1" &
      abs(p12-1e-5)<1e-12
  ]

  if (nrow(main_chr19)!=1L || nrow(main_chr8)!=1L)
    stop("无法唯一提取Step26主分析或chr8参照分析。")

  fwrite(main_chr19,file.path(abfdir,"Step26_chr19_primary_coloc.tsv"),sep="\t")
  fwrite(main_chr8,file.path(abfdir,"Step26_chr8_positive_control_coloc.tsv"),sep="\t")

  # 主分析SNP级H4 posterior
  save_top_h4 <- function(region_id,key_pattern,outfile) {
    keys <- names(abf_objects)
    k <- keys[
      grepl(paste0("^",region_id,"__MAF_GE_001__SDY_1__"),keys) &
      grepl("1e-05",tolower(keys),fixed=TRUE)
    ]
    if (length(k)!=1L) {
      # p12格式可能为1e-05或1e-05，直接用grid定位重建key
      cand <- keys[grepl(paste0("^",region_id,"__MAF_GE_001__SDY_1__"),keys)]
      pp <- vapply(cand,function(z) {
        x <- abf_objects[[z]]
        if (identical(x$status,"OK")) {
          # 通过对应abf表定位p12=1e-5的对象较麻烦；取对象结果和主分析一致的PP.H4
          as.numeric(x$summary["PP.H4.abf"])
        } else NA_real_
      },numeric(1))
      target <- if (region_id=="chr19_LOC2340_primary") main_chr19$PP.H4 else main_chr8$PP.H4
      k <- cand[which.min(abs(pp-target))]
    }
    if (length(k)!=1L || !identical(abf_objects[[k]]$status,"OK")) return(FALSE)
    rr <- as.data.table(abf_objects[[k]]$results)
    if (!"SNP.PP.H4" %in% names(rr)) return(FALSE)
    setorder(rr,-SNP.PP.H4)
    fwrite(rr,file.path(abfdir,outfile),sep="\t",na="NA")
    TRUE
  }

  # 更稳妥：直接重新运行两个main对象并写SNP结果
  get_main_obj <- function(region_id) {
    run_abf(harm[[region_id]],region_id,"MAF_GE_001","SDY_1",1e-5)
  }
  obj19 <- get_main_obj("chr19_LOC2340_primary")
  obj8  <- get_main_obj("chr8_positive_control")

  if (identical(obj19$status,"OK")) {
    rr <- as.data.table(obj19$results)
    setorder(rr,-SNP.PP.H4)
    fwrite(rr,file.path(abfdir,"Step26_chr19_primary_SNP_H4_posterior.tsv"),
           sep="\t",na="NA")
    fwrite(head(rr,20),
           file.path(abfdir,"Step26_chr19_primary_top20_H4_variants.tsv"),
           sep="\t",na="NA")
  }

  if (identical(obj8$status,"OK")) {
    rr <- as.data.table(obj8$results)
    setorder(rr,-SNP.PP.H4)
    fwrite(rr,file.path(abfdir,"Step26_chr8_control_SNP_H4_posterior.tsv"),
           sep="\t",na="NA")
    fwrite(head(rr,20),
           file.path(abfdir,"Step26_chr8_control_top20_H4_variants.tsv"),
           sep="\t",na="NA")
  }

  # ------------------------------------------------------------
  # 10. ABF稳健性摘要
  # ------------------------------------------------------------
  robustness <- abf[
    status=="OK",
    .(
      n_runs=.N,
      min_PP_H4=min(PP.H4,na.rm=TRUE),
      median_PP_H4=median(PP.H4,na.rm=TRUE),
      max_PP_H4=max(PP.H4,na.rm=TRUE),
      n_PP_H4_ge_080=sum(PP.H4>=0.80,na.rm=TRUE),
      n_PP_H4_ge_090=sum(PP.H4>=0.90,na.rm=TRUE)
    ),
    by=region_id
  ]
  fwrite(robustness,file.path(abfdir,"Step26_ABF_robustness_summary.tsv"),sep="\t")

  # ------------------------------------------------------------
  # 11. 可选SuSiE-RSS + coloc.susie（主区域和chr8参照）
  # ------------------------------------------------------------
  plink_bin <- path.expand("~/mr_tools/plink")
  ref_prefix <- file.path(base,"00_ld_reference","EUR")
  have_ref <- file.exists(plink_bin) &&
    all(file.exists(paste0(ref_prefix,c(".bed",".bim",".fam"))))

  susie_status <- list()

  run_susie_region <- function(region_id) {

    if (!have_ref)
      return(list(status="SKIPPED_NO_PLINK_REFERENCE"))

    d <- copy(harm[[region_id]])
    d <- d[pmin(MAF_fat,MAF_iron_aligned)>=0.01]

    # 仅保留唯一、可用rsID
    d <- d[
      !is.na(SNP_fat) &
      SNP_fat!="." &
      SNP_fat!="" &
      SNP_fat==SNP_iron
    ]
    d <- d[!duplicated(SNP_fat) & !duplicated(SNP_fat,fromLast=TRUE)]

    if (nrow(d)<100)
      return(list(status="SKIPPED_TOO_FEW_RSID_VARIANTS"))

    # 为避免与GWAS表中的CHR/BP字段在merge后自动变成CHR.x/CHR.y，
    # 在合并前给reference BIM字段使用明确前缀。
    bim <- fread(
      paste0(ref_prefix,".bim"),
      col.names=c("REF_CHR","SNP","REF_CM","REF_BP","REF_A1","REF_A2"),
      showProgress=FALSE
    )
    bim <- bim[!duplicated(SNP) & !duplicated(SNP,fromLast=TRUE)]
    setkey(bim,SNP)

    z <- merge(
      d,
      bim,
      by.x="SNP_fat",
      by.y="SNP",
      all=FALSE
    )

    required_merge_cols <- c(
      "SNP_fat","CHR","BP","REF_fat","ALT_fat",
      "BETA_fat","BETA_iron_aligned","VARBETA_fat","VARBETA_iron",
      "REF_CHR","REF_BP","REF_A1","REF_A2"
    )
    if (!all(required_merge_cols %in% names(z))) {
      return(list(
        status=paste0(
          "BIM_MERGE_COLUMNS_UNEXPECTED: ",
          paste(setdiff(required_merge_cols,names(z)),collapse=",")
        ),
        n_reference=nrow(z)
      ))
    }

    # GWAS harmonized table保留CHR/BP；reference使用REF_CHR/REF_BP。
    # 必须严格要求GRCh37染色体和坐标一致。
    z <- z[
      CHR == REF_CHR &
      BP == REF_BP
    ]

    # 在data.table外部构造逻辑向量时，必须显式引用z中的列；
    # 不能直接写REF_fat/REF_A1，否则R会把它们当成全局对象查找。
    allele_ok <- (
      (z$REF_fat == z$REF_A1 & z$ALT_fat == z$REF_A2) |
      (z$REF_fat == z$REF_A2 & z$ALT_fat == z$REF_A1) |
      (comp(z$REF_fat) == z$REF_A1 & comp(z$ALT_fat) == z$REF_A2) |
      (comp(z$REF_fat) == z$REF_A2 & comp(z$ALT_fat) == z$REF_A1)
    )
    z <- z[allele_ok]

    if (nrow(z)<100)
      return(list(status="SKIPPED_LOW_REFERENCE_COVERAGE",
                  n_reference=nrow(z)))

    region_dir <- file.path(susiedir,region_id)
    dir.create(region_dir,showWarnings=FALSE)

    extract_file <- file.path(region_dir,"extract.txt")
    fwrite(z[,.(SNP_fat)],extract_file,col.names=FALSE)

    subpref <- file.path(region_dir,"ref_subset")
    log1 <- file.path(region_dir,"plink_makebed.log")

    rc1 <- system2(
      plink_bin,
      c("--bfile",shQuote(ref_prefix),
        "--extract",shQuote(extract_file),
        "--make-bed",
        "--keep-allele-order",
        "--allow-no-sex",
        "--out",shQuote(subpref)),
      stdout=log1,stderr=log1
    )

    if (!identical(as.integer(rc1),0L) ||
        !file.exists(paste0(subpref,".bim")))
      return(list(status="PLINK_MAKEBED_FAILED"))

    sbim <- fread(
      paste0(subpref,".bim"),
      col.names=c("CHR","SNP","CM","BP","A1","A2")
    )
    if (nrow(sbim)<100)
      return(list(status="SKIPPED_TOO_FEW_SUBSET_VARIANTS",
                  n_reference=nrow(sbim)))

    # 以最终subset BIM顺序匹配GWAS
    setkey(z,SNP_fat)
    zz <- z[sbim$SNP]
    if (anyNA(zz$SNP_fat))
      return(list(status="SUBSET_ORDER_MATCH_FAILED"))

    # 对齐GWAS effect allele(ALT_fat)到PLINK A1
    sameA1 <- (
      zz$ALT_fat==sbim$A1 |
      comp(zz$ALT_fat)==sbim$A1
    )
    sameA2 <- (
      zz$ALT_fat==sbim$A2 |
      comp(zz$ALT_fat)==sbim$A2
    )
    if (any(!(sameA1 | sameA2)))
      return(list(status="ALLELE_TO_LD_ALIGNMENT_FAILED"))

    sign_ld <- ifelse(sameA1,1,-1)

    # LD
    ldpref <- file.path(region_dir,"ld")
    log2 <- file.path(region_dir,"plink_ld.log")
    rc2 <- system2(
      plink_bin,
      c("--bfile",shQuote(subpref),
        "--r","square","gz",
        "--allow-no-sex",
        "--out",shQuote(ldpref)),
      stdout=log2,stderr=log2
    )
    ldfile <- paste0(ldpref,".ld.gz")
    if (!identical(as.integer(rc2),0L) || !file.exists(ldfile))
      return(list(status="PLINK_LD_FAILED"))

    R <- as.matrix(read.table(gzfile(ldfile),header=FALSE,check.names=FALSE))
    storage.mode(R) <- "double"

    if (nrow(R)!=nrow(sbim) || ncol(R)!=nrow(sbim))
      return(list(status="LD_DIMENSION_MISMATCH"))

    if (any(!is.finite(R)))
      return(list(status="LD_NONFINITE"))

    # 数值对称化
    R <- (R+t(R))/2
    diag(R) <- 1

    # 对齐beta到LD A1
    beta1 <- zz$BETA_fat * sign_ld
    beta2 <- zz$BETA_iron_aligned * sign_ld
    vb1 <- zz$VARBETA_fat
    vb2 <- zz$VARBETA_iron

    D1 <- list(
      beta=beta1,
      varbeta=vb1,
      snp=sbim$SNP,
      position=sbim$BP,
      N=as.numeric(N["pancreas_fat"]),
      type="quant",
      sdY=1,
      LD=R
    )

    D2 <- list(
      beta=beta2,
      varbeta=vb2,
      snp=sbim$SNP,
      position=sbim$BP,
      N=as.numeric(N["pancreas_iron"]),
      type="quant",
      sdY=1,
      LD=R
    )

    s1 <- tryCatch(
      coloc::runsusie(D1,L=10,maxit=10000),
      error=function(e) e
    )
    if (inherits(s1,"error"))
      return(list(status=paste0("RUNSUSIE_FAT_ERROR: ",conditionMessage(s1)),
                  n_reference=nrow(sbim)))

    s2 <- tryCatch(
      coloc::runsusie(D2,L=10,maxit=10000),
      error=function(e) e
    )
    if (inherits(s2,"error"))
      return(list(status=paste0("RUNSUSIE_IRON_ERROR: ",conditionMessage(s2)),
                  n_reference=nrow(sbim)))

    cs <- tryCatch(
      coloc::coloc.susie(s1,s2,p1=1e-4,p2=1e-4,p12=1e-5),
      error=function(e) e
    )
    if (inherits(cs,"error"))
      return(list(status=paste0("COLOC_SUSIE_ERROR: ",conditionMessage(cs)),
                  n_reference=nrow(sbim)))

    saveRDS(
      list(D1=D1,D2=D2,susie1=s1,susie2=s2,coloc_susie=cs),
      file.path(region_dir,"susie_objects.rds")
    )

    if (!is.null(cs$summary)) {
      fwrite(as.data.table(cs$summary),
             file.path(region_dir,"coloc_susie_summary.tsv"),
             sep="\t",na="NA")
    }

    # credible set摘要
    cs1 <- tryCatch(susieR::susie_get_cs(s1),error=function(e) NULL)
    cs2 <- tryCatch(susieR::susie_get_cs(s2),error=function(e) NULL)

    if (!is.null(cs1)) saveRDS(cs1,file.path(region_dir,"fat_credible_sets.rds"))
    if (!is.null(cs2)) saveRDS(cs2,file.path(region_dir,"iron_credible_sets.rds"))

    list(
      status="OK",
      n_harmonized=nrow(d),
      n_reference=nrow(sbim),
      reference_coverage=nrow(sbim)/nrow(d),
      n_cs_fat=if (!is.null(cs1$cs)) length(cs1$cs) else NA_integer_,
      n_cs_iron=if (!is.null(cs2$cs)) length(cs2$cs) else NA_integer_
    )
  }

  for (rid in c("chr19_LOC2340_primary","chr8_positive_control")) {
    susie_status[[rid]] <- run_susie_region(rid)
  }

  susie_tab <- rbindlist(lapply(names(susie_status),function(rid) {
    x <- susie_status[[rid]]
    as.data.table(c(list(region_id=rid),x))
  }),fill=TRUE)

  fwrite(susie_tab,file.path(susiedir,"Step26_susie_status.tsv"),
         sep="\t",na="NA")

  # ------------------------------------------------------------
  # 12. 核心结果摘要
  # ------------------------------------------------------------
  summary <- data.table(
    item=c(
      "chr19_primary_PP_H4",
      "chr19_primary_PP_H3",
      "chr19_primary_n_variants",
      "chr8_control_PP_H4",
      "chr8_control_PP_H3",
      "chr8_control_n_variants",
      "chr19_sensitivity_min_PP_H4",
      "chr19_sensitivity_max_PP_H4"
    ),
    value=c(
      main_chr19$PP.H4,
      main_chr19$PP.H3,
      main_chr19$n_variants,
      main_chr8$PP.H4,
      main_chr8$PP.H3,
      main_chr8$n_variants,
      min(abf[
        grepl("^chr19",region_id) & status=="OK",
        PP.H4
      ],na.rm=TRUE),
      max(abf[
        grepl("^chr19",region_id) & status=="OK",
        PP.H4
      ],na.rm=TRUE)
    )
  )

  fwrite(summary,file.path(reports,"Step26_core_summary.tsv"),sep="\t")

  # ------------------------------------------------------------
  # 13. 自动但克制的解释标签
  # ------------------------------------------------------------
  classify_h4 <- function(pp) {
    if (!is.finite(pp)) return("NOT_ESTIMABLE")
    if (pp>=0.90) return("STRONG_H4_SUPPORT")
    if (pp>=0.80) return("MODERATE_TO_STRONG_H4_SUPPORT")
    if (pp>=0.50) return("INTERMEDIATE_H4_SUPPORT")
    return("LIMITED_H4_SUPPORT")
  }

  interpretation <- data.table(
    region=c("chr19_LOC2340_primary","chr8_positive_control"),
    PP_H4=c(main_chr19$PP.H4,main_chr8$PP.H4),
    evidence_label=c(
      classify_h4(main_chr19$PP.H4),
      classify_h4(main_chr8$PP.H4)
    ),
    caution=c(
      "chr19 was selected from exploratory Step25b local-rg sensitivity analysis; H4 does not imply causality.",
      "chr8 is a positive-control/reproducibility region previously analysed; H4 does not imply causality."
    )
  )
  fwrite(interpretation,file.path(reports,"Step26_interpretation_labels.tsv"),
         sep="\t")

  # ------------------------------------------------------------
  # 14. 方法与解释边界
  # ------------------------------------------------------------
  writeLines(c(
    "Step26 tests cross-trait colocalization between pancreatic fat and pancreatic iron.",
    "Primary candidate: LAVA LOC2340, GRCh37 chr19:33785836-34633274, identified in Step25b exploratory local-rg sensitivity analysis.",
    "Positive-control region: GRCh37 chr8:129047378-130047378, previously showing strong pancreatic imaging-trait colocalization.",
    "Raw regional GWAS summary statistics are extracted from the original VCFs without association-P filtering.",
    "Primary ABF analysis uses MAF>=0.01, sdY=1, p1=p2=1e-4 and p12=1e-5.",
    "Sensitivity analyses vary region boundaries, MAF restriction, sdY specification and p12.",
    "When the existing 1000G EUR PLINK reference is available, SuSiE-RSS and coloc.susie are additionally attempted for the chr19 primary and chr8 control regions.",
    "ABF PP.H4 is posterior support for a shared causal-variant model within the analysed region under the coloc assumptions and chosen priors.",
    "Neither PP.H4 nor local genetic correlation demonstrates that pancreatic fat causes pancreatic iron or vice versa.",
    "Because chr19 was prioritized from the same GWAS pair in Step25b, Step26 is mechanistic follow-up rather than independent replication."
  ),file.path(reports,"methods_and_interpretation.txt"))

  # ------------------------------------------------------------
  # 15. 完整性检查
  # ------------------------------------------------------------
  if (!identical(main_chr19$status,"OK"))
    stop("chr19主分析coloc.abf失败：",main_chr19$status)

  if (!identical(main_chr8$status,"OK"))
    stop("chr8阳性参照coloc.abf失败：",main_chr8$status)

  if (!is.finite(main_chr19$PP.H4) || !is.finite(main_chr8$PP.H4))
    stop("主分析PP.H4不是有限数值。")

  writeLines(
    "COMPLETED_STEP26_CROSS_TRAIT_COLOCALIZATION_REVIEW_RESULTS",
    file.path(reports,"completion_status.txt")
  )

  # ------------------------------------------------------------
  # 16. 打包
  # ------------------------------------------------------------
  zipfile <- file.path(out,"26_colocalization_review.zip")
  oldwd <- getwd()
  on.exit(setwd(oldwd),add=TRUE)
  setwd(out)

  files_to_zip <- c(
    list.files("reports",full.names=TRUE),
    list.files("coloc_abf",full.names=TRUE),
    list.files("coloc_susie",recursive=TRUE,full.names=TRUE),
    "regional_sumstats/extraction_QC.tsv"
  )
  files_to_zip <- files_to_zip[file.exists(files_to_zip)]

  utils::zip(zipfile,files=files_to_zip,flags="-q")

  # ------------------------------------------------------------
  # 17. 控制台摘要
  # ------------------------------------------------------------
  cat("\n===== 第26步完成 =====\n")
  cat("\nchr19 LOC2340主分析：\n")
  print(main_chr19)

  cat("\nchr8阳性参照：\n")
  print(main_chr8)

  cat("\nABF稳健性摘要：\n")
  print(robustness)

  cat("\nSuSiE状态：\n")
  print(susie_tab)

  cat("\n结果目录：",out,"\n")
  cat("复核包：",zipfile,"\n")
})()
