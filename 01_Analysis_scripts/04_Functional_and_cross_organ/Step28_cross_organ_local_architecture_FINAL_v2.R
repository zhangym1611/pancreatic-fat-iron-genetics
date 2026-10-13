# ============================================================
# 第28步：chr8 / chr19跨器官局部共享遗传结构整合
# 修正版v2：修复harmonize_two()中SNP_1/SNP_2/CHR/BP/REF_1/ALT_1在data.table外部未显式引用的问题。
# Cross-organ locus-specific colocalization integration
#
# 核心问题：
# 1) chr8核心共享位点是否主要体现胰腺脂肪-胰腺铁共享，
#    还是也广泛参与肝脂肪/肝铁/VAT/ASAT？
# 2) chr19次级共享位点是否为胰腺局部特异候选，
#    还是属于更广泛的跨器官代谢共享信号？
# 3) 与同器官“肝脂肪-肝铁”对照相比，胰腺脂肪-胰腺铁的局部共享
#    是否具有不同的遗传结构？
#
# 预设表型：
#   pancreas_fat   GCST90016675   N=25,617
#   pancreas_iron  GCST90016676   N=25,617
#   liver_fat      GCST90016673   N=32,858
#   liver_iron     GCST90016674   N=32,858
#   VAT            GCST90016671   N=32,860
#   ASAT           GCST90016672   N=32,860
#
# 预设pair：
#   pancreas_fat  vs pancreas_iron  （胰腺anchor）
#   liver_fat     vs liver_iron     （肝脏同器官对照）
#   pancreas_fat  vs liver_fat      （跨器官脂肪）
#   pancreas_iron vs liver_iron     （跨器官铁）
#   pancreas_fat  vs VAT            （内脏脂肪背景）
#   pancreas_fat  vs ASAT           （皮下脂肪背景）
#
# 区域：
#   chr8  : GRCh37 chr8:129047378-130047378
#   chr19 : GRCh37 chr19:33785836-34633274 (LAVA LOC2340)
#
# 方法：
# - 原始GWAS VCF按区域提取，不做P值筛选
# - biallelic SNV，FILTER PASS/.，有效ES/SE/AF
# - 主分析 MAF>=0.01
# - coloc.abf, quantitative trait, sdY=1
# - p1=p2=1e-4
# - p12主分析=1e-5
# - p12敏感性=1e-4,1e-5,1e-6
#
# 解释边界：
# - PP.H4支持“同一区域共享因果变异模型”，不代表因果方向；
# - 跨器官H4高不等于同一基因/同一细胞类型机制；
# - chr19来自Step25b探索性LAVA筛选，因此仍为secondary/exploratory locus；
# - 本步是局部位点整合，不替代Step22/23的全基因组rg结果。
# ============================================================

(function() {

  options(stringsAsFactors = FALSE)
  base <- path.expand("~/胰腺")

  # ------------------------------------------------------------
  # 0. R依赖
  # ------------------------------------------------------------
  lib <- path.expand("~/R_pancreas_step28")
  dir.create(lib, recursive=TRUE, showWarnings=FALSE)
  .libPaths(c(lib, .libPaths()))

  need <- function(pkg) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
      install.packages(pkg, lib=lib, repos="https://cloud.r-project.org")
    }
    if (!requireNamespace(pkg, quietly=TRUE)) {
      stop("无法加载R包：", pkg)
    }
  }

  for (p in c("data.table","coloc")) need(p)

  suppressPackageStartupMessages(library(data.table))

  cat("coloc version: ",
      as.character(utils::packageVersion("coloc")),
      "\n", sep="")

  # ------------------------------------------------------------
  # 1. 固定表型信息
  # ------------------------------------------------------------
  trait_info <- data.table(
    trait=c(
      "pancreas_fat",
      "pancreas_iron",
      "liver_fat",
      "liver_iron",
      "VAT",
      "ASAT"
    ),
    study_id=c(
      "ebi-a-GCST90016675",
      "ebi-a-GCST90016676",
      "ebi-a-GCST90016673",
      "ebi-a-GCST90016674",
      "ebi-a-GCST90016671",
      "ebi-a-GCST90016672"
    ),
    N=c(
      25617,
      25617,
      32858,
      32858,
      32860,
      32860
    ),
    category=c(
      "pancreas",
      "pancreas",
      "liver",
      "liver",
      "adipose",
      "adipose"
    )
  )

  trait_info[, vcf := file.path(
    base,
    paste0(study_id, ".vcf.gz")
  )]

  if (!all(file.exists(trait_info$vcf))) {
    stop(
      "缺少以下原始GWAS VCF：\n",
      paste(
        trait_info[!file.exists(vcf), vcf],
        collapse="\n"
      )
    )
  }

  # ------------------------------------------------------------
  # 2. 固定区域
  # ------------------------------------------------------------
  regions <- data.table(
    region_id=c(
      "chr8_primary_shared_locus",
      "chr19_secondary_shared_locus"
    ),
    CHR=c(8L,19L),
    START=c(
      129047378L,
      33785836L
    ),
    STOP=c(
      130047378L,
      34633274L
    ),
    role=c(
      "PRIMARY_SHARED_LOCUS",
      "EXPLORATORY_SECONDARY_SHARED_LOCUS"
    )
  )

  # ------------------------------------------------------------
  # 3. 固定pair
  # ------------------------------------------------------------
  pairs <- data.table(
    pair_id=c(
      "pancreas_fat__pancreas_iron",
      "liver_fat__liver_iron",
      "pancreas_fat__liver_fat",
      "pancreas_iron__liver_iron",
      "pancreas_fat__VAT",
      "pancreas_fat__ASAT"
    ),
    trait1=c(
      "pancreas_fat",
      "liver_fat",
      "pancreas_fat",
      "pancreas_iron",
      "pancreas_fat",
      "pancreas_fat"
    ),
    trait2=c(
      "pancreas_iron",
      "liver_iron",
      "liver_fat",
      "liver_iron",
      "VAT",
      "ASAT"
    ),
    comparison_type=c(
      "PANCREAS_WITHIN_ORGAN_ANCHOR",
      "LIVER_WITHIN_ORGAN_COMPARATOR",
      "CROSS_ORGAN_FAT",
      "CROSS_ORGAN_IRON",
      "PANCREAS_FAT_VISCERAL_ADIPOSE",
      "PANCREAS_FAT_SUBCUTANEOUS_ADIPOSE"
    )
  )

  # ------------------------------------------------------------
  # 4. 输出目录
  # ------------------------------------------------------------
  parent <- file.path(base, "28_cross_organ_local_architecture")
  dir.create(parent, recursive=TRUE, showWarnings=FALSE)

  out <- tempfile(
    "chr8_chr19_cross_organ_",
    tmpdir=parent
  )

  if (!dir.create(out)) {
    stop("无法创建Step28输出目录：", out)
  }

  reports <- file.path(out,"reports")
  extracted <- file.path(out,"regional_sumstats")
  colocdir <- file.path(out,"coloc")
  logs <- file.path(out,"logs")

  for (d in c(reports,extracted,colocdir,logs))
    dir.create(d,recursive=TRUE,showWarnings=FALSE)

  writeLines(
    "INCOMPLETE_CHECK_LOGS",
    file.path(reports,"completion_status.txt")
  )

  fwrite(
    trait_info[,.(trait,study_id,N,category)],
    file.path(reports,"Step28_trait_info.tsv"),
    sep="\t"
  )

  fwrite(
    regions,
    file.path(reports,"Step28_regions.tsv"),
    sep="\t"
  )

  fwrite(
    pairs,
    file.path(reports,"Step28_pairs.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 5. Python环境
  #    优先复用Step19固定Python；如果不存在则系统python3
  # ------------------------------------------------------------
  python_bin <- ""

  step19_cfg_candidates <- list.files(
    base,
    pattern="LDSC_paths\\.rds$",
    recursive=TRUE,
    full.names=TRUE
  )

  if (length(step19_cfg_candidates)) {

    # 优先路径中含19_LDSC_setup
    pref <- step19_cfg_candidates[
      grepl("19_LDSC_setup",step19_cfg_candidates,fixed=TRUE)
    ]

    cfgfile <- if (length(pref)) pref[1] else step19_cfg_candidates[1]

    cfg <- tryCatch(
      readRDS(cfgfile),
      error=function(e) NULL
    )

    if (!is.null(cfg) &&
        !is.null(cfg$python) &&
        file.exists(cfg$python)) {
      python_bin <- cfg$python
    }
  }

  if (!nzchar(python_bin)) {
    python_bin <- Sys.which("python3")
  }

  if (!nzchar(python_bin)) {
    python_bin <- Sys.which("python")
  }

  if (!nzchar(python_bin)) {
    stop("未找到可用Python。")
  }

  cat("Python: ",python_bin,"\n",sep="")

  # ------------------------------------------------------------
  # 6. Python纯标准库区域提取
  #    不依赖pandas，避免环境问题
  # ------------------------------------------------------------
  pyfile <- file.path(out,"extract_step28_regions.py")

  py <- c(
'import gzip, sys, math, csv, os, traceback',
'',
'base=sys.argv[1]',
'outdir=sys.argv[2]',
'trait_file=sys.argv[3]',
'region_file=sys.argv[4]',
'',
'def finite_float(x):',
'    try:',
'        v=float(x)',
'        return v if math.isfinite(v) else None',
'    except Exception:',
'        return None',
'',
'def parse_sample(fmt, sample):',
'    k=fmt.split(":")',
'    v=sample.split(":")',
'    if len(k)!=len(v):',
'        return None',
'    return dict(zip(k,v))',
'',
'# read traits',
'traits=[]',
'with open(trait_file,newline="") as f:',
'    r=csv.DictReader(f,delimiter="\\t")',
'    for row in r:',
'        traits.append(row)',
'',
'# read regions',
'regions=[]',
'with open(region_file,newline="") as f:',
'    r=csv.DictReader(f,delimiter="\\t")',
'    for row in r:',
'        row["CHR"]=int(row["CHR"])',
'        row["START"]=int(row["START"])',
'        row["STOP"]=int(row["STOP"])',
'        regions.append(row)',
'',
'qc=[]',
'',
'for tr in traits:',
'    trait=tr["trait"]',
'    vcf=tr["vcf"]',
'',
'    handles={}',
'    writers={}',
'    counts={}',
'',
'    for rg in regions:',
'        rid=rg["region_id"]',
'        path=os.path.join(outdir, trait+"__"+rid+".tsv")',
'        fh=open(path,"w",newline="")',
'        w=csv.writer(fh,delimiter="\\t")',
'        w.writerow(["trait","region_id","CHR","BP","SNP","REF","ALT","BETA","SE","VARBETA","AF","MAF"])',
'        handles[rid]=fh',
'        writers[rid]=w',
'        counts[rid]={"raw":0,"eligible":0}',
'',
'    with gzip.open(vcf,"rt") as f:',
'        for line in f:',
'            if line.startswith("#"):',
'                continue',
'            p=line.rstrip("\\n").split("\\t")',
'            if len(p)<10:',
'                continue',
'',
'            chrom=p[0].replace("chr","")',
'            if not chrom.isdigit():',
'                continue',
'            ch=int(chrom)',
'            try:',
'                bp=int(p[1])',
'            except Exception:',
'                continue',
'',
'            hits=[rg for rg in regions if rg["CHR"]==ch and rg["START"]<=bp<=rg["STOP"]]',
'            if not hits:',
'                continue',
'',
'            for rg in hits:',
'                counts[rg["region_id"]]["raw"]+=1',
'',
'            ref=p[3].upper()',
'            alt=p[4].upper()',
'            filt=p[6]',
'',
'            if len(ref)!=1 or len(alt)!=1 or ref not in "ACGT" or alt not in "ACGT":',
'                continue',
'            if filt not in ("PASS","."):',
'                continue',
'',
'            sm=parse_sample(p[8],p[9])',
'            if sm is None:',
'                continue',
'            if not all(k in sm for k in ("ES","SE","AF")):',
'                continue',
'',
'            beta=finite_float(sm["ES"])',
'            se=finite_float(sm["SE"])',
'            af=finite_float(sm["AF"])',
'',
'            if beta is None or se is None or af is None:',
'                continue',
'            if se<=0 or not (0<af<1):',
'                continue',
'',
'            maf=min(af,1-af)',
'            rsid=p[2]',
'',
'            for rg in hits:',
'                rid=rg["region_id"]',
'                writers[rid].writerow([',
'                    trait,rid,ch,bp,rsid,ref,alt,beta,se,se*se,af,maf',
'                ])',
'                counts[rid]["eligible"]+=1',
'',
'    for rg in regions:',
'        rid=rg["region_id"]',
'        handles[rid].close()',
'        qc.append([',
'            trait,rid,counts[rid]["raw"],counts[rid]["eligible"]',
'        ])',
'',
'with open(os.path.join(outdir,"Step28_extraction_QC.tsv"),"w",newline="") as f:',
'    w=csv.writer(f,delimiter="\\t")',
'    w.writerow(["trait","region_id","raw_in_region","eligible"])',
'    w.writerows(qc)'
  )

  writeLines(py,pyfile,useBytes=TRUE)

  trait_py_file <- file.path(out,"traits_for_python.tsv")
  fwrite(
    trait_info[,.(trait,vcf)],
    trait_py_file,
    sep="\t"
  )

  region_py_file <- file.path(out,"regions_for_python.tsv")
  fwrite(
    regions[,.(region_id,CHR,START,STOP)],
    region_py_file,
    sep="\t"
  )

  pylog <- file.path(logs,"regional_extraction.log")

  rc <- system2(
    python_bin,
    args=vapply(
      c(pyfile,base,extracted,trait_py_file,region_py_file),
      shQuote,
      character(1)
    ),
    stdout=pylog,
    stderr=pylog,
    env=c(
      "PYTHONNOUSERSITE=1",
      "PYTHONPATH="
    )
  )

  if (!identical(as.integer(rc),0L)) {
    stop("Step28区域提取失败，请查看：",pylog)
  }

  eqc_file <- file.path(extracted,"Step28_extraction_QC.tsv")
  if (!file.exists(eqc_file)) {
    stop("区域提取QC文件缺失。")
  }

  extraction_qc <- fread(eqc_file)
  fwrite(
    extraction_qc,
    file.path(reports,"Step28_extraction_QC.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 7. 读区域数据
  # ------------------------------------------------------------
  regional <- list()

  for (tr in trait_info$trait) {
    for (rid in regions$region_id) {

      f <- file.path(
        extracted,
        paste0(tr,"__",rid,".tsv")
      )

      if (!file.exists(f)) {
        stop("缺少区域文件：",f)
      }

      x <- fread(f)

      if (!nrow(x)) {
        stop("区域文件没有合格变异：",f)
      }

      regional[[paste(tr,rid,sep="__")]] <- x
    }
  }

  # ------------------------------------------------------------
  # 8. harmonization
  # ------------------------------------------------------------
  comp <- function(a) {
    chartr("ACGT","TGCA",toupper(a))
  }

  harmonize_two <- function(a,b,trait1,trait2,region_id) {

    # 坐标重复全部剔除
    da <- a[, .N, by=.(CHR,BP)][N>1,.(CHR,BP)]
    db <- b[, .N, by=.(CHR,BP)][N>1,.(CHR,BP)]

    if (nrow(da)) a <- a[!da,on=.(CHR,BP)]
    if (nrow(db)) b <- b[!db,on=.(CHR,BP)]

    m <- merge(
      a,b,
      by=c("CHR","BP"),
      suffixes=c("_1","_2"),
      all=FALSE
    )

    if (!nrow(m)) return(m)

    same <- (
      m$REF_1==m$REF_2 &
      m$ALT_1==m$ALT_2
    )

    swap <- (
      m$REF_1==m$ALT_2 &
      m$ALT_1==m$REF_2
    )

    comp_same <- (
      comp(m$REF_1)==m$REF_2 &
      comp(m$ALT_1)==m$ALT_2
    )

    comp_swap <- (
      comp(m$REF_1)==m$ALT_2 &
      comp(m$ALT_1)==m$REF_2
    )

    ok <- same | swap | comp_same | comp_swap

    m <- m[ok]
    same <- same[ok]
    swap <- swap[ok]
    comp_same <- comp_same[ok]
    comp_swap <- comp_swap[ok]

    if (!nrow(m)) return(m)

    flip <- swap | comp_swap

    m[, BETA_2_aligned :=
      ifelse(flip,-BETA_2,BETA_2)
    ]

    m[, AF_2_aligned :=
      ifelse(flip,1-AF_2,AF_2)
    ]

    m[, MAF_2_aligned :=
      pmin(AF_2_aligned,1-AF_2_aligned)
    ]

    # 为coloc生成统一唯一SNP ID
    # 注意：这里位于data.table的[]环境之外，必须显式使用m$列名。
    rs_same <- (
      !is.na(m$SNP_1) &
      m$SNP_1!="." &
      m$SNP_1!="" &
      m$SNP_1==m$SNP_2
    )

    snp <- ifelse(
      rs_same,
      m$SNP_1,
      paste0(
        "chr",m$CHR,":",m$BP,":",m$REF_1,":",m$ALT_1
      )
    )

    # 最终ID必须唯一
    dup <- duplicated(snp) | duplicated(snp,fromLast=TRUE)

    if (any(dup)) {
      snp[dup] <- paste0(
        "chr",
        m$CHR[dup],":",
        m$BP[dup],":",
        m$REF_1[dup],":",
        m$ALT_1[dup],":",
        seq_len(sum(dup))
      )
    }

    m[, SNP_UNIQUE := snp]

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

    m[, common_MAF_ge_001 :=
      MAF_1>=0.01 &
      MAF_2_aligned>=0.01
    ]

    m[, trait1 := trait1]
    m[, trait2 := trait2]
    m[, region_id := region_id]

    setorder(m,BP)

    m
  }

  harm_list <- list()
  harm_qc <- list()

  for (i in seq_len(nrow(pairs))) {

    pr <- pairs[i]

    for (rid in regions$region_id) {

      a <- regional[[paste(pr$trait1,rid,sep="__")]]
      b <- regional[[paste(pr$trait2,rid,sep="__")]]

      h <- harmonize_two(
        a,b,
        pr$trait1,
        pr$trait2,
        rid
      )

      key <- paste(pr$pair_id,rid,sep="__")
      harm_list[[key]] <- h

      harm_qc[[key]] <- data.table(
        pair_id=pr$pair_id,
        comparison_type=pr$comparison_type,
        region_id=rid,
        trait1=pr$trait1,
        trait2=pr$trait2,
        n_trait1=nrow(a),
        n_trait2=nrow(b),
        n_common_harmonized=nrow(h),
        n_common_MAF_ge_001=if (nrow(h))
          sum(h$common_MAF_ge_001) else 0L,
        n_same=if (nrow(h))
          sum(h$harmonization=="same") else 0L,
        n_swap=if (nrow(h))
          sum(h$harmonization=="swap") else 0L,
        n_comp_same=if (nrow(h))
          sum(h$harmonization=="complement_same") else 0L,
        n_comp_swap=if (nrow(h))
          sum(h$harmonization=="complement_swap") else 0L
      )

      if (nrow(h)) {
        fwrite(
          h,
          file.path(
            extracted,
            paste0(
              "harmonized__",
              pr$pair_id,"__",
              rid,".tsv.gz"
            )
          ),
          sep="\t",
          na="NA"
        )
      }
    }
  }

  harm_qc_tab <- rbindlist(harm_qc,fill=TRUE)

  fwrite(
    harm_qc_tab,
    file.path(reports,"Step28_harmonization_QC.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 9. coloc.abf
  # ------------------------------------------------------------
  getN <- function(tr) {
    trait_info[trait==tr,N][1]
  }

  run_coloc <- function(h,trait1,trait2,p12) {

    x <- h[
      common_MAF_ge_001==TRUE
    ]

    if (nrow(x)<100) {
      return(list(
        status="TOO_FEW_COMMON_VARIANTS",
        n=nrow(x)
      ))
    }

    D1 <- list(
      beta=x$BETA_1,
      varbeta=x$VARBETA_1,
      snp=x$SNP_UNIQUE,
      position=x$BP,
      MAF=x$MAF_1,
      N=as.numeric(getN(trait1)),
      type="quant",
      sdY=1
    )

    D2 <- list(
      beta=x$BETA_2_aligned,
      varbeta=x$VARBETA_2,
      snp=x$SNP_UNIQUE,
      position=x$BP,
      MAF=x$MAF_2_aligned,
      N=as.numeric(getN(trait2)),
      type="quant",
      sdY=1
    )

    res <- tryCatch(
      coloc::coloc.abf(
        dataset1=D1,
        dataset2=D2,
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
        n=nrow(x)
      ))
    }

    list(
      status="OK",
      n=nrow(x),
      summary=res$summary,
      results=res$results
    )
  }

  grid <- CJ(
    pair_id=pairs$pair_id,
    region_id=regions$region_id,
    p12=c(1e-4,1e-5,1e-6),
    unique=TRUE
  )

  coloc_rows <- list()
  coloc_objects <- list()

  for (i in seq_len(nrow(grid))) {

    g <- grid[i]

    pr <- pairs[pair_id==g$pair_id][1]

    key <- paste(
      g$pair_id,
      g$region_id,
      sep="__"
    )

    rr <- run_coloc(
      harm_list[[key]],
      pr$trait1,
      pr$trait2,
      g$p12
    )

    objkey <- paste(
      g$pair_id,
      g$region_id,
      format(g$p12,scientific=TRUE),
      sep="__"
    )

    coloc_objects[[objkey]] <- rr

    row <- data.table(
      pair_id=g$pair_id,
      comparison_type=pr$comparison_type,
      trait1=pr$trait1,
      trait2=pr$trait2,
      region_id=g$region_id,
      p1=1e-4,
      p2=1e-4,
      p12=g$p12,
      status=rr$status,
      n_common_variants=if (!is.null(rr$n)) rr$n else NA_integer_,
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

  coloc_all <- rbindlist(
    coloc_rows,
    fill=TRUE
  )

  fwrite(
    coloc_all,
    file.path(
      colocdir,
      "Step28_cross_organ_coloc_all_sensitivity.tsv"
    ),
    sep="\t",
    na="NA"
  )

  # ------------------------------------------------------------
  # 10. 主分析 p12=1e-5
  # ------------------------------------------------------------
  main <- coloc_all[
    abs(p12-1e-5)<1e-12
  ]

  fwrite(
    main,
    file.path(
      colocdir,
      "Step28_cross_organ_coloc_main.tsv"
    ),
    sep="\t",
    na="NA"
  )

  # ------------------------------------------------------------
  # 11. 主分析PP.H4 matrix
  # ------------------------------------------------------------
  main_matrix <- dcast(
    main,
    pair_id + comparison_type + trait1 + trait2 ~ region_id,
    value.var="PP.H4"
  )

  fwrite(
    main_matrix,
    file.path(
      reports,
      "Step28_main_PP_H4_matrix.tsv"
    ),
    sep="\t",
    na="NA"
  )

  # ------------------------------------------------------------
  # 12. p12稳健性
  # ------------------------------------------------------------
  robustness <- coloc_all[
    status=="OK",
    .(
      min_PP_H4=min(PP.H4,na.rm=TRUE),
      median_PP_H4=median(PP.H4,na.rm=TRUE),
      max_PP_H4=max(PP.H4,na.rm=TRUE),
      n_PP_H4_ge_080=sum(PP.H4>=0.80,na.rm=TRUE),
      n_PP_H4_ge_090=sum(PP.H4>=0.90,na.rm=TRUE),
      max_PP_H3=max(PP.H3,na.rm=TRUE)
    ),
    by=.(
      pair_id,
      comparison_type,
      trait1,
      trait2,
      region_id
    )
  ]

  fwrite(
    robustness,
    file.path(
      colocdir,
      "Step28_cross_organ_coloc_robustness.tsv"
    ),
    sep="\t",
    na="NA"
  )

  # ------------------------------------------------------------
  # 13. 主分析SNP-level H4 posterior
  # ------------------------------------------------------------
  post_rows <- list()

  for (i in seq_len(nrow(main))) {

    if (main$status[i]!="OK") next

    pr <- pairs[pair_id==main$pair_id[i]][1]

    key <- paste(
      main$pair_id[i],
      main$region_id[i],
      sep="__"
    )

    rr <- run_coloc(
      harm_list[[key]],
      pr$trait1,
      pr$trait2,
      1e-5
    )

    if (!identical(rr$status,"OK")) next

    sr <- as.data.table(rr$results)

    if (!"SNP.PP.H4" %in% names(sr)) next

    sr[, pair_id := main$pair_id[i]]
    sr[, comparison_type := pr$comparison_type]
    sr[, trait1 := pr$trait1]
    sr[, trait2 := pr$trait2]
    sr[, region_id := main$region_id[i]]

    setorder(sr,-SNP.PP.H4)

    post_rows[[paste(main$pair_id[i],
                     main$region_id[i],sep="__")]] <- sr
  }

  post <- if (length(post_rows))
    rbindlist(post_rows,fill=TRUE) else data.table()

  if (nrow(post)) {

    fwrite(
      post,
      file.path(
        colocdir,
        "Step28_main_SNP_H4_posterior.tsv.gz"
      ),
      sep="\t",
      na="NA"
    )

    top20 <- post[
      order(-SNP.PP.H4),
      head(.SD,20),
      by=.(
        pair_id,
        comparison_type,
        trait1,
        trait2,
        region_id
      )
    ]

    fwrite(
      top20,
      file.path(
        colocdir,
        "Step28_main_top20_H4_variants.tsv"
      ),
      sep="\t",
      na="NA"
    )
  }

  # ------------------------------------------------------------
  # 14. 核心候选SNP核对
  # ------------------------------------------------------------
  known <- data.table(
    region_id=c(
      rep("chr8_primary_shared_locus",5),
      rep("chr19_secondary_shared_locus",2)
    ),
    SNP=c(
      "rs4733612",
      "rs1561927",
      "rs7824269",
      "rs7819154",
      "rs12676584",
      "rs33824",
      "rs430752"
    )
  )

  if (nrow(post)) {

    known_post <- merge(
      known,
      post[
        ,
        .(
          pair_id,
          comparison_type,
          trait1,
          trait2,
          region_id,
          snp,
          SNP.PP.H4
        )
      ],
      by.x=c("region_id","SNP"),
      by.y=c("region_id","snp"),
      all.x=TRUE
    )

    fwrite(
      known_post,
      file.path(
        reports,
        "Step28_known_candidate_SNP_H4_posteriors.tsv"
      ),
      sep="\t",
      na="NA"
    )
  }

  # ------------------------------------------------------------
  # 15. 区域级整合解释
  # ------------------------------------------------------------
  classify_pair <- function(h4,h3) {

    if (!is.finite(h4)) return("NOT_ESTIMABLE")

    if (h4>=0.90)
      return("STRONG_SHARED_SIGNAL_SUPPORT")

    if (h4>=0.80)
      return("MODERATE_TO_STRONG_SHARED_SIGNAL_SUPPORT")

    if (h4>=0.50)
      return("INTERMEDIATE_SHARED_SIGNAL_SUPPORT")

    if (is.finite(h3) && h3>=0.80)
      return("DISTINCT_SIGNAL_MODEL_DOMINANT")

    "LIMITED_SHARED_SIGNAL_SUPPORT"
  }

  main[, evidence_label :=
    mapply(
      classify_pair,
      PP.H4,
      PP.H3
    )
  ]

  fwrite(
    main,
    file.path(
      reports,
      "Step28_main_evidence_labels.tsv"
    ),
    sep="\t",
    na="NA"
  )

  # 按region汇总：
  # anchor = pancreas fat-pancreas iron
  # cross-organ/background = 其余5个pairs
  region_summary <- rbindlist(
    lapply(regions$region_id, function(rid) {

      x <- main[
        region_id==rid &
        status=="OK"
      ]

      anchor <- x[
        pair_id=="pancreas_fat__pancreas_iron"
      ]

      cross <- x[
        pair_id!="pancreas_fat__pancreas_iron"
      ]

      data.table(
        region_id=rid,
        anchor_PP_H4=if (nrow(anchor))
          anchor$PP.H4[1] else NA_real_,
        anchor_PP_H3=if (nrow(anchor))
          anchor$PP.H3[1] else NA_real_,
        n_non_anchor_pairs=nrow(cross),
        n_non_anchor_H4_ge_050=sum(cross$PP.H4>=0.50,na.rm=TRUE),
        n_non_anchor_H4_ge_080=sum(cross$PP.H4>=0.80,na.rm=TRUE),
        n_non_anchor_H4_ge_090=sum(cross$PP.H4>=0.90,na.rm=TRUE),
        max_non_anchor_PP_H4=if (nrow(cross))
          max(cross$PP.H4,na.rm=TRUE) else NA_real_,
        max_non_anchor_pair=if (
          nrow(cross) &&
          any(is.finite(cross$PP.H4))
        ) {
          cross$pair_id[
            which.max(cross$PP.H4)
          ]
        } else NA_character_
      )
    }),
    fill=TRUE
  )

  region_summary <- merge(
    region_summary,
    regions[,.(region_id,role)],
    by="region_id",
    all.x=TRUE
  )

  # descriptive label，不做因果断言
  region_summary[, locus_pattern :=
    fifelse(
      anchor_PP_H4>=0.80 &
        n_non_anchor_H4_ge_080==0,
      "PANCREAS_ENRICHED_LOCAL_SHARING_PATTERN",
      fifelse(
        n_non_anchor_H4_ge_080>=1,
        "CROSS_ORGAN_SHARED_LOCAL_PATTERN",
        "NO_CLEAR_HIGH_H4_CROSS_ORGAN_PATTERN"
      )
    )
  ]

  fwrite(
    region_summary,
    file.path(
      reports,
      "Step28_region_specificity_summary.tsv"
    ),
    sep="\t",
    na="NA"
  )

  # ------------------------------------------------------------
  # 16. 与Step26 anchor结果一致性核对（若存在）
  # ------------------------------------------------------------
  step26_files <- list.files(
    file.path(base,"26_cross_trait_colocalization"),
    pattern="Step26_chr19_primary_coloc\\.tsv$|Step26_chr8_positive_control_coloc\\.tsv$",
    recursive=TRUE,
    full.names=TRUE
  )

  anchor_check <- list()

  if (length(step26_files)) {

    f19 <- step26_files[
      grepl("chr19_primary",step26_files)
    ]

    f8 <- step26_files[
      grepl("chr8_positive_control",step26_files)
    ]

    if (length(f19)) {
      old19 <- fread(f19[which.max(file.info(f19)$mtime)])
      oldh4 <- if ("PP.H4" %in% names(old19))
        as.numeric(old19$PP.H4[1]) else NA_real_

      newh4 <- main[
        region_id=="chr19_secondary_shared_locus" &
        pair_id=="pancreas_fat__pancreas_iron",
        PP.H4
      ][1]

      anchor_check[[1]] <- data.table(
        region_id="chr19_secondary_shared_locus",
        prior_Step26_PP_H4=oldh4,
        Step28_PP_H4=newh4,
        absolute_difference=abs(oldh4-newh4)
      )
    }

    if (length(f8)) {
      old8 <- fread(f8[which.max(file.info(f8)$mtime)])
      oldh4 <- if ("PP.H4" %in% names(old8))
        as.numeric(old8$PP.H4[1]) else NA_real_

      newh4 <- main[
        region_id=="chr8_primary_shared_locus" &
        pair_id=="pancreas_fat__pancreas_iron",
        PP.H4
      ][1]

      anchor_check[[2]] <- data.table(
        region_id="chr8_primary_shared_locus",
        prior_Step26_PP_H4=oldh4,
        Step28_PP_H4=newh4,
        absolute_difference=abs(oldh4-newh4)
      )
    }
  }

  if (length(anchor_check)) {

    anchor_check_tab <- rbindlist(
      anchor_check,
      fill=TRUE
    )

    fwrite(
      anchor_check_tab,
      file.path(
        reports,
        "Step28_anchor_reproducibility_check.tsv"
      ),
      sep="\t",
      na="NA"
    )
  }

  # ------------------------------------------------------------
  # 17. QC
  # ------------------------------------------------------------
  summary <- data.table(
    metric=c(
      "n_traits",
      "n_regions",
      "n_pairs",
      "n_main_coloc_expected",
      "n_main_coloc_OK",
      "n_sensitivity_coloc_expected",
      "n_sensitivity_coloc_OK",
      "n_main_H4_ge_050",
      "n_main_H4_ge_080",
      "n_main_H4_ge_090"
    ),
    value=c(
      nrow(trait_info),
      nrow(regions),
      nrow(pairs),
      nrow(regions)*nrow(pairs),
      main[status=="OK",.N],
      nrow(regions)*nrow(pairs)*3L,
      coloc_all[status=="OK",.N],
      main[status=="OK" & PP.H4>=0.50,.N],
      main[status=="OK" & PP.H4>=0.80,.N],
      main[status=="OK" & PP.H4>=0.90,.N]
    )
  )

  fwrite(
    summary,
    file.path(reports,"Step28_summary.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 18. 方法与解释说明
  # ------------------------------------------------------------
  writeLines(
    c(
      "Step28 performs locus-specific cross-organ colocalization at the established chr8 primary shared locus and the exploratory chr19 secondary shared locus.",
      "Six prespecified phenotype pairs are analysed: pancreatic fat-pancreatic iron, liver fat-liver iron, pancreatic fat-liver fat, pancreatic iron-liver iron, pancreatic fat-VAT, and pancreatic fat-ASAT.",
      "All GWAS are MRI-derived traits from the same abdominal imaging GWAS framework and are analysed using raw regional summary statistics without association-P filtering.",
      "Only autosomal biallelic SNVs with valid ES, SE and AF are retained; primary analysis further requires MAF>=0.01 in both traits.",
      "Alleles are harmonized within each pair using same, swapped and strand-complement configurations; the second trait beta is aligned to the first trait ALT allele.",
      "Primary coloc.abf uses quantitative traits, sdY=1, p1=p2=1e-4 and p12=1e-5.",
      "Sensitivity analyses use p12=1e-4, 1e-5 and 1e-6.",
      "PP.H4 supports a shared causal-variant model within the analysed region under the specified coloc assumptions and priors; it does not establish causal direction or a shared effector gene.",
      "High cross-organ PP.H4 is interpreted as a cross-organ shared local genetic pattern, not as proof of the same biological mechanism across tissues.",
      "Low cross-organ PP.H4 does not establish pancreas specificity in an absolute sense because statistical power differs among traits.",
      "The chr19 locus remains exploratory because it was originally prioritized in the Step25b sensitivity analysis.",
      "Step28 complements, but does not replace, the genome-wide LDSC rg results from Steps22-23."
    ),
    file.path(
      reports,
      "methods_and_interpretation.txt"
    )
  )

  # ------------------------------------------------------------
  # 19. 完成判定
  # ------------------------------------------------------------
  review_flags <- character()

  if (main[status=="OK",.N] != nrow(regions)*nrow(pairs)) {
    review_flags <- c(
      review_flags,
      "NOT_ALL_MAIN_COLOC_COMPLETED"
    )
  }

  if (coloc_all[status=="OK",.N] != nrow(regions)*nrow(pairs)*3L) {
    review_flags <- c(
      review_flags,
      "NOT_ALL_SENSITIVITY_COLOC_COMPLETED"
    )
  }

  if (any(harm_qc_tab$n_common_MAF_ge_001 < 100)) {
    review_flags <- c(
      review_flags,
      "AT_LEAST_ONE_PAIR_REGION_HAS_LT100_COMMON_MAF001_VARIANTS"
    )
  }

  if (length(review_flags)) {

    writeLines(
      "COMPLETED_STEP28_WITH_REVIEW_FLAGS",
      file.path(reports,"completion_status.txt")
    )

    writeLines(
      review_flags,
      file.path(reports,"review_flags.txt")
    )

  } else {

    writeLines(
      "COMPLETED_STEP28_CROSS_ORGAN_LOCAL_ARCHITECTURE",
      file.path(reports,"completion_status.txt")
    )
  }

  # ------------------------------------------------------------
  # 20. 精简review包
  # ------------------------------------------------------------
  zipfile <- file.path(
    out,
    "28_cross_organ_local_architecture_review.zip"
  )

  oldwd <- getwd()
  on.exit(setwd(oldwd),add=TRUE)
  setwd(out)

  files_to_zip <- c(
    list.files(
      "reports",
      recursive=TRUE,
      full.names=TRUE
    ),
    list.files(
      "coloc",
      recursive=TRUE,
      full.names=TRUE
    )
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
  # 21. 控制台输出
  # ------------------------------------------------------------
  cat("\n========================================\n")
  cat("第28步完成\n")
  cat("========================================\n\n")

  cat("主分析结果（p12=1e-5）：\n")
  print(
    main[
      ,
      .(
        region_id,
        pair_id,
        comparison_type,
        PP.H3,
        PP.H4,
        evidence_label
      )
    ]
  )

  cat("\n区域特异性/跨器官模式摘要：\n")
  print(region_summary)

  cat("\nStep28 summary：\n")
  print(summary)

  cat("\n结果目录：",out,"\n")
  cat("复核包：",zipfile,"\n")

})()
