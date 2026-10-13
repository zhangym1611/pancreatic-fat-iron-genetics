# ============================================================
# STEP 29b
# Corrected MVMR with GLOBAL UNION LD CLUMPING
#
# Purpose
# -------
# Correct the Step29 MVMR workflow by re-clumping the UNION of
# exposure-specific candidate instruments within each MVMR model.
#
# Why this correction is needed
# -----------------------------
# Exposure-specific clumping guarantees independence only within each
# exposure's IV set. After taking the union across exposures, SNPs from
# different exposure sets may still be in LD. This script therefore:
#
#   1) selects exposure-specific candidate IVs;
#   2) takes the UNION across exposures;
#   3) performs a SECOND / GLOBAL LD clumping on that union using
#      EUR reference, r2 < 0.001, 10,000 kb;
#   4) re-harmonizes all retained SNPs across exposures and outcome;
#   5) re-estimates MVMR by weighted least squares without intercept;
#   6) compares corrected results with the original Step29 MVMR results
#      if those files are available.
#
# Models
# ------
# A) pancreas_fat + VAT -> pancreas_iron
# B) pancreas_fat + ASAT -> pancreas_iron
# C) pancreas_fat + VAT + ASAT -> pancreas_iron
#
# Instrument thresholds
# ---------------------
# pancreas_fat : GWS P < 5e-8
# VAT          : GWS P < 5e-8
# ASAT         : relaxed P < 5e-6 because strict GWS clumping yielded
#                only 2 independent instruments in Step29.
#
# NOTE:
# The ASAT-containing models are therefore sensitivity / exploratory
# MVMR models, not equally strong primary analyses.
#
# LD clumping
# -----------
# r2 < 0.001
# window = 10,000 kb
# global union-level clumping within each MVMR model
#
# Important interpretation boundary
# ---------------------------------
# All traits arise from the same UK Biobank abdominal MRI GWAS framework,
# so substantial sample overlap is expected. MVMR estimates should be
# interpreted as directional multivariable genetic evidence, not as
# independent causal validation.
#
# Output
# ------
# ~/胰腺/29b_MVMR_global_clumping/<run_id>/
#   reports/
#   union_clumping/
#   harmonized/
#   mvmr/
#   logs/
#   29b_MVMR_global_clumping_review.zip
# ============================================================

(function() {

  options(stringsAsFactors = FALSE)
  base <- path.expand("~/胰腺")

  # ------------------------------------------------------------
  # 0. packages
  # ------------------------------------------------------------
  lib <- path.expand("~/R_pancreas_step29b")
  dir.create(lib, recursive=TRUE, showWarnings=FALSE)
  .libPaths(c(lib, .libPaths()))

  need <- function(pkg) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
      install.packages(pkg, lib=lib, repos="https://cloud.r-project.org")
    }
    if (!requireNamespace(pkg, quietly=TRUE)) {
      stop("无法加载R包: ", pkg)
    }
  }

  need("data.table")
  suppressPackageStartupMessages(library(data.table))

  # ------------------------------------------------------------
  # 1. trait definitions
  # ------------------------------------------------------------
  traits <- data.table(
    trait=c("pancreas_fat","pancreas_iron","VAT","ASAT"),
    study_id=c(
      "ebi-a-GCST90016675",
      "ebi-a-GCST90016676",
      "ebi-a-GCST90016671",
      "ebi-a-GCST90016672"
    ),
    N=c(25617,25617,32860,32860)
  )

  traits[, vcf := file.path(base, paste0(study_id, ".vcf.gz"))]

  if (!all(file.exists(traits$vcf))) {
    stop(
      "缺少以下GWAS VCF:\n",
      paste(traits[!file.exists(vcf), vcf], collapse="\n")
    )
  }

  # ------------------------------------------------------------
  # 2. output directories
  # ------------------------------------------------------------
  parent <- file.path(base, "29b_MVMR_global_clumping")
  dir.create(parent, recursive=TRUE, showWarnings=FALSE)

  out <- tempfile("MVMR_global_clump_", tmpdir=parent)
  if (!dir.create(out)) stop("无法创建Step29b输出目录")

  reports <- file.path(out,"reports")
  union_dir <- file.path(out,"union_clumping")
  harm_dir <- file.path(out,"harmonized")
  mvmr_dir <- file.path(out,"mvmr")
  logs <- file.path(out,"logs")

  for (d in c(reports,union_dir,harm_dir,mvmr_dir,logs))
    dir.create(d, recursive=TRUE, showWarnings=FALSE)

  writeLines(
    "INCOMPLETE_CHECK_LOGS",
    file.path(reports,"completion_status.txt")
  )

  fwrite(
    traits[,.(trait,study_id,N)],
    file.path(reports,"Step29b_trait_info.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 3. locate PLINK
  # ------------------------------------------------------------
  plink <- Sys.which("plink")
  plink2 <- Sys.which("plink2")

  if (!nzchar(plink) && !nzchar(plink2)) {
    cand <- list.files(
      c(path.expand("~/mr_tools"), base),
      pattern="^(plink|plink2)$",
      recursive=TRUE,
      full.names=TRUE
    )
    cand <- cand[file.exists(cand)]

    if (length(cand)) {
      if (any(basename(cand)=="plink")) {
        plink <- cand[basename(cand)=="plink"][1]
      } else {
        plink2 <- cand[1]
      }
    }
  }

  if (!nzchar(plink) && !nzchar(plink2)) {
    stop("未找到PLINK/plink2。")
  }

  plink_bin <- if (nzchar(plink)) plink else plink2
  plink_mode <- if (grepl("plink2$",basename(plink_bin))) "plink2" else "plink"

  cat("PLINK: ", plink_bin, "\n", sep="")

  # ------------------------------------------------------------
  # 4. locate EUR LD reference
  # ------------------------------------------------------------
  find_ldref <- function() {

    roots <- c(path.expand("~/mr_tools"), base)

    beds <- unlist(lapply(
      roots,
      function(r) {
        if (!dir.exists(r)) return(character())
        list.files(
          r,
          pattern="\\.bed$",
          recursive=TRUE,
          full.names=TRUE
        )
      }
    ))

    if (!length(beds)) return(NULL)

    stems <- sub("\\.bed$","",beds)

    ok <- file.exists(paste0(stems,".bim")) &
      file.exists(paste0(stems,".fam"))

    stems <- stems[ok]

    if (!length(stems)) return(NULL)

    score <- rep(0,length(stems))
    score <- score + 10*grepl("00_ld_reference",stems,fixed=TRUE)
    score <- score + 5*grepl("EUR|eur",stems)
    score <- score + 4*grepl("1000G|1000g|1KG|1kg",stems)
    score <- score + 2*grepl("503",stems)
    score <- score - 5*grepl(
      "tmp|temp|output|result",
      stems,
      ignore.case=TRUE
    )

    stems[order(score,decreasing=TRUE)][1]
  }

  ldref <- find_ldref()

  if (is.null(ldref) || !nzchar(ldref)) {
    stop(
      "未自动找到EUR PLINK LD reference (.bed/.bim/.fam)。"
    )
  }

  cat("LD reference: ", ldref, "\n", sep="")

  fwrite(
    data.table(
      plink_bin=plink_bin,
      plink_mode=plink_mode,
      ld_reference=ldref,
      clump_r2=0.001,
      clump_kb=10000
    ),
    file.path(reports,"Step29b_environment.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 5. Python parser for GWAS VCF
  # ------------------------------------------------------------
  python_bin <- Sys.which("python3")
  if (!nzchar(python_bin)) python_bin <- Sys.which("python")
  if (!nzchar(python_bin)) stop("未找到Python")

  pyfile <- file.path(out,"extract_gwas_step29b.py")

  py <- c(
'import gzip,sys,math,csv',
'vcf=sys.argv[1]',
'trait=sys.argv[2]',
'outfile=sys.argv[3]',
'',
'def ff(x):',
'    try:',
'        v=float(x)',
'        return v if math.isfinite(v) else None',
'    except:',
'        return None',
'',
'with gzip.open(vcf,"rt") as f, open(outfile,"w",newline="") as o:',
'    w=csv.writer(o,delimiter="\\t")',
'    w.writerow(["trait","CHR","BP","SNP","REF","ALT","BETA","SE","P","AF","MAF"])',
'    for line in f:',
'        if line.startswith("#"):',
'            continue',
'        p=line.rstrip("\\n").split("\\t")',
'        if len(p)<10:',
'            continue',
'        ch=p[0].replace("chr","")',
'        if not ch.isdigit():',
'            continue',
'        ch=int(ch)',
'        if ch<1 or ch>22:',
'            continue',
'        try:',
'            bp=int(p[1])',
'        except:',
'            continue',
'        ref=p[3].upper()',
'        alt=p[4].upper()',
'        if len(ref)!=1 or len(alt)!=1:',
'            continue',
'        if ref not in "ACGT" or alt not in "ACGT":',
'            continue',
'        if p[6] not in ("PASS","."):',
'            continue',
'        keys=p[8].split(":")',
'        vals=p[9].split(":")',
'        if len(keys)!=len(vals):',
'            continue',
'        d=dict(zip(keys,vals))',
'        if not all(k in d for k in ("ES","SE","AF")):',
'            continue',
'        beta=ff(d["ES"])',
'        se=ff(d["SE"])',
'        af=ff(d["AF"])',
'        if beta is None or se is None or af is None:',
'            continue',
'        if se<=0 or not (0<af<1):',
'            continue',
'        pv=None',
'        if "LP" in d:',
'            lp=ff(d["LP"])',
'            if lp is not None:',
'                try:',
'                    pv=10.0**(-lp)',
'                except:',
'                    pv=None',
'        if pv is None:',
'            z=abs(beta/se)',
'            pv=math.erfc(z/math.sqrt(2.0))',
'        maf=min(af,1-af)',
'        rsid=p[2]',
'        if rsid in ("","."):',
'            rsid=f"chr{ch}:{bp}:{ref}:{alt}"',
'        w.writerow([trait,ch,bp,rsid,ref,alt,beta,se,pv,af,maf])'
  )

  writeLines(py, pyfile, useBytes=TRUE)

  gwas <- list()

  for (i in seq_len(nrow(traits))) {

    tr <- traits[i]

    tmp <- file.path(out,paste0(tr$trait,"_full.tsv"))

    logf <- file.path(
      logs,
      paste0(tr$trait,"_extract.log")
    )

    rc <- system2(
      python_bin,
      c(
        shQuote(pyfile),
        shQuote(tr$vcf),
        shQuote(tr$trait),
        shQuote(tmp)
      ),
      stdout=logf,
      stderr=logf
    )

    if (!identical(as.integer(rc),0L)) {
      stop("GWAS提取失败: ", tr$trait)
    }

    x <- fread(tmp)
    unlink(tmp)

    setorder(x,P)
    x <- x[!duplicated(SNP)]

    gwas[[tr$trait]] <- x
  }

  # ------------------------------------------------------------
  # 6. fixed instrument-threshold policy for Step29b
  # ------------------------------------------------------------
  threshold_policy <- data.table(
    trait=c("pancreas_fat","VAT","ASAT"),
    threshold_label=c(
      "GWS_5e8",
      "GWS_5e8",
      "RELAXED_5e6"
    ),
    p_threshold=c(
      5e-8,
      5e-8,
      5e-6
    ),
    role=c(
      "PRIMARY",
      "PRIMARY",
      "SENSITIVITY_BECAUSE_ONLY_2_GWS_IVS"
    )
  )

  fwrite(
    threshold_policy,
    file.path(reports,"Step29b_threshold_policy.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 7. models
  # ------------------------------------------------------------
  models <- list(
    PF_VAT=c("pancreas_fat","VAT"),
    PF_ASAT=c("pancreas_fat","ASAT"),
    PF_VAT_ASAT=c("pancreas_fat","VAT","ASAT")
  )

  # ------------------------------------------------------------
  # 8. helper: allele complement
  # ------------------------------------------------------------
  comp <- function(a) {
    chartr("ACGT","TGCA",toupper(a))
  }

  # ------------------------------------------------------------
  # 9. helper: exposure candidate set
  # ------------------------------------------------------------
  get_candidates <- function(tr) {

    pol <- threshold_policy[trait==tr]

    if (!nrow(pol)) {
      stop("没有阈值策略: ",tr)
    }

    pthr <- pol$p_threshold[1]

    x <- copy(gwas[[tr]][
      P < pthr &
      MAF >= 0.01 &
      grepl("^rs[0-9]+$",SNP)
    ])

    setorder(x,P)
    x <- x[!duplicated(SNP)]

    x[, source_exposure := tr]
    x[, threshold_label := pol$threshold_label[1]]
    x[, p_threshold := pthr]

    x
  }

  # ------------------------------------------------------------
  # 10. GLOBAL union clumping for each model
  # ------------------------------------------------------------
  global_clump <- function(model_name, exposures) {

    cand_list <- lapply(exposures,get_candidates)
    names(cand_list) <- exposures

    # union metadata: one SNP can be candidate for >1 exposure
    long <- rbindlist(cand_list,fill=TRUE)

    source_map <- long[
      ,
      .(
        source_exposures=paste(
          sort(unique(source_exposure)),
          collapse=";"
        ),
        source_thresholds=paste(
          sort(unique(threshold_label)),
          collapse=";"
        )
      ),
      by=SNP
    ]

    # rank SNP by its smallest association P across included exposures
    union <- long[
      ,
      .(
        CLUMP_P=min(P,na.rm=TRUE),
        CHR=CHR[which.min(P)],
        BP=BP[which.min(P)]
      ),
      by=SNP
    ]

    union <- merge(
      union,
      source_map,
      by="SNP",
      all.x=TRUE
    )

    setorder(union,CLUMP_P)

    fwrite(
      union,
      file.path(
        union_dir,
        paste0(model_name,"__union_before_global_clump.tsv")
      ),
      sep="\t"
    )

    clump_input <- file.path(
      union_dir,
      paste0(model_name,"__global_clump_input.tsv")
    )

    fwrite(
      union[,.(SNP,P=CLUMP_P)],
      clump_input,
      sep="\t"
    )

    prefix <- file.path(
      union_dir,
      paste0(model_name,"__global")
    )

    args <- if (plink_mode=="plink") {
      c(
        "--bfile",shQuote(ldref),
        "--clump",shQuote(clump_input),
        "--clump-snp-field","SNP",
        "--clump-field","P",
        "--clump-p1","1",
        "--clump-p2","1",
        "--clump-r2","0.001",
        "--clump-kb","10000",
        "--out",shQuote(prefix)
      )
    } else {
      c(
        "--bfile",shQuote(ldref),
        "--clump",shQuote(clump_input),
        "--clump-id-field","SNP",
        "--clump-p-field","P",
        "--clump-p1","1",
        "--clump-p2","1",
        "--clump-r2","0.001",
        "--clump-kb","10000",
        "--out",shQuote(prefix)
      )
    }

    logf <- file.path(
      logs,
      paste0(model_name,"__global_clump.log")
    )

    rc <- system2(
      plink_bin,
      args=args,
      stdout=logf,
      stderr=logf
    )

    cands <- c(
      paste0(prefix,".clumped"),
      paste0(prefix,".clumps")
    )

    cf <- cands[file.exists(cands)][1]

    kept <- character()

    if (length(cf) && !is.na(cf)) {
      cc <- tryCatch(
        fread(cf,fill=TRUE),
        error=function(e) NULL
      )

      if (!is.null(cc) && nrow(cc)) {
        snpcol <- intersect(
          c("SNP","ID","ID_A"),
          names(cc)
        )

        if (length(snpcol)) {
          kept <- as.character(cc[[snpcol[1]]])
        }
      }
    }

    after <- union[SNP %in% kept]

    # if PLINK returned no SNPs, stop instead of silently continuing
    if (!nrow(after)) {
      stop(
        "GLOBAL clumping没有保留任何SNP: ",
        model_name,
        "。请查看日志: ",
        logf
      )
    }

    fwrite(
      after,
      file.path(
        union_dir,
        paste0(model_name,"__union_after_global_clump.tsv")
      ),
      sep="\t"
    )

    list(
      before=union,
      after=after,
      n_before=nrow(union),
      n_after=nrow(after),
      log=logf
    )
  }

  # ------------------------------------------------------------
  # 11. helper: aligned association extraction
  # ------------------------------------------------------------
  get_assoc_aligned <- function(tr, snps, anchor_dt) {

    x <- copy(gwas[[tr]][SNP %in% snps])

    setorder(x,P)
    x <- x[!duplicated(SNP)]

    keep_cols <- intersect(
      c("SNP","REF","ALT","BETA","SE","AF","P","MAF","CHR","BP"),
      names(x)
    )

    x <- x[, ..keep_cols]

    rename_map <- c(
      REF="REF_x",
      ALT="ALT_x",
      BETA="BETA_x",
      SE="SE_x",
      AF="AF_x",
      P="P_x",
      MAF="MAF_x",
      CHR="CHR_x",
      BP="BP_x"
    )

    present_old <- intersect(
      names(rename_map),
      names(x)
    )

    if (length(present_old)) {
      setnames(
        x,
        old=present_old,
        new=unname(rename_map[present_old])
      )
    }

    m <- merge(
      anchor_dt,
      x,
      by="SNP",
      all.x=TRUE
    )

    required <- c(
      "REF_anchor","ALT_anchor",
      "REF_x","ALT_x",
      "BETA_x","SE_x","AF_x"
    )

    miss <- setdiff(required,names(m))

    if (length(miss)) {
      stop(
        "get_assoc_aligned缺少必要列: ",
        paste(miss,collapse=", "),
        "; trait=",tr
      )
    }

    same <- (
      m$REF_anchor==m$REF_x &
      m$ALT_anchor==m$ALT_x
    )

    swap_raw <- (
      m$REF_anchor==m$ALT_x &
      m$ALT_anchor==m$REF_x
    )

    csame_raw <- (
      comp(m$REF_anchor)==m$REF_x &
      comp(m$ALT_anchor)==m$ALT_x
    )

    cswap_raw <- (
      comp(m$REF_anchor)==m$ALT_x &
      comp(m$ALT_anchor)==m$REF_x
    )

    # mutually exclusive precedence
    swap <- (!same) & swap_raw
    csame <- (!same) & (!swap) & csame_raw
    cswap <- (!same) & (!swap) & (!csame) & cswap_raw

    ok <- (
      !is.na(m$BETA_x) &
      !is.na(m$SE_x) &
      !is.na(m$AF_x) &
      (same | swap | csame | cswap)
    )

    flip <- swap | cswap

    m[, beta_aligned := ifelse(
      flip,
      -BETA_x,
      BETA_x
    )]

    m[, se_aligned := SE_x]

    m[, af_aligned := ifelse(
      flip,
      1-AF_x,
      AF_x
    )]

    m[, align_ok := ok]

    m[, align_class := fifelse(
      same,"same",
      fifelse(
        swap,"swap",
        fifelse(
          csame,"complement_same",
          fifelse(
            cswap,"complement_swap",
            "incompatible"
          )
        )
      )
    )]

    m
  }

  # ------------------------------------------------------------
  # 12. corrected MVMR
  # ------------------------------------------------------------
  run_mvmr_corrected <- function(
    model_name,
    exposures,
    outcome="pancreas_iron"
  ) {

    gc <- global_clump(model_name,exposures)

    snps <- gc$after$SNP

    # anchor allele:
    # among included exposures choose the smallest P representation per SNP
    allrows <- rbindlist(
      lapply(
        exposures,
        function(tr) {
          gwas[[tr]][
            SNP %in% snps,
            .(SNP,REF,ALT,P,source_exposure=tr)
          ]
        }
      ),
      fill=TRUE
    )

    setorder(allrows,P)

    anchor <- allrows[
      !duplicated(SNP),
      .(SNP,REF_anchor=REF,ALT_anchor=ALT)
    ]

    # exposure matrices
    mats <- list()

    for (tr in exposures) {

      a <- get_assoc_aligned(
        tr,
        snps,
        anchor
      )

      mats[[tr]] <- a[
        ,
        .(
          SNP,
          beta_aligned,
          se_aligned,
          af_aligned,
          align_ok,
          align_class
        )
      ]

      setnames(
        mats[[tr]],
        c(
          "beta_aligned",
          "se_aligned",
          "af_aligned",
          "align_ok",
          "align_class"
        ),
        c(
          paste0("bx_",tr),
          paste0("se_",tr),
          paste0("af_",tr),
          paste0("ok_",tr),
          paste0("align_",tr)
        )
      )
    }

    dat <- Reduce(
      function(x,y) {
        merge(x,y,by="SNP",all=TRUE)
      },
      c(list(anchor),mats)
    )

    # outcome
    oy <- get_assoc_aligned(
      outcome,
      snps,
      anchor
    )

    oy <- oy[
      ,
      .(
        SNP,
        beta_aligned,
        se_aligned,
        af_aligned,
        align_ok,
        align_class
      )
    ]

    setnames(
      oy,
      c(
        "beta_aligned",
        "se_aligned",
        "af_aligned",
        "align_ok",
        "align_class"
      ),
      c(
        "by",
        "se_y",
        "af_y",
        "ok_y",
        "align_y"
      )
    )

    dat <- merge(
      dat,
      oy,
      by="SNP",
      all.x=TRUE
    )

    # ----------------------------------------------------------
    # complete-data QC
    # ----------------------------------------------------------
    okcols <- c(
      paste0("ok_",exposures),
      "ok_y"
    )

    betacols <- c(
      paste0("bx_",exposures),
      "by"
    )

    secols <- c(
      paste0("se_",exposures),
      "se_y"
    )

    keep <- rep(TRUE,nrow(dat))

    for (cc in c(okcols,betacols,secols)) {
      keep <- keep & !is.na(dat[[cc]])
    }

    keep <- keep & dat$ok_y

    for (cc in paste0("ok_",exposures)) {
      keep <- keep & dat[[cc]]
    }

    dat <- dat[keep]

    # ----------------------------------------------------------
    # palindromic SNP handling
    # ----------------------------------------------------------
    pal <- (
      (dat$REF_anchor=="A" & dat$ALT_anchor=="T") |
      (dat$REF_anchor=="T" & dat$ALT_anchor=="A") |
      (dat$REF_anchor=="C" & dat$ALT_anchor=="G") |
      (dat$REF_anchor=="G" & dat$ALT_anchor=="C")
    )

    # Use anchor AF from the first exposure with available AF:
    afcols <- paste0("af_",exposures)

    anchor_af <- rep(NA_real_,nrow(dat))

    for (cc in afcols) {
      idx <- is.na(anchor_af) & !is.na(dat[[cc]])
      anchor_af[idx] <- dat[[cc]][idx]
    }

    ambiguous_pal <- (
      pal &
      !is.na(anchor_af) &
      anchor_af>0.42 &
      anchor_af<0.58
    )

    dat[, palindromic := pal]
    dat[, anchor_AF := anchor_af]
    dat[, ambiguous_palindromic := ambiguous_pal]

    # Only ambiguous palindromic SNPs are removed.
    dat <- dat[ambiguous_palindromic == FALSE]

    if (nrow(dat) <= length(exposures)) {
      return(list(
        status="TOO_FEW_MVMR_INSTRUMENTS_AFTER_HARMONIZATION",
        global_clump=gc,
        harmonized=dat
      ))
    }

    # ----------------------------------------------------------
    # weighted multivariable IVW
    # ----------------------------------------------------------
    bx_cols <- paste0("bx_", exposures)

    X <- as.matrix(
      dat[, ..bx_cols]
    )

    y <- dat$by
    w <- 1/(dat$se_y^2)

    XtW <- t(X) %*% (X*w)
    XtWy <- t(X) %*% (y*w)

    beta <- tryCatch(
      as.numeric(solve(XtW,XtWy)),
      error=function(e) {
        rep(NA_real_,ncol(X))
      }
    )

    vc <- tryCatch(
      solve(XtW),
      error=function(e) {
        matrix(
          NA_real_,
          nrow=ncol(X),
          ncol=ncol(X)
        )
      }
    )

    se <- sqrt(diag(vc))
    z <- beta/se
    p <- 2*pnorm(abs(z),lower.tail=FALSE)

    res <- data.table(
      model=model_name,
      exposure=exposures,
      outcome=outcome,
      nsnp=nrow(dat),
      beta=beta,
      se=se,
      p=p
    )

    # attach threshold labels
    res <- merge(
      res,
      threshold_policy[
        ,
        .(
          exposure=trait,
          threshold_label,
          p_threshold,
          threshold_role=role
        )
      ],
      by="exposure",
      all.x=TRUE,
      sort=FALSE
    )

    # ----------------------------------------------------------
    # residual heterogeneity
    # ----------------------------------------------------------
    resid <- y - as.numeric(X %*% beta)

    Q <- sum(w*resid^2)
    Q_df <- nrow(dat)-length(exposures)

    Q_p <- if (Q_df>0) {
      pchisq(
        Q,
        Q_df,
        lower.tail=FALSE
      )
    } else {
      NA_real_
    }

    qtab <- data.table(
      model=model_name,
      nsnp=nrow(dat),
      Q=Q,
      Q_df=Q_df,
      Q_p=Q_p
    )

    # ----------------------------------------------------------
    # descriptive F proxies
    # ----------------------------------------------------------
    ftab <- rbindlist(
      lapply(
        exposures,
        function(tr) {

          z2 <- (
            dat[[paste0("bx_",tr)]] /
            dat[[paste0("se_",tr)]]
          )^2

          data.table(
            model=model_name,
            exposure=tr,
            nsnp=nrow(dat),
            mean_F_proxy=mean(z2,na.rm=TRUE),
            median_F_proxy=median(z2,na.rm=TRUE),
            min_F_proxy=min(z2,na.rm=TRUE),
            n_F_proxy_lt10=sum(z2<10,na.rm=TRUE)
          )
        }
      )
    )

    # ----------------------------------------------------------
    # model collinearity diagnostics
    # ----------------------------------------------------------
    cor_mat <- tryCatch(
      cor(X,use="pairwise.complete.obs"),
      error=function(e) matrix(NA_real_,ncol(X),ncol(X))
    )

    coldiag <- data.table(
      model=model_name,
      metric=c(
        "n_instruments",
        "condition_number_XtWX"
      ),
      value=c(
        nrow(dat),
        tryCatch(
          kappa(XtW),
          error=function(e) NA_real_
        )
      )
    )

    # ----------------------------------------------------------
    # leave-one-out MVMR
    # ----------------------------------------------------------
    loo <- data.table()

    if (nrow(dat) > length(exposures)+1) {

      loo_list <- vector("list",nrow(dat))

      for (k in seq_len(nrow(dat))) {

        Xk <- X[-k,,drop=FALSE]
        yk <- y[-k]
        wk <- w[-k]

        XtWk <- t(Xk) %*% (Xk*wk)
        XtWyk <- t(Xk) %*% (yk*wk)

        bk <- tryCatch(
          as.numeric(solve(XtWk,XtWyk)),
          error=function(e) rep(NA_real_,ncol(Xk))
        )

        loo_list[[k]] <- data.table(
          model=model_name,
          omitted_SNP=dat$SNP[k],
          exposure=exposures,
          beta=bk
        )
      }

      loo <- rbindlist(loo_list,fill=TRUE)
    }

    list(
      status="OK",
      global_clump=gc,
      harmonized=dat,
      result=res,
      heterogeneity=qtab,
      F_proxy=ftab,
      cor_matrix=cor_mat,
      collinearity=coldiag,
      leave_one_out=loo
    )
  }

  # ------------------------------------------------------------
  # 13. run all models
  # ------------------------------------------------------------
  results <- list()
  hets <- list()
  fproxies <- list()
  coldiags <- list()
  loos <- list()
  qc <- list()

  for (nm in names(models)) {

    exposures <- models[[nm]]

    cat(
      "\nRunning corrected MVMR model: ",
      nm,
      " [",
      paste(exposures,collapse=" + "),
      "]\n",
      sep=""
    )

    rr <- run_mvmr_corrected(
      nm,
      exposures,
      outcome="pancreas_iron"
    )

    qc[[nm]] <- data.table(
      model=nm,
      exposures=paste(exposures,collapse="+"),
      status=rr$status,
      n_union_before=if (!is.null(rr$global_clump))
        rr$global_clump$n_before else NA_integer_,
      n_union_after_global_clump=if (!is.null(rr$global_clump))
        rr$global_clump$n_after else NA_integer_,
      n_after_harmonization=if (!is.null(rr$harmonized))
        nrow(rr$harmonized) else NA_integer_
    )

    if (!is.null(rr$harmonized)) {
      fwrite(
        rr$harmonized,
        file.path(
          harm_dir,
          paste0(
            "Step29b_",
            nm,
            "__harmonized_after_global_clump.tsv"
          )
        ),
        sep="\t"
      )
    }

    if (identical(rr$status,"OK")) {

      results[[nm]] <- rr$result
      hets[[nm]] <- rr$heterogeneity
      fproxies[[nm]] <- rr$F_proxy
      coldiags[[nm]] <- rr$collinearity

      if (nrow(rr$leave_one_out)) {
        loos[[nm]] <- rr$leave_one_out
      }

      # save exposure correlation matrix
      cm <- as.data.table(
        rr$cor_matrix,
        keep.rownames="exposure1"
      )

      setnames(
        cm,
        old=setdiff(names(cm),"exposure1"),
        new=models[[nm]]
      )

      fwrite(
        cm,
        file.path(
          mvmr_dir,
          paste0(
            "Step29b_",
            nm,
            "__exposure_beta_correlation.tsv"
          )
        ),
        sep="\t"
      )
    }
  }

  qc_tab <- rbindlist(qc,fill=TRUE)

  fwrite(
    qc_tab,
    file.path(reports,"Step29b_MVMR_QC.tsv"),
    sep="\t"
  )

  if (length(results)) {
    res_tab <- rbindlist(results,fill=TRUE)

    fwrite(
      res_tab,
      file.path(
        mvmr_dir,
        "Step29b_corrected_MVMR_results.tsv"
      ),
      sep="\t"
    )
  } else {
    res_tab <- data.table()
  }

  if (length(hets)) {
    fwrite(
      rbindlist(hets,fill=TRUE),
      file.path(
        mvmr_dir,
        "Step29b_corrected_MVMR_heterogeneity.tsv"
      ),
      sep="\t"
    )
  }

  if (length(fproxies)) {
    fwrite(
      rbindlist(fproxies,fill=TRUE),
      file.path(
        mvmr_dir,
        "Step29b_corrected_MVMR_F_proxies.tsv"
      ),
      sep="\t"
    )
  }

  if (length(coldiags)) {
    fwrite(
      rbindlist(coldiags,fill=TRUE),
      file.path(
        mvmr_dir,
        "Step29b_corrected_MVMR_collinearity.tsv"
      ),
      sep="\t"
    )
  }

  if (length(loos)) {
    fwrite(
      rbindlist(loos,fill=TRUE),
      file.path(
        mvmr_dir,
        "Step29b_corrected_MVMR_leave_one_out.tsv"
      ),
      sep="\t"
    )
  }

  # ------------------------------------------------------------
  # 14. compare with original Step29 MVMR if available
  # ------------------------------------------------------------
  old_files <- list.files(
    file.path(base,"29_directional_MR_MVMR"),
    pattern="Step29_MVMR_results\\.tsv$",
    recursive=TRUE,
    full.names=TRUE
  )

  comparison <- data.table()

  if (length(old_files) && nrow(res_tab)) {

    oldf <- old_files[
      which.max(file.info(old_files)$mtime)
    ]

    old <- fread(oldf)

    keep_old <- intersect(
      c(
        "model","exposure","outcome","nsnp",
        "beta","se","p","threshold_used"
      ),
      names(old)
    )

    old <- old[, ..keep_old]

    setnames(
      old,
      intersect(c("nsnp","beta","se","p"),names(old)),
      paste0(
        intersect(c("nsnp","beta","se","p"),names(old)),
        "_Step29_original"
      )
    )

    new <- copy(res_tab)

    setnames(
      new,
      c("nsnp","beta","se","p"),
      c(
        "nsnp_Step29b_corrected",
        "beta_Step29b_corrected",
        "se_Step29b_corrected",
        "p_Step29b_corrected"
      )
    )

    comparison <- merge(
      old,
      new,
      by=c("model","exposure","outcome"),
      all=TRUE
    )

    if (
      all(
        c(
          "beta_Step29_original",
          "beta_Step29b_corrected"
        ) %in% names(comparison)
      )
    ) {
      comparison[
        ,
        beta_change :=
          beta_Step29b_corrected -
          beta_Step29_original
      ]
    }

    fwrite(
      comparison,
      file.path(
        reports,
        "Step29_vs_Step29b_MVMR_comparison.tsv"
      ),
      sep="\t"
    )
  }

  # ------------------------------------------------------------
  # 15. summary
  # ------------------------------------------------------------
  summary <- data.table(
    metric=c(
      "n_models",
      "n_models_OK",
      "total_union_before_clump",
      "total_union_after_global_clump",
      "total_after_harmonization"
    ),
    value=c(
      length(models),
      sum(qc_tab$status=="OK"),
      sum(qc_tab$n_union_before,na.rm=TRUE),
      sum(qc_tab$n_union_after_global_clump,na.rm=TRUE),
      sum(qc_tab$n_after_harmonization,na.rm=TRUE)
    )
  )

  fwrite(
    summary,
    file.path(reports,"Step29b_summary.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 16. methods / interpretation
  # ------------------------------------------------------------
  writeLines(
    c(
      "Step29b corrects the original multivariable MR workflow by applying an additional global LD clumping step to the union of exposure-specific candidate instruments within each MVMR model.",
      "This is necessary because SNP sets clumped independently for different exposures can still be correlated with one another after union.",
      "Global clumping uses the European PLINK reference with r2<0.001 and a 10,000-kb window.",
      "Pancreatic fat and VAT use genome-wide significant candidate instruments (P<5e-8).",
      "ASAT uses the relaxed P<5e-6 candidate threshold because only two independent genome-wide significant ASAT instruments were available in Step29; therefore ASAT-containing models are sensitivity/exploratory analyses.",
      "The global clumping ranking P value for each SNP is the smallest association P across the exposures included in that MVMR model.",
      "After global clumping, all retained variants are harmonized to a common effect allele across exposures and pancreatic iron.",
      "Palindromic SNPs are retained only when not frequency-ambiguous; ambiguous palindromic SNPs with anchor AF between 0.42 and 0.58 are excluded.",
      "MVMR is estimated using weighted least squares without an intercept, weighted by the inverse variance of the pancreatic-iron association.",
      "Residual Cochran-type Q, leave-one-out estimates, exposure-association correlation matrices and matrix condition numbers are reported as model diagnostics.",
      "The reported F statistics are descriptive per-exposure association-strength proxies and are NOT exact Sanderson-Windmeijer conditional F statistics because covariance information for the exposure association estimates is unavailable.",
      "All traits originate from the same UK Biobank abdominal MRI GWAS framework, so substantial sample overlap is expected. Results should therefore be interpreted as directional multivariable genetic evidence rather than independent causal validation."
    ),
    file.path(
      reports,
      "methods_and_interpretation.txt"
    )
  )

  # ------------------------------------------------------------
  # 17. completion flags
  # ------------------------------------------------------------
  flags <- character()

  if (any(qc_tab$status!="OK")) {
    flags <- c(
      flags,
      "AT_LEAST_ONE_MVMR_MODEL_NOT_ESTIMABLE"
    )
  }

  if (
    any(
      qc_tab$n_union_after_global_clump <
      qc_tab$n_union_before,
      na.rm=TRUE
    )
  ) {
    flags <- c(
      flags,
      "GLOBAL_CLUMPING_REMOVED_CORRELATED_UNION_INSTRUMENTS_EXPECTED"
    )
  }

  if (
    length(fproxies)
  ) {
    ft <- rbindlist(fproxies,fill=TRUE)

    if (any(ft$mean_F_proxy < 10,na.rm=TRUE)) {
      flags <- c(
        flags,
        "AT_LEAST_ONE_EXPOSURE_MEAN_F_PROXY_LT10"
      )
    }
  }

  if (length(flags)) {

    writeLines(
      "COMPLETED_STEP29B_WITH_REVIEW_FLAGS",
      file.path(
        reports,
        "completion_status.txt"
      )
    )

    writeLines(
      flags,
      file.path(
        reports,
        "review_flags.txt"
      )
    )

  } else {

    writeLines(
      "COMPLETED_STEP29B_GLOBAL_UNION_CLUMPED_MVMR",
      file.path(
        reports,
        "completion_status.txt"
      )
    )
  }

  # ------------------------------------------------------------
  # 18. review zip
  # ------------------------------------------------------------
  zipfile <- file.path(
    out,
    "29b_MVMR_global_clumping_review.zip"
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
      "union_clumping",
      pattern="\\.tsv$",
      recursive=TRUE,
      full.names=TRUE
    ),
    list.files(
      "harmonized",
      recursive=TRUE,
      full.names=TRUE
    ),
    list.files(
      "mvmr",
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
  # 19. console
  # ------------------------------------------------------------
  cat("\n========================================\n")
  cat("STEP29b completed\n")
  cat("========================================\n\n")

  cat("Global-clumping QC:\n")
  print(qc_tab)

  if (nrow(res_tab)) {
    cat("\nCorrected MVMR results:\n")
    print(res_tab)
  }

  if (nrow(comparison)) {
    cat("\nStep29 vs Step29b comparison:\n")
    print(comparison)
  }

  cat("\nSummary:\n")
  print(summary)

  cat(
    "\nOutput directory:\n",
    out,
    "\n",
    sep=""
  )

  cat(
    "Review zip:\n",
    zipfile,
    "\n",
    sep=""
  )

})()
