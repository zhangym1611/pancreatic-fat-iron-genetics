# ============================================================
# STEP 29
# 修正版v6：修复关键的palindromic allele harmonization优先级错误。A/T或C/G SNP若REF/ALT已完全一致，必须优先判定same而不能同时因complement_swap被反向；UVMR和MVMR均已修复。
# Directional causal integration: bidirectional MR + MVMR
# Pancreatic fat, pancreatic iron, VAT and ASAT
#
# Primary scientific questions:
# 1) Is pancreatic fat genetically associated with pancreatic iron in a
#    directional MR framework?
# 2) Is the reverse direction supported?
# 3) Does the pancreatic-fat -> pancreatic-iron signal persist after
#    accounting for visceral adiposity (VAT)?
# 4) Is the result specific to VAT rather than general adiposity
#    (ASAT comparator)?
#
# Important interpretation boundary:
# - All abdominal MRI traits originate from the same UK Biobank imaging
#   GWAS framework, so substantial sample overlap is expected.
# - Therefore this is NOT an independent two-sample MR validation.
# - MR estimates should be interpreted as directional genetic evidence
#   under IV assumptions, not definitive causality.
#
# GWAS (GRCh37):
# pancreas_fat  GCST90016675, N=25,617
# pancreas_iron GCST90016676, N=25,617
# VAT           GCST90016671, N=32,860
# ASAT          GCST90016672, N=32,860
#
# Instrument strategy:
# - Main threshold: P < 5e-8
# - Sensitivity threshold: P < 5e-6 (only if needed)
# - LD clumping: r2 < 0.001, 10,000 kb
# - Prefer existing EUR PLINK reference on server; script auto-discovers
# - Palindromic variants with ambiguous frequency are excluded
#
# UVMR:
# - IVW
# - weighted median (>=3 IVs)
# - MR-Egger (>=3 IVs)
# - Cochran Q
# - Egger intercept
# - leave-one-out IVW
#
# MVMR:
# A) pancreas_fat + VAT -> pancreas_iron
# B) pancreas_fat + ASAT -> pancreas_iron
# C) pancreas_fat + VAT + ASAT -> pancreas_iron
#
# MVMR estimator:
# - multivariable IVW using weighted least squares without intercept
# - instruments = union of exposure-specific clumped GWS instruments
# - exposure associations harmonised to one common effect allele
# - weak-instrument diagnostics: conditional-style F proxy from
#   exposure association strength (descriptive only; exact Sanderson-
#   Windmeijer F requires covariance information not available here)
#
# Outputs:
# ~/胰腺/29_directional_MR_MVMR/<run_id>/
#   reports/
#   instruments/
#   harmonized/
#   uvmr/
#   mvmr/
#   logs/
#   29_MR_MVMR_review.zip
# ============================================================

(function() {

  options(stringsAsFactors = FALSE)
  base <- path.expand("~/胰腺")

  # ------------------------------------------------------------
  # 0. Packages
  # ------------------------------------------------------------
  lib <- path.expand("~/R_pancreas_step29")
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

  for (p in c("data.table")) need(p)
  suppressPackageStartupMessages(library(data.table))

  # ------------------------------------------------------------
  # 1. Trait definitions
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
  # 2. Output dirs
  # ------------------------------------------------------------
  parent <- file.path(base, "29_directional_MR_MVMR")
  dir.create(parent, recursive=TRUE, showWarnings=FALSE)

  out <- tempfile("MR_MVMR_", tmpdir=parent)
  if (!dir.create(out)) stop("无法创建Step29输出目录")

  reports <- file.path(out,"reports")
  instruments_dir <- file.path(out,"instruments")
  harmonized_dir <- file.path(out,"harmonized")
  uvmr_dir <- file.path(out,"uvmr")
  mvmr_dir <- file.path(out,"mvmr")
  logs <- file.path(out,"logs")

  for (d in c(reports,instruments_dir,harmonized_dir,uvmr_dir,mvmr_dir,logs))
    dir.create(d, recursive=TRUE, showWarnings=FALSE)

  writeLines("INCOMPLETE_CHECK_LOGS",
             file.path(reports,"completion_status.txt"))

  fwrite(traits[,.(trait,study_id,N)],
         file.path(reports,"Step29_trait_info.tsv"), sep="\t")

  # ------------------------------------------------------------
  # 3. Find PLINK
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
    stop(
      "未找到PLINK/plink2。Step29需要LD clumping。\n",
      "请先确认此前使用的EUR LD参考和PLINK路径。"
    )
  }

  plink_bin <- if (nzchar(plink)) plink else plink2
  plink_mode <- if (grepl("plink2$",basename(plink_bin))) "plink2" else "plink"

  cat("PLINK: ",plink_bin,"\n",sep="")

  # ------------------------------------------------------------
  # 4. Discover EUR LD reference
  # ------------------------------------------------------------
  find_ldref <- function() {

    roots <- c(
      path.expand("~/mr_tools"),
      base
    )

    beds <- unlist(lapply(
      roots,
      function(r) {
        if (!dir.exists(r)) return(character())
        list.files(r, pattern="\\.bed$", recursive=TRUE, full.names=TRUE)
      }
    ))

    if (!length(beds)) return(NULL)

    stems <- sub("\\.bed$","",beds)
    ok <- file.exists(paste0(stems,".bim")) & file.exists(paste0(stems,".fam"))
    stems <- stems[ok]

    if (!length(stems)) return(NULL)

    # Prefer EUR / 1KG / 1000G-like refs; deprioritize GWAS temp outputs
    score <- rep(0,length(stems))
    score <- score + 5*grepl("EUR|eur",stems)
    score <- score + 4*grepl("1000G|1kg|1KG|1000g",stems)
    score <- score + 3*grepl("503",stems)
    score <- score - 5*grepl("tmp|temp|output|result",stems,ignore.case=TRUE)

    stems[order(score,decreasing=TRUE)][1]
  }

  ldref <- find_ldref()

  if (is.null(ldref) || !nzchar(ldref)) {
    stop(
      "未自动找到可用于clumping的EUR PLINK参考(.bed/.bim/.fam)。\n",
      "请把此前503 EUR参考放在~/mr_tools或~/胰腺目录下。"
    )
  }

  cat("LD reference: ",ldref,"\n",sep="")

  fwrite(
    data.table(
      plink_bin=plink_bin,
      plink_mode=plink_mode,
      ld_reference=ldref
    ),
    file.path(reports,"Step29_environment.tsv"),
    sep="\t"
  )

  # ------------------------------------------------------------
  # 5. Python standard-library parser for GWAS VCF
  # ------------------------------------------------------------
  python_bin <- Sys.which("python3")
  if (!nzchar(python_bin)) python_bin <- Sys.which("python")
  if (!nzchar(python_bin)) stop("未找到Python")

  pyfile <- file.path(out,"extract_gwas_for_mr.py")

  py <- c(
'import gzip,sys,math,csv,os',
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
'        ref=p[3].upper(); alt=p[4].upper()',
'        if len(ref)!=1 or len(alt)!=1 or ref not in "ACGT" or alt not in "ACGT":',
'            continue',
'        if p[6] not in ("PASS","."):',
'            continue',
'        keys=p[8].split(":"); vals=p[9].split(":")',
'        if len(keys)!=len(vals):',
'            continue',
'        d=dict(zip(keys,vals))',
'        if not all(k in d for k in ("ES","SE","AF")):',
'            continue',
'        beta=ff(d["ES"]); se=ff(d["SE"]); af=ff(d["AF"])',
'        if beta is None or se is None or af is None or se<=0 or not (0<af<1):',
'            continue',
'        # Prefer LP if present; LP=-log10(P)',
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

  # extract each trait
  gwas <- list()

  for (i in seq_len(nrow(traits))) {
    tr <- traits[i]
    outf <- file.path(out, paste0(tr$trait,"_full.tsv.gz"))
    tmp <- sub("\\.gz$","",outf)

    rc <- system2(
      python_bin,
      c(shQuote(pyfile),shQuote(tr$vcf),shQuote(tr$trait),shQuote(tmp)),
      stdout=file.path(logs,paste0(tr$trait,"_extract.log")),
      stderr=file.path(logs,paste0(tr$trait,"_extract.log"))
    )

    if (!identical(as.integer(rc),0L)) {
      stop("GWAS提取失败: ",tr$trait)
    }

    x <- fread(tmp)
    fwrite(x,outf,sep="\t")
    unlink(tmp)

    gwas[[tr$trait]] <- x
  }

  # ------------------------------------------------------------
  # 6. Build GWS candidates and clump
  # ------------------------------------------------------------
  thresholds <- data.table(
    threshold_label=c("GWS_5e8","RELAXED_5e6"),
    p_threshold=c(5e-8,5e-6)
  )

  clump_one <- function(x, trait, pthr, label) {

    cand <- x[
      P < pthr &
      MAF >= 0.01 &
      grepl("^rs[0-9]+$",SNP)
    ]

    if (!nrow(cand)) {
      return(list(
        candidates=cand,
        clumped=cand,
        status="NO_CANDIDATES"
      ))
    }

    # remove duplicate rsIDs, keep smallest P
    setorder(cand,P)
    cand <- cand[!duplicated(SNP)]

    infile <- file.path(
      instruments_dir,
      paste0(trait,"__",label,"__clump_input.tsv")
    )

    fwrite(
      cand[,.(SNP,P)],
      infile,
      sep="\t"
    )

    prefix <- file.path(
      instruments_dir,
      paste0(trait,"__",label)
    )

    args <- if (plink_mode=="plink") {
      c(
        "--bfile",shQuote(ldref),
        "--clump",shQuote(infile),
        "--clump-snp-field","SNP",
        "--clump-field","P",
        "--clump-p1",format(pthr,scientific=TRUE),
        "--clump-p2","1",
        "--clump-r2","0.001",
        "--clump-kb","10000",
        "--out",shQuote(prefix)
      )
    } else {
      c(
        "--bfile",shQuote(ldref),
        "--clump",shQuote(infile),
        "--clump-id-field","SNP",
        "--clump-p-field","P",
        "--clump-p1",format(pthr,scientific=TRUE),
        "--clump-p2","1",
        "--clump-r2","0.001",
        "--clump-kb","10000",
        "--out",shQuote(prefix)
      )
    }

    rc <- system2(
      plink_bin,
      args=args,
      stdout=file.path(logs,paste0(trait,"__",label,"__plink.log")),
      stderr=file.path(logs,paste0(trait,"__",label,"__plink.log"))
    )

    # Parse output robustly
    clumped_rs <- character()

    cands <- c(
      paste0(prefix,".clumped"),
      paste0(prefix,".clumps")
    )
    cf <- cands[file.exists(cands)][1]

    if (length(cf) && !is.na(cf)) {
      cc <- tryCatch(fread(cf,fill=TRUE),error=function(e) NULL)
      if (!is.null(cc) && nrow(cc)) {
        snpcol <- intersect(c("SNP","ID","ID_A"),names(cc))
        if (length(snpcol)) clumped_rs <- as.character(cc[[snpcol[1]]])
      }
    }

    cl <- cand[SNP %in% clumped_rs]

    if (!nrow(cl) && length(clumped_rs)==0) {
      status <- "CLUMP_OUTPUT_EMPTY"
    } else {
      status <- "OK"
    }

    list(
      candidates=cand,
      clumped=cl,
      status=status
    )
  }

  clump_results <- list()
  clump_qc <- list()

  for (tr in traits$trait) {
    for (j in seq_len(nrow(thresholds))) {

      lab <- thresholds$threshold_label[j]
      pthr <- thresholds$p_threshold[j]

      rr <- clump_one(gwas[[tr]],tr,pthr,lab)
      key <- paste(tr,lab,sep="__")
      clump_results[[key]] <- rr

      fwrite(
        rr$candidates,
        file.path(instruments_dir,paste0(key,"__candidates.tsv")),
        sep="\t"
      )

      fwrite(
        rr$clumped,
        file.path(instruments_dir,paste0(key,"__clumped.tsv")),
        sep="\t"
      )

      clump_qc[[key]] <- data.table(
        trait=tr,
        threshold_label=lab,
        p_threshold=pthr,
        n_candidates=nrow(rr$candidates),
        n_clumped=nrow(rr$clumped),
        status=rr$status
      )
    }
  }

  clump_qc_tab <- rbindlist(clump_qc)
  fwrite(clump_qc_tab,
         file.path(reports,"Step29_clumping_QC.tsv"),sep="\t")

  # ------------------------------------------------------------
  # 7. Harmonization utilities
  # ------------------------------------------------------------
  comp <- function(a) chartr("ACGT","TGCA",toupper(a))

  harmonize_exposure_outcome <- function(exp, out, expname, outname) {

    # merge by rsID first
    a <- copy(exp)
    b <- copy(out)

    setorder(a,P); a <- a[!duplicated(SNP)]
    setorder(b,P); b <- b[!duplicated(SNP)]

    m <- merge(
      a,b,
      by="SNP",
      suffixes=c("_exp","_out"),
      all=FALSE
    )

    if (!nrow(m)) return(m)

    same <- m$REF_exp==m$REF_out & m$ALT_exp==m$ALT_out
    swap_raw <- m$REF_exp==m$ALT_out & m$ALT_exp==m$REF_out
    csame_raw <- comp(m$REF_exp)==m$REF_out & comp(m$ALT_exp)==m$ALT_out
    cswap_raw <- comp(m$REF_exp)==m$ALT_out & comp(m$ALT_exp)==m$REF_out

    # IMPORTANT: palindromic SNPs can satisfy more than one relation at once.
    # Use mutually exclusive precedence:
    # same > swap > complement_same > complement_swap.
    # This prevents A/T or C/G SNPs with identical allele coding from being
    # erroneously flipped merely because they also satisfy complement_swap.
    swap <- (!same) & swap_raw
    csame <- (!same) & (!swap) & csame_raw
    cswap <- (!same) & (!swap) & (!csame) & cswap_raw

    ok <- same | swap | csame | cswap
    m <- m[ok]

    if (!nrow(m)) return(m)

    same <- same[ok]; swap <- swap[ok]; csame <- csame[ok]; cswap <- cswap[ok]
    flip <- swap | cswap

    # 显式使用m$列名，避免data.table作用域歧义
    m[, beta_out_aligned := ifelse(flip,-m$BETA_out,m$BETA_out)]
    m[, af_out_aligned := ifelse(flip,1-m$AF_out,m$AF_out)]

    # ambiguous palindromic check
    pal <- (
      (m$REF_exp=="A" & m$ALT_exp=="T") |
      (m$REF_exp=="T" & m$ALT_exp=="A") |
      (m$REF_exp=="C" & m$ALT_exp=="G") |
      (m$REF_exp=="G" & m$ALT_exp=="C")
    )

    ambig <- pal &
      m$AF_exp>0.42 & m$AF_exp<0.58

    m[, palindromic := pal]
    m[, ambiguous_palindromic := ambig]
    # data.table新版不建议DT[!col]这种单符号写法；显式按逻辑列筛选
    m <- m[ambiguous_palindromic == FALSE]

    m[, exposure := expname]
    m[, outcome := outname]
    m
  }

  # ------------------------------------------------------------
  # 8. UVMR estimators
  # ------------------------------------------------------------
  ivw_fit <- function(bx,by,se_y) {
    w <- 1/(se_y^2)
    beta <- sum(w*bx*by)/sum(w*bx^2)
    se <- sqrt(1/sum(w*bx^2))
    z <- beta/se
    p <- 2*pnorm(abs(z),lower.tail=FALSE)
    list(beta=beta,se=se,p=p)
  }

  q_ivw <- function(bx,by,se_y,beta_ivw) {
    w <- 1/(se_y^2)
    q <- sum(w*(by-beta_ivw*bx)^2)
    df <- length(bx)-1
    p <- if (df>0) pchisq(q,df,lower.tail=FALSE) else NA_real_
    c(Q=q,Q_df=df,Q_p=p)
  }

  egger_fit <- function(bx,by,se_y) {
    w <- 1/(se_y^2)
    fit <- lm(by ~ bx, weights=w)
    co <- summary(fit)$coef
    list(
      intercept=co["(Intercept)","Estimate"],
      intercept_se=co["(Intercept)","Std. Error"],
      intercept_p=co["(Intercept)","Pr(>|t|)"],
      beta=co["bx","Estimate"],
      se=co["bx","Std. Error"],
      p=co["bx","Pr(>|t|)"]
    )
  }

  weighted_median <- function(beta_ratio,se_ratio) {
    keep <- is.finite(beta_ratio) & is.finite(se_ratio) & se_ratio>0
    beta_ratio <- beta_ratio[keep]
    se_ratio <- se_ratio[keep]
    if (!length(beta_ratio)) return(NA_real_)

    w <- 1/(se_ratio^2)
    ord <- order(beta_ratio)
    b <- beta_ratio[ord]
    w <- w[ord]/sum(w)
    cw <- cumsum(w)
    b[which(cw>=0.5)[1]]
  }

  run_uvmr <- function(exp_trait,out_trait,threshold_label) {

    inst <- clump_results[[paste(exp_trait,threshold_label,sep="__")]]$clumped
    outg <- gwas[[out_trait]]

    h <- harmonize_exposure_outcome(inst,outg,exp_trait,out_trait)

    outfile <- file.path(
      harmonized_dir,
      paste0(exp_trait,"__",out_trait,"__",threshold_label,"__harmonized.tsv")
    )
    fwrite(h,outfile,sep="\t")

    if (!nrow(h)) {
      return(list(
        results=data.table(
          exposure=exp_trait,outcome=out_trait,
          threshold_label=threshold_label,
          method="NO_DATA",nsnp=0,
          beta=NA_real_,se=NA_real_,p=NA_real_
        ),
        het=data.table(),
        egger=data.table(),
        loo=data.table(),
        harmonized=h
      ))
    }

    bx <- h$BETA_exp
    by <- h$beta_out_aligned
    sy <- h$SE_out

    outrows <- list()

    if (length(bx)==1) {
      beta <- by/bx
      se <- sy/abs(bx)
      p <- 2*pnorm(abs(beta/se),lower.tail=FALSE)
      outrows[[1]] <- data.table(
        exposure=exp_trait,outcome=out_trait,
        threshold_label=threshold_label,
        method="Wald_ratio",nsnp=1,
        beta=beta,se=se,p=p
      )
    } else {
      ivw <- ivw_fit(bx,by,sy)
      outrows[[1]] <- data.table(
        exposure=exp_trait,outcome=out_trait,
        threshold_label=threshold_label,
        method="IVW",nsnp=length(bx),
        beta=ivw$beta,se=ivw$se,p=ivw$p
      )

      # weighted median
      ratio <- by/bx
      se_ratio <- sy/abs(bx)
      wm <- weighted_median(ratio,se_ratio)

      # approximate bootstrap-free SE using weighted MAD / sqrt(n)
      wm_se <- if (length(ratio)>=3) {
        1.4826*median(abs(ratio-wm),na.rm=TRUE)/sqrt(length(ratio))
      } else NA_real_

      wm_p <- if (is.finite(wm_se) && wm_se>0)
        2*pnorm(abs(wm/wm_se),lower.tail=FALSE) else NA_real_

      outrows[[2]] <- data.table(
        exposure=exp_trait,outcome=out_trait,
        threshold_label=threshold_label,
        method="Weighted_median_approx",
        nsnp=length(bx),
        beta=wm,se=wm_se,p=wm_p
      )
    }

    het <- data.table()
    eg <- data.table()
    loo <- data.table()

    if (length(bx)>=2) {
      ivw <- ivw_fit(bx,by,sy)
      q <- q_ivw(bx,by,sy,ivw$beta)
      het <- data.table(
        exposure=exp_trait,outcome=out_trait,
        threshold_label=threshold_label,
        nsnp=length(bx),
        Q=q["Q"],Q_df=q["Q_df"],Q_p=q["Q_p"]
      )

      if (length(bx)>=3) {
        e <- egger_fit(bx,by,sy)
        eg <- data.table(
          exposure=exp_trait,outcome=out_trait,
          threshold_label=threshold_label,
          nsnp=length(bx),
          egger_beta=e$beta,
          egger_se=e$se,
          egger_p=e$p,
          intercept=e$intercept,
          intercept_se=e$intercept_se,
          intercept_p=e$intercept_p
        )

        outrows[[length(outrows)+1]] <- data.table(
          exposure=exp_trait,outcome=out_trait,
          threshold_label=threshold_label,
          method="MR_Egger",
          nsnp=length(bx),
          beta=e$beta,se=e$se,p=e$p
        )
      }

      if (length(bx)>=3) {
        loo <- rbindlist(lapply(seq_along(bx),function(k) {
          fit <- ivw_fit(bx[-k],by[-k],sy[-k])
          data.table(
            exposure=exp_trait,outcome=out_trait,
            threshold_label=threshold_label,
            omitted_SNP=h$SNP[k],
            nsnp=length(bx)-1,
            beta=fit$beta,se=fit$se,p=fit$p
          )
        }))
      }
    }

    list(
      results=rbindlist(outrows,fill=TRUE),
      het=het,
      egger=eg,
      loo=loo,
      harmonized=h
    )
  }

  # Primary bidirectional MR
  uvmr_specs <- data.table(
    exposure=c("pancreas_fat","pancreas_iron"),
    outcome=c("pancreas_iron","pancreas_fat")
  )

  uvmr_all <- list()
  uvmr_het <- list()
  uvmr_egger <- list()
  uvmr_loo <- list()

  for (i in seq_len(nrow(uvmr_specs))) {
    for (lab in thresholds$threshold_label) {

      rr <- run_uvmr(
        uvmr_specs$exposure[i],
        uvmr_specs$outcome[i],
        lab
      )

      key <- paste(
        uvmr_specs$exposure[i],
        uvmr_specs$outcome[i],
        lab,
        sep="__"
      )

      uvmr_all[[key]] <- rr$results
      if (nrow(rr$het)) uvmr_het[[key]] <- rr$het
      if (nrow(rr$egger)) uvmr_egger[[key]] <- rr$egger
      if (nrow(rr$loo)) uvmr_loo[[key]] <- rr$loo
    }
  }

  uvmr_res <- rbindlist(uvmr_all,fill=TRUE)
  fwrite(uvmr_res,file.path(uvmr_dir,"Step29_UVMR_results.tsv"),sep="\t")

  if (length(uvmr_het))
    fwrite(rbindlist(uvmr_het,fill=TRUE),
           file.path(uvmr_dir,"Step29_UVMR_heterogeneity.tsv"),sep="\t")

  if (length(uvmr_egger))
    fwrite(rbindlist(uvmr_egger,fill=TRUE),
           file.path(uvmr_dir,"Step29_UVMR_Egger_intercept.tsv"),sep="\t")

  if (length(uvmr_loo))
    fwrite(rbindlist(uvmr_loo,fill=TRUE),
           file.path(uvmr_dir,"Step29_UVMR_leave_one_out.tsv"),sep="\t")

  # ------------------------------------------------------------
  # 9. MVMR harmonization
  # ------------------------------------------------------------
  mvmr_models <- list(
    PF_VAT=c("pancreas_fat","VAT"),
    PF_ASAT=c("pancreas_fat","ASAT"),
    PF_VAT_ASAT=c("pancreas_fat","VAT","ASAT")
  )

  # Choose threshold for each exposure:
  # primary GWS; relaxed only if exposure has <3 instruments
  choose_inst <- function(tr) {
    main <- clump_results[[paste(tr,"GWS_5e8",sep="__")]]$clumped
    if (nrow(main)>=3) {
      attr(main,"threshold_used") <- "GWS_5e8"
      return(main)
    }
    rel <- clump_results[[paste(tr,"RELAXED_5e6",sep="__")]]$clumped
    attr(rel,"threshold_used") <- "RELAXED_5e6"
    rel
  }

  # Pull association for rsIDs and align to anchor effect allele
  get_assoc_aligned <- function(tr,snps,anchor_dt) {

    x <- copy(gwas[[tr]][SNP %in% snps])
    setorder(x,P)
    x <- x[!duplicated(SNP)]

    # anchor_dt本身已经使用REF_anchor/ALT_anchor命名，因此merge()不会
    # 自动给x侧REF/ALT/BETA/SE/AF添加"_x"后缀。这里在merge前显式重命名，
    # 保证后续列名固定，不依赖merge的suffix行为。
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

    present_old <- intersect(names(rename_map), names(x))
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
      "REF_x","ALT_x","BETA_x","SE_x","AF_x"
    )

    miss <- setdiff(required,names(m))
    if (length(miss)) {
      stop(
        "get_assoc_aligned缺少必要列: ",
        paste(miss,collapse=", "),
        "；trait=",tr
      )
    }

    same <- m$REF_anchor==m$REF_x & m$ALT_anchor==m$ALT_x
    swap_raw <- m$REF_anchor==m$ALT_x & m$ALT_anchor==m$REF_x
    csame_raw <- comp(m$REF_anchor)==m$REF_x & comp(m$ALT_anchor)==m$ALT_x
    cswap_raw <- comp(m$REF_anchor)==m$ALT_x & comp(m$ALT_anchor)==m$REF_x

    # Mutually exclusive allele-orientation precedence.
    swap <- (!same) & swap_raw
    csame <- (!same) & (!swap) & csame_raw
    cswap <- (!same) & (!swap) & (!csame) & cswap_raw

    ok <- !is.na(m$BETA_x) & !is.na(m$SE_x) &
      (same|swap|csame|cswap)
    flip <- swap|cswap

    m[, beta_aligned := ifelse(flip,-BETA_x,BETA_x)]
    m[, se_aligned := SE_x]
    m[, af_aligned := ifelse(flip,1-AF_x,AF_x)]
    m[, align_ok := ok]

    m
  }

  run_mvmr <- function(model_name,exposures,outcome="pancreas_iron") {

    inst_list <- lapply(exposures,choose_inst)
    names(inst_list) <- exposures

    threshold_used <- sapply(inst_list,function(z) attr(z,"threshold_used"))

    union_snps <- unique(unlist(lapply(inst_list,function(z) z$SNP)))
    union_snps <- union_snps[grepl("^rs[0-9]+$",union_snps)]

    if (!length(union_snps)) {
      return(list(status="NO_INSTRUMENTS"))
    }

    # anchor allele = first available allele among exposures
    allrows <- rbindlist(lapply(exposures,function(tr) {
      gwas[[tr]][SNP %in% union_snps,
                 .(SNP,REF,ALT,P)]
    }),fill=TRUE)

    setorder(allrows,P)
    anchor <- allrows[!duplicated(SNP),.(SNP,REF,ALT)]
    setnames(anchor,c("REF","ALT"),c("REF_anchor","ALT_anchor"))

    # exposure matrices
    mats <- list()

    for (tr in exposures) {
      a <- get_assoc_aligned(tr,union_snps,anchor)
      mats[[tr]] <- a[,.(SNP,beta_aligned,se_aligned,align_ok)]
      setnames(
        mats[[tr]],
        c("beta_aligned","se_aligned","align_ok"),
        c(
          paste0("bx_",tr),
          paste0("se_",tr),
          paste0("ok_",tr)
        )
      )
    }

    dat <- Reduce(
      function(x,y) merge(x,y,by="SNP",all=TRUE),
      c(list(anchor),mats)
    )

    # outcome
    oy <- get_assoc_aligned(outcome,union_snps,anchor)
    oy <- oy[,.(SNP,beta_aligned,se_aligned,align_ok)]
    setnames(oy,
             c("beta_aligned","se_aligned","align_ok"),
             c("by","se_y","ok_y"))

    dat <- merge(dat,oy,by="SNP",all.x=TRUE)

    okcols <- c(paste0("ok_",exposures),"ok_y")
    betacols <- c(paste0("bx_",exposures),"by")
    secols <- c(paste0("se_",exposures),"se_y")

    keep <- rep(TRUE,nrow(dat))
    for (cc in c(okcols,betacols,secols)) {
      keep <- keep & !is.na(dat[[cc]])
    }

    keep <- keep & dat$ok_y
    for (cc in paste0("ok_",exposures)) keep <- keep & dat[[cc]]

    dat <- dat[keep]

    # remove ambiguous anchor palindromic variants
    pal <- (
      (dat$REF_anchor=="A" & dat$ALT_anchor=="T") |
      (dat$REF_anchor=="T" & dat$ALT_anchor=="A") |
      (dat$REF_anchor=="C" & dat$ALT_anchor=="G") |
      (dat$REF_anchor=="G" & dat$ALT_anchor=="C")
    )

    # no common anchor frequency here -> conservatively remove all palindromic
    dat <- dat[!pal, ]

    if (nrow(dat) <= length(exposures)) {
      return(list(
        status="TOO_FEW_MVMR_INSTRUMENTS",
        harmonized=dat,
        threshold_used=threshold_used
      ))
    }

    bx_cols <- paste0("bx_", exposures)
    X <- as.matrix(dat[, ..bx_cols])
    y <- dat$by
    w <- 1/(dat$se_y^2)

    # weighted regression without intercept
    XtW <- t(X) %*% (X*w)
    XtWy <- t(X) %*% (y*w)

    beta <- tryCatch(
      as.numeric(solve(XtW,XtWy)),
      error=function(e) rep(NA_real_,ncol(X))
    )

    vc <- tryCatch(
      solve(XtW),
      error=function(e) matrix(NA_real_,ncol(X),ncol(X))
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
      p=p,
      threshold_used=threshold_used[exposures]
    )

    # Cochran-type residual Q
    resid <- y - as.numeric(X %*% beta)
    Q <- sum(w*resid^2)
    Q_df <- nrow(dat)-length(exposures)
    Q_p <- if (Q_df>0) pchisq(Q,Q_df,lower.tail=FALSE) else NA_real_

    qtab <- data.table(
      model=model_name,
      nsnp=nrow(dat),
      Q=Q,
      Q_df=Q_df,
      Q_p=Q_p
    )

    # Simple instrument-strength descriptive diagnostics
    # mean per-SNP F = mean((beta/se)^2) for each exposure
    ftab <- rbindlist(lapply(exposures,function(tr) {
      z2 <- (dat[[paste0("bx_",tr)]]/dat[[paste0("se_",tr)]])^2
      data.table(
        model=model_name,
        exposure=tr,
        nsnp=nrow(dat),
        mean_F=mean(z2,na.rm=TRUE),
        median_F=median(z2,na.rm=TRUE),
        min_F=min(z2,na.rm=TRUE)
      )
    }))

    list(
      status="OK",
      result=res,
      heterogeneity=qtab,
      F_proxy=ftab,
      harmonized=dat,
      threshold_used=threshold_used
    )
  }

  mvmr_res <- list()
  mvmr_q <- list()
  mvmr_f <- list()
  mvmr_qc <- list()

  for (nm in names(mvmr_models)) {

    rr <- run_mvmr(nm,mvmr_models[[nm]],"pancreas_iron")

    mvmr_qc[[nm]] <- data.table(
      model=nm,
      exposures=paste(mvmr_models[[nm]],collapse="+"),
      status=rr$status,
      nsnp=if (!is.null(rr$harmonized)) nrow(rr$harmonized) else NA_integer_
    )

    if (!is.null(rr$harmonized)) {
      fwrite(
        rr$harmonized,
        file.path(harmonized_dir,paste0("MVMR_",nm,"__harmonized.tsv")),
        sep="\t"
      )
    }

    if (identical(rr$status,"OK")) {
      mvmr_res[[nm]] <- rr$result
      mvmr_q[[nm]] <- rr$heterogeneity
      mvmr_f[[nm]] <- rr$F_proxy
    }
  }

  mvmr_qc_tab <- rbindlist(mvmr_qc,fill=TRUE)
  fwrite(mvmr_qc_tab,file.path(reports,"Step29_MVMR_QC.tsv"),sep="\t")

  if (length(mvmr_res))
    fwrite(rbindlist(mvmr_res,fill=TRUE),
           file.path(mvmr_dir,"Step29_MVMR_results.tsv"),sep="\t")

  if (length(mvmr_q))
    fwrite(rbindlist(mvmr_q,fill=TRUE),
           file.path(mvmr_dir,"Step29_MVMR_heterogeneity.tsv"),sep="\t")

  if (length(mvmr_f))
    fwrite(rbindlist(mvmr_f,fill=TRUE),
           file.path(mvmr_dir,"Step29_MVMR_F_proxy.tsv"),sep="\t")

  # ------------------------------------------------------------
  # 10. Summary
  # ------------------------------------------------------------
  summary <- data.table(
    metric=c(
      "n_traits",
      "n_uvmr_directions",
      "n_uvmr_thresholds",
      "n_mvmr_models",
      "n_mvmr_OK"
    ),
    value=c(
      nrow(traits),
      nrow(uvmr_specs),
      nrow(thresholds),
      length(mvmr_models),
      sum(mvmr_qc_tab$status=="OK")
    )
  )

  fwrite(summary,file.path(reports,"Step29_summary.tsv"),sep="\t")

  # ------------------------------------------------------------
  # 11. Methods / interpretation
  # ------------------------------------------------------------
  writeLines(
    c(
      "Step29 evaluates directional genetic evidence between pancreatic fat and pancreatic iron and performs multivariable MR adjusting for VAT and/or ASAT.",
      "Genome-wide significant instruments are selected at P<5e-8 and LD-clumped at r2<0.001 within 10,000 kb using an available European PLINK LD reference.",
      "A relaxed P<5e-6 threshold is generated only as a sensitivity instrument set and should not replace the primary genome-wide-significant analysis.",
      "Univariable MR includes IVW, weighted-median approximation, MR-Egger when instrument count permits, Cochran Q, Egger intercept and leave-one-out IVW.",
      "MVMR models evaluate pancreatic fat and VAT, pancreatic fat and ASAT, and pancreatic fat with both VAT and ASAT as simultaneous exposures for pancreatic iron.",
      "The MVMR estimator is weighted least squares without intercept using the union of exposure-specific clumped instruments.",
      "Because exact covariance information among exposure-SNP associations is unavailable, reported MVMR F statistics are descriptive per-exposure strength proxies rather than exact Sanderson-Windmeijer conditional F statistics.",
      "All traits are derived from the same UK Biobank abdominal MRI GWAS framework, so substantial sample overlap is expected.",
      "Therefore MR/MVMR findings are interpreted as directional genetic evidence under instrumental-variable assumptions, not as independent causal validation.",
      "The reverse pancreatic-iron to pancreatic-fat analysis may be underpowered if only one independent iron instrument is available."
    ),
    file.path(reports,"methods_and_interpretation.txt")
  )

  # ------------------------------------------------------------
  # 12. Completion flags
  # ------------------------------------------------------------
  flags <- character()

  if (clump_qc_tab[
      threshold_label=="GWS_5e8" &
      trait=="pancreas_fat", n_clumped] < 3) {
    flags <- c(flags,"PANCREAS_FAT_HAS_LT3_GWS_INSTRUMENTS")
  }

  if (clump_qc_tab[
      threshold_label=="GWS_5e8" &
      trait=="pancreas_iron", n_clumped] < 3) {
    flags <- c(flags,"PANCREAS_IRON_HAS_LT3_GWS_INSTRUMENTS")
  }

  if (any(mvmr_qc_tab$status!="OK")) {
    flags <- c(flags,"AT_LEAST_ONE_MVMR_MODEL_NOT_ESTIMABLE")
  }

  if (length(flags)) {
    writeLines("COMPLETED_STEP29_WITH_REVIEW_FLAGS",
               file.path(reports,"completion_status.txt"))
    writeLines(flags,file.path(reports,"review_flags.txt"))
  } else {
    writeLines("COMPLETED_STEP29_DIRECTIONAL_MR_MVMR",
               file.path(reports,"completion_status.txt"))
  }

  # ------------------------------------------------------------
  # 13. Review zip
  # ------------------------------------------------------------
  zipfile <- file.path(out,"29_MR_MVMR_review.zip")

  oldwd <- getwd()
  on.exit(setwd(oldwd),add=TRUE)
  setwd(out)

  files_to_zip <- c(
    list.files("reports",recursive=TRUE,full.names=TRUE),
    list.files("uvmr",recursive=TRUE,full.names=TRUE),
    list.files("mvmr",recursive=TRUE,full.names=TRUE),
    list.files("instruments",pattern="__clumped\\.tsv$",recursive=TRUE,full.names=TRUE),
    list.files("harmonized",recursive=TRUE,full.names=TRUE)
  )

  files_to_zip <- files_to_zip[file.exists(files_to_zip) & !dir.exists(files_to_zip)]

  utils::zip(zipfile,files=files_to_zip,flags="-q")

  # ------------------------------------------------------------
  # 14. Console
  # ------------------------------------------------------------
  cat("\n========================================\n")
  cat("STEP29 completed\n")
  cat("========================================\n\n")

  cat("Clumping QC:\n")
  print(clump_qc_tab)

  cat("\nUVMR results:\n")
  print(uvmr_res)

  cat("\nMVMR QC:\n")
  print(mvmr_qc_tab)

  if (length(mvmr_res)) {
    cat("\nMVMR results:\n")
    print(rbindlist(mvmr_res,fill=TRUE))
  }

  cat("\nOutput directory:\n",out,"\n",sep="")
  cat("Review zip:\n",zipfile,"\n",sep="")

})()
