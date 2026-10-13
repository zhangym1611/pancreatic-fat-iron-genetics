# PUBLICATION CURATION NOTE: Original analytical script numbered Step31 / Step31b.
# This copy is renamed Step30a / Step30b because historical Step30 was
# for manuscript figures/tables, not an additional statistical analysis.
# Only output directory labels, report filenames, completion markers,
# root selection and strict Step30a completion validation were updated.
# Original source version is preserved under 02_Provenance.
# ============================================================
# STEP 31
# Prioritized-locus glycemic-trait colocalization follow-up
#
# Goal:
# Test whether the two prioritized pancreatic fat–iron loci also
# share local genetic signals with fasting glucose, fasting insulin,
# HbA1c, and type 2 diabetes.
#
# Regions (GRCh37):
# chr8:129047378-130047378
# chr19:33785836-34633274
#
# Glycemic GWAS:
# fasting glucose  GCST90002232 (MAGIC 2021, EUR, N=200622)
# fasting insulin  GCST90002238 (MAGIC 2021, EUR, N=151013)
# HbA1c            GCST90002244 (MAGIC 2021, EUR, N=146806)
# type 2 diabetes  GCST006867   (Xue 2018, predominantly EUR, N=655666)
#
# Main coloc priors:
# p1=1e-4, p2=1e-4, p12=1e-5
# p12 sensitivity: 1e-6, 1e-5, 1e-4
#
# Any meaningful ABF result should undergo Step30b SuSiE confirmation.
# v2 fix: legacy GWAS Catalog studies may have harmonised-only TSV files;
#         downloader now prefers submitted files and falls back to harmonised files.
# v3 fix: resume mode reuses previously downloaded glycemic GWAS files from
#         earlier Step30a run folders and downloads only missing datasets.
# v4 fix: corrected Linux path detection for prior raw_glycemic files.
# v5 fix: legacy GWAS Catalog records are searched inside harmonised/;
#         GRCh37 formatted *.f.tsv.gz files are preferred and GRCh38
#         harmonised *.h.tsv.gz files are blocked for GRCh37 loci.
# v6 fix: harmonization now explicitly renames all pancreatic/glycemic
#         columns before merging; no dependence on automatic merge suffixes.
# v7 fix: allele-orientation labels are stored row-wise before compatible-row
#         filtering, preventing vector-length mismatch after harmonization QC.
# ============================================================

(function() {

  options(stringsAsFactors = FALSE)
  base <- path.expand(Sys.getenv("PANCREAS_PROJECT_ROOT", unset="~/胰腺"))

  lib <- path.expand("~/R_pancreas_step30a")
  dir.create(lib, recursive=TRUE, showWarnings=FALSE)
  .libPaths(c(lib, .libPaths()))

  need <- function(pkg) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
      install.packages(pkg, lib=lib, repos="https://cloud.r-project.org")
    }
    if (!requireNamespace(pkg, quietly=TRUE)) stop("Cannot load/install: ", pkg)
  }
  need("data.table")
  need("coloc")
  suppressPackageStartupMessages({
    library(data.table)
    library(coloc)
  })

  parent <- file.path(base, "30a_glycaemic_locus_followup")
  dir.create(parent, recursive=TRUE, showWarnings=FALSE)

  # ------------------------------------------------------------
  # Resume-aware cache
  # ------------------------------------------------------------
  # A stable cache avoids re-downloading large GWAS files across reruns.
  cache_dir <- file.path(parent, "shared_raw_cache")
  dir.create(cache_dir, recursive=TRUE, showWarnings=FALSE)

  # Recover raw files from any previous Step30a run folders.
  previous_raw <- list.files(
    parent,
    pattern="\\.(tsv|txt)(\\.gz)?$",
    recursive=TRUE,
    full.names=TRUE
  )

  previous_raw <- previous_raw[
    grepl("/raw_glycemic/", previous_raw, fixed=TRUE)
  ]

  if (length(previous_raw)) {
    for (f in previous_raw) {
      dest <- file.path(cache_dir, basename(f))
      if (!file.exists(dest) || file.info(dest)$size < file.info(f)$size) {
        ok <- file.copy(f, dest, overwrite=TRUE)
        if (ok) {
          message("[Step30a-v3] recovered cached GWAS: ", basename(f))
        }
      }
    }
  }

  out <- tempfile("glycemic_followup_", tmpdir=parent)
  dir.create(out)

  raw_dir <- file.path(out,"raw_glycemic")
  ext_dir <- file.path(out,"extracted")
  harm_dir <- file.path(out,"harmonized")
  coloc_dir <- file.path(out,"coloc")
  rep_dir <- file.path(out,"reports")
  log_dir <- file.path(out,"logs")

  for (d in c(raw_dir,ext_dir,harm_dir,coloc_dir,rep_dir,log_dir)) {
    dir.create(d, recursive=TRUE)
  }

  writeLines("INCOMPLETE", file.path(rep_dir,"completion_status.txt"))

  regions <- data.table(
    region_id=c("chr8_primary_shared_locus","chr19_secondary_shared_locus"),
    CHR=c(8L,19L),
    START=c(129047378L,33785836L),
    STOP=c(130047378L,34633274L)
  )
  fwrite(regions, file.path(rep_dir,"Step30a_regions.tsv"), sep="\t")

  glycemic <- data.table(
    trait=c("fasting_glucose","fasting_insulin","HbA1c","T2D"),
    accession=c("GCST90002232","GCST90002238","GCST90002244","GCST006867"),
    N=c(200622,151013,146806,655666),
    type=c("quant","quant","quant","cc"),
    ncase=c(NA,NA,NA,61714),
    ncontrol=c(NA,NA,NA,593952),
    ancestry=c("European","European","European","Predominantly European"),
    source=c("MAGIC Chen 2021","MAGIC Chen 2021","MAGIC Chen 2021","Xue 2018")
  )
  glycemic[, s := ifelse(type=="cc", ncase/(ncase+ncontrol), NA_real_)]
  fwrite(glycemic, file.path(rep_dir,"Step30a_glycemic_GWAS_metadata.tsv"), sep="\t")

  pancreatic <- data.table(
    trait=c("pancreas_fat","pancreas_iron"),
    accession=c("GCST90016675","GCST90016676"),
    N=c(25617,25617)
  )

  find_pan_vcf <- function(accession) {
    candidates <- c(
      file.path(base,paste0("ebi-a-",accession,".vcf.gz")),
      file.path(base,paste0(accession,".vcf.gz"))
    )
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) return(hit[1])
    rec <- list.files(base, pattern=paste0(accession,".*\\.vcf\\.gz$"),
                      recursive=TRUE, full.names=TRUE)
    if (length(rec)) return(rec[1])
    stop("Cannot locate pancreatic VCF for ", accession)
  }
  pancreatic[, vcf := vapply(accession, find_pan_vcf, character(1))]
  fwrite(pancreatic, file.path(rep_dir,"Step30a_pancreatic_GWAS_files.tsv"), sep="\t")

  curl <- Sys.which("curl")
  if (!nzchar(curl)) stop("curl not found")
  py <- Sys.which("python3")
  if (!nzchar(py)) py <- Sys.which("python")
  if (!nzchar(py)) stop("python not found")

  gwas_parent_bin <- function(acc) {
    nstr <- sub("^GCST","",acc)
    n <- as.numeric(nstr)
    width <- nchar(nstr)
    lo <- floor((n-1)/1000)*1000 + 1
    hi <- ceiling(n/1000)*1000
    paste0("GCST",sprintf(paste0("%0",width,"d"),lo),
           "-GCST",sprintf(paste0("%0",width,"d"),hi))
  }
  gwas_base_url <- function(acc) {
    paste0("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
           gwas_parent_bin(acc),"/",acc,"/")
  }

  discover_summary_file <- function(acc) {

    url <- gwas_base_url(acc)

    fetch_index <- function(index_url) {

      tmp <- tempfile(fileext=".html")

      rc <- system2(
        curl,
        c(
          "-L","--fail","--retry","3",
          "-o",shQuote(tmp),
          shQuote(index_url)
        )
      )

      if (!identical(as.integer(rc),0L)) {
        stop("Cannot access GWAS Catalog directory: ", index_url)
      }

      html <- paste(readLines(tmp,warn=FALSE),collapse="\n")
      unlink(tmp)

      hrefs <- unique(unlist(regmatches(
        html,
        gregexpr('href="[^"]+"',html,perl=TRUE)
      )))

      hrefs <- sub('^href="','',hrefs)
      hrefs <- sub('"$','',hrefs)

      hrefs
    }

    hrefs <- fetch_index(url)

    # ----------------------------------------------------------
    # 1) Prefer directly deposited structured TSV/TXT files
    # ----------------------------------------------------------
    direct_sumstats <- hrefs[
      grepl("\\.(tsv|txt)(\\.gz)?$",hrefs,ignore.case=TRUE) &
      !grepl("readme|metadata|yaml",hrefs,ignore.case=TRUE)
    ]

    if (length(direct_sumstats)) {

      raw_cand <- direct_sumstats[
        !grepl("harmonis|harmoniz|\\.h\\.tsv",direct_sumstats,ignore.case=TRUE)
      ]

      if (length(raw_cand)) {

        score <- as.integer(grepl(acc,raw_cand,fixed=TRUE))*10 +
          as.integer(grepl("GRCh37|buildGRCh37|build37|hg19",
                           raw_cand,ignore.case=TRUE))*5 +
          as.integer(grepl("\\.gz$",raw_cand,ignore.case=TRUE))*2

        raw_cand <- raw_cand[order(score,decreasing=TRUE)]

        message(
          "[Step30a-v5] using submitted/formatted top-level summary statistics for ",
          acc,": ",raw_cand[1]
        )

        return(list(
          url=paste0(url,raw_cand[1]),
          file_type="TOP_LEVEL_STRUCTURED_SUMSTATS"
        ))
      }
    }

    # ----------------------------------------------------------
    # 2) Some older GWAS Catalog records deposit an author .gz
    #    at top level and place standardized files in harmonised/.
    #
    #    IMPORTANT:
    #    - *.f.tsv.gz = formatted file, original build retained
    #    - *.h.tsv.gz = harmonised file, coordinates standardized
    #                   by GWAS Catalog to GRCh38
    #
    #    Our pancreatic loci are GRCh37, so prefer build37.f.tsv.gz.
    # ----------------------------------------------------------
    has_harmonised_dir <- any(
      hrefs == "harmonised/" |
      grepl("^harmonised/?$",hrefs,ignore.case=TRUE)
    )

    if (has_harmonised_dir) {

      hurl <- paste0(url,"harmonised/")
      hh <- fetch_index(hurl)

      formatted37 <- hh[
        grepl("\\.f\\.tsv\\.gz$",hh,ignore.case=TRUE) &
        grepl("build37|GRCh37|hg19",hh,ignore.case=TRUE)
      ]

      if (length(formatted37)) {

        score <- as.integer(grepl(acc,formatted37,fixed=TRUE))*10 +
          as.integer(grepl("build37",formatted37,ignore.case=TRUE))*5

        formatted37 <- formatted37[order(score,decreasing=TRUE)]

        message(
          "[Step30a-v5] no top-level structured TSV detected for ",acc,
          "; using GRCh37 GWAS Catalog formatted file: ",
          formatted37[1]
        )

        return(list(
          url=paste0(hurl,formatted37[1]),
          file_type="GWAS_CATALOG_FORMATTED_GRCh37"
        ))
      }

      # If no explicit build37 formatted file is present, stop rather
      # than silently use a GRCh38 harmonised file with GRCh37 regions.
      harmonised38 <- hh[
        grepl("\\.h\\.tsv\\.gz$",hh,ignore.case=TRUE)
      ]

      if (length(harmonised38)) {
        stop(
          "Only GRCh38 harmonised file(s) were found for ",acc,
          " but the current Step30a regions are GRCh37. ",
          "Automatic use was blocked to prevent coordinate-build mismatch."
        )
      }
    }

    # ----------------------------------------------------------
    # 3) Do not automatically use opaque author .gz files because
    #    their header/format may be study-specific. Report clearly.
    # ----------------------------------------------------------
    opaque_gz <- hrefs[
      grepl("\\.gz$",hrefs,ignore.case=TRUE) &
      !grepl("\\.(tsv|txt)\\.gz$",hrefs,ignore.case=TRUE)
    ]

    stop(
      "No directly usable GRCh37 structured summary-statistics file found for ",
      acc," in ",url,
      if (length(opaque_gz)) paste0(
        "\nOpaque author-contributed .gz file(s) present: ",
        paste(opaque_gz,collapse=", ")
      ) else "",
      "\nAvailable hrefs were:\n",
      paste(hrefs,collapse="\n")
    )
  }

  dl_qc <- list()

  for (i in seq_len(nrow(glycemic))) {

    tr <- glycemic$trait[i]
    acc <- glycemic$accession[i]

    # ----------------------------------------------------------
    # First try to reuse a previously downloaded file
    # ----------------------------------------------------------
    cached <- list.files(
      cache_dir,
      pattern=paste0("^",acc,".*\\.(tsv|txt)(\\.gz)?$"),
      full.names=TRUE,
      ignore.case=TRUE
    )

    cached <- cached[
      file.exists(cached) &
      file.info(cached)$size > 1000
    ]

    if (length(cached)) {

      # Prefer the largest matching cached file.
      cached <- cached[
        order(file.info(cached)$size, decreasing=TRUE)
      ][1]

      dest <- file.path(raw_dir, basename(cached))

      if (!file.exists(dest)) {
        file.copy(cached, dest, overwrite=TRUE)
      }

      message(
        "[Step30a-v3] reusing cached ",tr,
        " file: ",basename(cached)
      )

      dl_qc[[tr]] <- data.table(
        trait=tr,
        accession=acc,
        source_file_type="REUSED_PREVIOUS_DOWNLOAD",
        remote_url=NA_character_,
        local_file=dest,
        file_size_bytes=file.info(dest)$size
      )

      next
    }

    # ----------------------------------------------------------
    # Only missing datasets are discovered/downloaded
    # ----------------------------------------------------------
    message("[Step30a-v3] no cached file for ",tr,"; discovering ",acc)

    discovered <- discover_summary_file(acc)
    remote <- discovered$url
    source_file_type <- discovered$file_type

    dest <- file.path(raw_dir, basename(sub("\\?.*$","",remote)))

    if (!file.exists(dest) || file.info(dest)$size<1000) {

      message("[Step30a-v3] downloading ",remote)

      rc <- system2(
        curl,
        c(
          "-L","--fail","--retry","5","--retry-delay","5",
          "-o",shQuote(dest),
          shQuote(remote)
        )
      )

      if (!identical(as.integer(rc),0L)) {
        stop("Download failed for ",tr)
      }
    }

    # Save new download into persistent shared cache.
    cache_copy <- file.path(cache_dir, basename(dest))
    if (!file.exists(cache_copy) ||
        file.info(cache_copy)$size < file.info(dest)$size) {
      file.copy(dest, cache_copy, overwrite=TRUE)
    }

    dl_qc[[tr]] <- data.table(
      trait=tr,
      accession=acc,
      source_file_type=source_file_type,
      remote_url=remote,
      local_file=dest,
      file_size_bytes=file.info(dest)$size
    )
  }
  dl <- rbindlist(dl_qc)
  fwrite(dl,file.path(rep_dir,"Step30a_download_QC.tsv"),sep="\t")
  fwrite(data.table(cache_file=list.files(cache_dir,full.names=TRUE)),
         file.path(rep_dir,"Step30a_resume_cache_files.tsv"),sep="\t")

  pan_parser <- file.path(out,"extract_pancreatic_regions.py")
  writeLines(c(
'import gzip,sys,csv,math',
'vcf=sys.argv[1]; trait=sys.argv[2]; outfile=sys.argv[3]',
'regions=[(8,129047378,130047378,"chr8_primary_shared_locus"),(19,33785836,34633274,"chr19_secondary_shared_locus")]',
'def ff(x):',
'    try:',
'        v=float(x); return v if math.isfinite(v) else None',
'    except: return None',
'with gzip.open(vcf,"rt") as f, open(outfile,"w",newline="") as o:',
'    w=csv.writer(o,delimiter="\\t")',
'    w.writerow(["trait","region_id","CHR","BP","SNP","REF","ALT","BETA","SE","P","EAF","MAF"])',
'    for line in f:',
'        if line.startswith("#"): continue',
'        p=line.rstrip("\\n").split("\\t")',
'        if len(p)<10: continue',
'        ch=p[0].replace("chr","")',
'        if not ch.isdigit(): continue',
'        ch=int(ch)',
'        try: bp=int(p[1])',
'        except: continue',
'        rid=None',
'        for c,s,e,r in regions:',
'            if ch==c and s<=bp<=e: rid=r; break',
'        if rid is None: continue',
'        ref=p[3].upper(); alt=p[4].upper()',
'        if len(ref)!=1 or len(alt)!=1 or ref not in "ACGT" or alt not in "ACGT": continue',
'        keys=p[8].split(":"); vals=p[9].split(":")',
'        if len(keys)!=len(vals): continue',
'        d=dict(zip(keys,vals))',
'        beta=ff(d.get("ES")); se=ff(d.get("SE")); af=ff(d.get("AF"))',
'        if beta is None or se is None or se<=0 or af is None or not (0<af<1): continue',
'        lp=ff(d.get("LP")); pv=None',
'        if lp is not None:',
'            try: pv=10.0**(-lp)',
'            except: pv=None',
'        if pv is None: pv=math.erfc(abs(beta/se)/math.sqrt(2.0))',
'        snp=p[2] if p[2] not in ("",".") else f"chr{ch}:{bp}:{ref}:{alt}"',
'        w.writerow([trait,rid,ch,bp,snp,ref,alt,beta,se,pv,af,min(af,1-af)])'
  ),pan_parser,useBytes=TRUE)

  gly_parser <- file.path(out,"extract_glycemic_regions.py")
  writeLines(c(
'import gzip,sys,csv,math,re',
'fn=sys.argv[1]; trait=sys.argv[2]; outfile=sys.argv[3]',
'regions=[(8,129047378,130047378,"chr8_primary_shared_locus"),(19,33785836,34633274,"chr19_secondary_shared_locus")]',
'def norm(x): return re.sub(r"[^a-z0-9]+","_",x.strip().lower()).strip("_")',
'def pick(cols,aliases):',
'    for a in aliases:',
'        if a in cols: return cols.index(a)',
'    return None',
'def ff(x):',
'    try:',
'        v=float(x); return v if math.isfinite(v) else None',
'    except: return None',
'opener=gzip.open if fn.endswith(".gz") else open',
'with opener(fn,"rt",errors="replace") as f:',
'    header=f.readline().rstrip("\\n\\r")',
'    delim="\\t" if "\\t" in header else ","',
'    h=header.split(delim); cols=[norm(x) for x in h]',
'    i_chr=pick(cols,["chromosome","chr","chromosome_name"])',
'    i_bp=pick(cols,["base_pair_location","position","pos","bp","base_pair_position"])',
'    i_snp=pick(cols,["variant_id","snp","rsid","rs_id","markername","marker_name"])',
'    i_ea=pick(cols,["effect_allele","ea","a1","tested_allele","alt"])',
'    i_oa=pick(cols,["other_allele","non_effect_allele","oa","a2","reference_allele","ref"])',
'    i_beta=pick(cols,["beta","effect","effect_size","estimate","b"])',
'    i_or=pick(cols,["odds_ratio","or"])',
'    i_se=pick(cols,["standard_error","se","stderr","standarderror"])',
'    i_p=pick(cols,["p_value","pvalue","p_value_text","p","pval"])',
'    i_eaf=pick(cols,["effect_allele_frequency","eaf","effect_allele_freq","af","allele_frequency"])',
'    required={"CHR":i_chr,"BP":i_bp,"EA":i_ea,"OA":i_oa,"SE":i_se,"P":i_p}',
'    missing=[k for k,v in required.items() if v is None]',
'    if missing or (i_beta is None and i_or is None):',
'        sys.stderr.write("HEADER="+header+"\\n"); sys.stderr.write("NORMALIZED="+str(cols)+"\\n"); sys.exit(12)',
'    with open(outfile,"w",newline="") as o:',
'        w=csv.writer(o,delimiter="\\t")',
'        w.writerow(["trait","region_id","CHR","BP","SNP","EA","OA","BETA","SE","P","EAF","MAF"])',
'        reader=csv.reader(f,delimiter=delim)',
'        for p in reader:',
'            if len(p)<len(h): continue',
'            ch=str(p[i_chr]).replace("chr","").replace("CHR","")',
'            if not ch.isdigit(): continue',
'            ch=int(ch)',
'            try: bp=int(float(p[i_bp]))',
'            except: continue',
'            rid=None',
'            for c,s,e,r in regions:',
'                if ch==c and s<=bp<=e: rid=r; break',
'            if rid is None: continue',
'            ea=p[i_ea].upper(); oa=p[i_oa].upper()',
'            if len(ea)!=1 or len(oa)!=1 or ea not in "ACGT" or oa not in "ACGT": continue',
'            beta=ff(p[i_beta]) if i_beta is not None else None',
'            if beta is None and i_or is not None:',
'                orv=ff(p[i_or])',
'                if orv is not None and orv>0: beta=math.log(orv)',
'            se=ff(p[i_se]); pv=ff(p[i_p])',
'            if beta is None or se is None or se<=0 or pv is None or pv<0 or pv>1: continue',
'            eaf=ff(p[i_eaf]) if i_eaf is not None else None',
'            maf=min(eaf,1-eaf) if eaf is not None and 0<eaf<1 else None',
'            snp=p[i_snp] if i_snp is not None and p[i_snp] not in ("",".","NA") else f"chr{ch}:{bp}:{oa}:{ea}"',
'            w.writerow([trait,rid,ch,bp,snp,ea,oa,beta,se,pv,eaf if eaf is not None else "",maf if maf is not None else ""])'
  ),gly_parser,useBytes=TRUE)

  pan_dt <- list()
  for (i in seq_len(nrow(pancreatic))) {
    tr <- pancreatic$trait[i]
    fout <- file.path(ext_dir,paste0(tr,"_regions.tsv"))
    flog <- file.path(log_dir,paste0(tr,"_extract.log"))
    rc <- system2(py,c(shQuote(pan_parser),shQuote(pancreatic$vcf[i]),shQuote(tr),shQuote(fout)),
                  stdout=flog,stderr=flog)
    if (!identical(as.integer(rc),0L)) stop("Pancreatic extraction failed: ",tr)
    x <- fread(fout)
    x <- x[MAF>=0.01]
    setorder(x,region_id,P)
    x <- x[!duplicated(paste(region_id,SNP))]
    pan_dt[[tr]] <- x
  }

  gly_dt <- list()
  extract_qc <- list()
  for (i in seq_len(nrow(glycemic))) {
    tr <- glycemic$trait[i]
    fin <- dl[trait==tr,local_file][1]
    fout <- file.path(ext_dir,paste0(tr,"_regions.tsv"))
    flog <- file.path(log_dir,paste0(tr,"_extract.log"))
    message("[Step30a] extracting ",tr)
    rc <- system2(py,c(shQuote(gly_parser),shQuote(fin),shQuote(tr),shQuote(fout)),
                  stdout=flog,stderr=flog)
    if (!identical(as.integer(rc),0L)) stop("Glycemic extraction failed for ",tr,". See ",flog)
    x <- fread(fout)
    x <- x[is.na(MAF) | MAF>=0.01]
    setorder(x,region_id,P)
    x <- x[!duplicated(paste(region_id,SNP))]
    gly_dt[[tr]] <- x
    extract_qc[[tr]] <- x[,.(n_variants=.N,n_with_EAF=sum(!is.na(EAF)),min_p=min(P,na.rm=TRUE)),by=region_id][,trait:=tr]
  }
  fwrite(rbindlist(extract_qc,fill=TRUE),file.path(rep_dir,"Step30a_extraction_QC.tsv"),sep="\t")

  comp <- function(a) chartr("ACGT","TGCA",toupper(a))

  harmonize_pair <- function(pan,gly,pan_trait,gly_trait,reg_id) {

    a <- copy(pan[region_id==reg_id])
    b <- copy(gly[region_id==reg_id])

    if (!nrow(a) || !nrow(b)) return(data.table())

    # ----------------------------------------------------------
    # Explicitly rename columns BEFORE merge.
    # Do not rely on data.table/base merge suffix behavior.
    # ----------------------------------------------------------
    a2 <- a[, .(
      CHR_pan = CHR,
      BP_pan = BP,
      SNP_pan = SNP,
      REF_pan = REF,
      ALT_pan = ALT,
      BETA_pan = BETA,
      SE_pan = SE,
      P_pan = P,
      EAF_pan = EAF,
      MAF_pan = MAF
    )]

    b2 <- b[, .(
      CHR_gly = CHR,
      BP_gly = BP,
      SNP_gly = SNP,
      EA_gly = EA,
      OA_gly = OA,
      BETA_gly = BETA,
      SE_gly = SE,
      P_gly = P,
      EAF_gly = EAF,
      MAF_gly = MAF
    )]

    # ----------------------------------------------------------
    # 1) Primary merge by rsID when both sides have rsIDs
    # ----------------------------------------------------------
    ar <- a2[grepl("^rs[0-9]+$", SNP_pan)]
    br <- b2[grepl("^rs[0-9]+$", SNP_gly)]

    m1 <- data.table()

    if (nrow(ar) && nrow(br)) {

      m1 <- merge(
        ar,
        br,
        by.x="SNP_pan",
        by.y="SNP_gly",
        allow.cartesian=FALSE
      )

      if (nrow(m1)) {
        m1 <- m1[
          CHR_pan==CHR_gly & BP_pan==BP_gly
        ]

        if (nrow(m1)) {
          m1 <- m1[, .(
            CHR=CHR_pan,
            BP=BP_pan,
            SNP=SNP_pan,
            REF_pan,
            ALT_pan,
            BETA_pan,
            SE_pan,
            P_pan,
            EAF_pan,
            MAF_pan,
            EA_gly,
            OA_gly,
            BETA_gly,
            SE_gly,
            P_gly,
            EAF_gly,
            MAF_gly
          )]
        }
      }
    }

    # ----------------------------------------------------------
    # 2) Secondary merge by genomic position
    # ----------------------------------------------------------
    m2 <- merge(
      a2,
      b2,
      by.x=c("CHR_pan","BP_pan"),
      by.y=c("CHR_gly","BP_gly"),
      allow.cartesian=FALSE
    )

    if (nrow(m2)) {

      m2 <- m2[, .(
        CHR=CHR_pan,
        BP=BP_pan,
        SNP=ifelse(
          grepl("^rs[0-9]+$",SNP_pan),
          SNP_pan,
          SNP_gly
        ),
        REF_pan,
        ALT_pan,
        BETA_pan,
        SE_pan,
        P_pan,
        EAF_pan,
        MAF_pan,
        EA_gly,
        OA_gly,
        BETA_gly,
        SE_gly,
        P_gly,
        EAF_gly,
        MAF_gly
      )]
    }

    # Combine rsID and position matches.
    m <- rbindlist(list(m1,m2),fill=TRUE)

    if (!nrow(m)) return(m)

    # Remove duplicate variant matches.
    setorder(m,CHR,BP,SNP)
    m <- m[!duplicated(paste(CHR,BP,SNP,sep=":"))]

    # ----------------------------------------------------------
    # Allele harmonization
    #
    # Pancreatic VCF:
    #   REF_pan = non-effect/reference allele
    #   ALT_pan = effect allele
    #   BETA_pan = ALT effect
    #
    # Glycemic GWAS:
    #   EA_gly = effect allele
    #   OA_gly = other allele
    # ----------------------------------------------------------
    m[, REF_pan := toupper(REF_pan)]
    m[, ALT_pan := toupper(ALT_pan)]
    m[, EA_gly := toupper(EA_gly)]
    m[, OA_gly := toupper(OA_gly)]

    m <- m[
      REF_pan %chin% c("A","C","G","T") &
      ALT_pan %chin% c("A","C","G","T") &
      EA_gly %chin% c("A","C","G","T") &
      OA_gly %chin% c("A","C","G","T")
    ]

    if (!nrow(m)) return(m)

    # ----------------------------------------------------------
    # Compute orientation row-wise and STORE IT IN m BEFORE FILTERING.
    # This prevents vector-length mismatch after compatible rows are removed.
    # ----------------------------------------------------------
    m[, same_match := (
      ALT_pan == EA_gly &
      REF_pan == OA_gly
    )]

    m[, swap_raw := (
      ALT_pan == OA_gly &
      REF_pan == EA_gly
    )]

    m[, csame_raw := (
      comp(ALT_pan) == EA_gly &
      comp(REF_pan) == OA_gly
    )]

    m[, cswap_raw := (
      comp(ALT_pan) == OA_gly &
      comp(REF_pan) == EA_gly
    )]

    # Mutually exclusive precedence:
    # same > swap > complement_same > complement_swap
    m[, swap_match := (!same_match) & swap_raw]

    m[, csame_match := (
      (!same_match) &
      (!swap_match) &
      csame_raw
    )]

    m[, cswap_match := (
      (!same_match) &
      (!swap_match) &
      (!csame_match) &
      cswap_raw
    )]

    m[, compatible := (
      same_match |
      swap_match |
      csame_match |
      cswap_match
    )]

    m[, orientation := fifelse(
      same_match,
      "same",
      fifelse(
        swap_match,
        "swap",
        fifelse(
          csame_match,
          "complement_same",
          fifelse(
            cswap_match,
            "complement_swap",
            NA_character_
          )
        )
      )
    )]

    m[, flip_gly := orientation %chin% c(
      "swap",
      "complement_swap"
    )]

    # Now filter, after row-specific orientation is stored.
    m <- m[compatible==TRUE & !is.na(orientation)]

    if (!nrow(m)) return(m)

    # Remove temporary matching columns after filtering.
    m[, c(
      "same_match",
      "swap_raw",
      "csame_raw",
      "cswap_raw",
      "swap_match",
      "csame_match",
      "cswap_match"
    ) := NULL]

    m[, BETA_gly_aligned := ifelse(
      flip_gly,
      -BETA_gly,
      BETA_gly
    )]

    m[, EAF_gly_aligned := ifelse(
      is.na(EAF_gly),
      NA_real_,
      ifelse(flip_gly,1-EAF_gly,EAF_gly)
    )]

    # ----------------------------------------------------------
    # Palindromic ambiguity
    # ----------------------------------------------------------
    m[, palindromic := (
      (ALT_pan=="A" & REF_pan=="T") |
      (ALT_pan=="T" & REF_pan=="A") |
      (ALT_pan=="C" & REF_pan=="G") |
      (ALT_pan=="G" & REF_pan=="C")
    )]

    # Primary rule uses pancreatic MAF, which is consistently available.
    m[, ambiguous_pal := palindromic & MAF_pan>0.42]

    m <- m[ambiguous_pal==FALSE]

    # ----------------------------------------------------------
    # Frequency discrepancy QC where glycemic EAF is available
    # ----------------------------------------------------------
    m[, freq_diff := ifelse(
      is.na(EAF_gly_aligned),
      NA_real_,
      abs(EAF_pan-EAF_gly_aligned)
    )]

    # Do not auto-remove modest frequency differences.
    # Flag only extreme mismatch >0.20 for removal.
    m[, extreme_freq_mismatch := !is.na(freq_diff) & freq_diff>0.20]

    m <- m[extreme_freq_mismatch==FALSE]

    m[, pan_trait := pan_trait]
    m[, gly_trait := gly_trait]
    m[, region_id := reg_id]

    setorder(m,BP)

    m
  }

  main_results <- list()
  sens_results <- list()
  harm_qc <- list()
  p12_grid <- c(1e-6,1e-5,1e-4)

  for (reg in regions$region_id) {
    for (pt in pancreatic$trait) {
      for (gt in glycemic$trait) {

        message("[Step30a] ",reg," : ",pt," vs ",gt)
        h <- harmonize_pair(pan_dt[[pt]],gly_dt[[gt]],pt,gt,reg)

        fwrite(h,file.path(harm_dir,paste0(reg,"__",pt,"__",gt,".tsv")),sep="\t")

        harm_qc[[paste(reg,pt,gt,sep="__")]] <- data.table(
          region_id=reg,pancreatic_trait=pt,glycemic_trait=gt,
          n_harmonized=nrow(h),n_gly_EAF_available=sum(!is.na(h$EAF_gly_aligned))
        )

        if (nrow(h)<50) {
          for (p12 in p12_grid) {
            rr <- data.table(region_id=reg,pancreatic_trait=pt,glycemic_trait=gt,
                             p12=p12,status="LT50_COMMON_VARIANTS",n_common=nrow(h),
                             PP.H0=NA_real_,PP.H1=NA_real_,PP.H2=NA_real_,PP.H3=NA_real_,PP.H4=NA_real_)
            sens_results[[length(sens_results)+1]] <- rr
            if (identical(p12,1e-5)) main_results[[length(main_results)+1]] <- rr
          }
          next
        }

        d1 <- list(
          beta=h$BETA_pan,varbeta=h$SE_pan^2,type="quant",
          N=pancreatic[trait==pt,N][1],MAF=h$MAF_pan,
          snp=h$SNP,position=h$BP,sdY=1
        )

        gm <- glycemic[trait==gt]
        gly_maf <- ifelse(!is.na(h$EAF_gly_aligned),
                          pmin(h$EAF_gly_aligned,1-h$EAF_gly_aligned),
                          h$MAF_pan)

        if (gm$type[1]=="quant") {
          d2 <- list(
            beta=h$BETA_gly_aligned,varbeta=h$SE_gly^2,type="quant",
            N=gm$N[1],MAF=gly_maf,snp=h$SNP,position=h$BP
          )
        } else {
          d2 <- list(
            beta=h$BETA_gly_aligned,varbeta=h$SE_gly^2,type="cc",
            N=gm$N[1],s=gm$s[1],MAF=gly_maf,snp=h$SNP,position=h$BP
          )
        }

        for (p12 in p12_grid) {
          fit <- tryCatch(
            coloc::coloc.abf(d1,d2,p1=1e-4,p2=1e-4,p12=p12),
            error=function(e)e
          )

          if (inherits(fit,"error")) {
            rr <- data.table(region_id=reg,pancreatic_trait=pt,glycemic_trait=gt,
                             p12=p12,status=paste0("ERROR: ",conditionMessage(fit)),
                             n_common=nrow(h),
                             PP.H0=NA_real_,PP.H1=NA_real_,PP.H2=NA_real_,PP.H3=NA_real_,PP.H4=NA_real_)
          } else {
            ss <- fit$summary
            rr <- data.table(region_id=reg,pancreatic_trait=pt,glycemic_trait=gt,
                             p12=p12,status="OK",n_common=nrow(h),
                             PP.H0=as.numeric(ss["PP.H0.abf"]),
                             PP.H1=as.numeric(ss["PP.H1.abf"]),
                             PP.H2=as.numeric(ss["PP.H2.abf"]),
                             PP.H3=as.numeric(ss["PP.H3.abf"]),
                             PP.H4=as.numeric(ss["PP.H4.abf"]))
            if (identical(p12,1e-5)) {
              fwrite(as.data.table(fit$results),
                     file.path(coloc_dir,paste0(reg,"__",pt,"__",gt,"__SNP_H4_posterior.tsv")),
                     sep="\t")
            }
          }

          sens_results[[length(sens_results)+1]] <- rr
          if (identical(p12,1e-5)) main_results[[length(main_results)+1]] <- rr
        }
      }
    }
  }

  main <- rbindlist(main_results,fill=TRUE)
  sens <- rbindlist(sens_results,fill=TRUE)
  hqc <- rbindlist(harm_qc,fill=TRUE)

  fwrite(main,file.path(coloc_dir,"Step30a_glycemic_coloc_main.tsv"),sep="\t")
  fwrite(sens,file.path(coloc_dir,"Step30a_glycemic_coloc_p12_sensitivity.tsv"),sep="\t")
  fwrite(hqc,file.path(rep_dir,"Step30a_harmonization_QC.tsv"),sep="\t")

  main[,evidence_label:=fifelse(
    status!="OK","NOT_EVALUABLE",
    fifelse(PP.H4>=0.80,"STRONG_SHARED_SIGNAL",
    fifelse(PP.H4>=0.50,"MODERATE_SHARED_SIGNAL",
    fifelse(PP.H4>=0.20,"WEAK_SUGGESTIVE_SHARED_SIGNAL",
            "LITTLE_SHARED_SIGNAL_SUPPORT")))
  )]
  fwrite(main,file.path(rep_dir,"Step30a_main_evidence_labels.tsv"),sep="\t")

  susie_candidates <- main[
    status=="OK" & PP.H4>=0.50,
    .(region_id,pancreatic_trait,glycemic_trait,PP.H3,PP.H4,
      reason="PP.H4 >= 0.50 in main ABF; signal-level SuSiE confirmation recommended")
  ]
  fwrite(susie_candidates,file.path(rep_dir,"Step30a_SuSiE_followup_candidates.tsv"),sep="\t")

  summary_tab <- data.table(
    metric=c("n_regions","n_glycemic_traits","n_pancreatic_traits",
             "n_planned_pairwise_main_colocs","n_main_colocs_OK",
             "n_main_PP_H4_ge_0.5","n_main_PP_H4_ge_0.8"),
    value=c(nrow(regions),nrow(glycemic),nrow(pancreatic),
            16,sum(main$status=="OK"),
            sum(main$status=="OK"&main$PP.H4>=0.5,na.rm=TRUE),
            sum(main$status=="OK"&main$PP.H4>=0.8,na.rm=TRUE))
  )
  fwrite(summary_tab,file.path(rep_dir,"Step30a_summary.tsv"),sep="\t")

  writeLines(c(
    "Step30a is a targeted metabolic follow-up of two pre-prioritized pancreatic fat–iron loci.",
    "Glycemic traits were restricted a priori to fasting glucose, fasting insulin, HbA1c, and type 2 diabetes.",
    "Main priors: p1=p2=1e-4, p12=1e-5. Sensitivity: p12=1e-6 and 1e-4.",
    "Quantitative MAGIC traits use coloc sdY estimation from beta, varbeta, MAF and N.",
    "HbA1c must be interpreted cautiously because its genetics can reflect erythrocytic as well as glycemic biology.",
    "Fasting glucose and fasting insulin in Chen 2021 are BMI-adjusted.",
    "Any meaningful ABF signal should undergo Step30b SuSiE confirmation before biological interpretation."
  ),file.path(rep_dir,"methods_and_interpretation.txt"))

  if (nrow(main)==16 && all(main$status=="OK")) {
    writeLines("COMPLETED_STEP30A_ABF_GLYCEMIC_LOCUS_FOLLOWUP",
               file.path(rep_dir,"completion_status.txt"))
  } else {
    writeLines("COMPLETED_STEP30A_WITH_REVIEW_FLAGS",
               file.path(rep_dir,"completion_status.txt"))
  }

  zipfile <- file.path(out,"30a_glycaemic_locus_followup_review.zip")
  oldwd <- getwd()
  on.exit(setwd(oldwd),add=TRUE)
  setwd(out)
  files_to_zip <- c(
    list.files("reports",recursive=TRUE,full.names=TRUE),
    list.files("coloc",recursive=TRUE,full.names=TRUE),
    list.files("harmonized",pattern="\\.tsv$",recursive=TRUE,full.names=TRUE)
  )
  files_to_zip <- files_to_zip[file.exists(files_to_zip)&!dir.exists(files_to_zip)]
  utils::zip(zipfile,files=files_to_zip,flags="-q")

  cat("\n========================================\n")
  cat("STEP30A completed\n")
  cat("========================================\n\n")
  cat("Main coloc results (p12=1e-5):\n")
  print(main)
  cat("\nSummary:\n")
  print(summary_tab)
  if (nrow(susie_candidates)) {
    cat("\nRecommended Step30b SuSiE follow-up:\n")
    print(susie_candidates)
  } else {
    cat("\nNo main ABF comparison reached PP.H4 >= 0.50.\n")
  }
  cat("\nOutput directory:\n",out,"\n",sep="")
  cat("Review zip:\n",zipfile,"\n",sep="")

})()
