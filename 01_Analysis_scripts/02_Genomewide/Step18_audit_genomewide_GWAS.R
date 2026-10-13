# 第18步：8个GWAS全基因组输入审计（不进行MR、共定位或遗传相关估计）。
# 分块顺序扫描原VCF，内存不随全文件记录数增长；不按P值筛选。
# 网络仅用于可选GWAS Catalog元数据核查；失败不妨碍本地文件扫描。
(function() {
  if(!requireNamespace('data.table',quietly=TRUE)) stop('Missing R package: data.table')
  base <- path.expand('~/胰腺')
  registry <- data.frame(
    accession=paste0('GCST',c('90016666','90016669','90016671','90016672','90016673','90016674','90016675','90016676')),
    trait=c('liver_volume','pancreas_volume','visceral_adipose_tissue_volume',
      'abdominal_subcutaneous_adipose_tissue_volume','liver_fat','liver_iron','pancreas_fat','pancreas_iron'),
    phenotype_cn=c('肝脏体积','胰腺体积','内脏脂肪组织体积','腹部皮下脂肪组织体积',
      '肝脏脂肪','肝脏铁','胰腺脂肪','胰腺铁'),stringsAsFactors=FALSE)
  registry$id <- paste0('ebi-a-',registry$accession)
  registry$file <- file.path(base,paste0(registry$id,'.vcf.gz'))
  if(any(!file.exists(registry$file))) stop('Missing VCF files:\n',paste(registry$file[!file.exists(registry$file)],collapse='\n'))
  parent <- file.path(base,'18_genomewide_input_audit');dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile('GWAS_audit_',tmpdir=parent);if(!dir.create(out)) stop('Cannot create output directory')
  for(n in c('headers','catalog_metadata','per_study')) dir.create(file.path(out,n))
  sink(file.path(out,'run_log.txt'),split=TRUE);on.exit(sink(),add=TRUE)
  logs <- character();done <- FALSE
  on.exit({
    writeLines(if(length(logs)) logs else 'No captured warnings/errors.',file.path(out,'warnings_and_errors.txt'))
    writeLines(if(done) 'COMPLETED_REVIEW_REQUIRED_BEFORE_LDSC' else 'INCOMPLETE',file.path(out,'completion_status.txt'))
    writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  },add=TRUE)
  wt <- function(x,n) data.table::fwrite(x,file.path(out,n),sep='\t',na='NA')
  wt(registry,'study_registry.tsv')
  cat('Output:',out,'\nFull-genome audit; 8 files will be read sequentially.\n')
  # 环境清单：不打印环境变量或任何令牌，也不自动安装软件。
  commands <- c('python3','python','conda','mamba','micromamba','git','gzip','bcftools','tabix','plink','plink2')
  paths <- Sys.which(commands)
  wt(data.frame(command=commands,path=unname(paths),on_PATH=nzchar(paths)),'command_availability.tsv')
  pkgs <- c('data.table','Rsamtools','GenomicRanges','IRanges','httr','jsonlite','TwoSampleMR','ieugwasr','coloc')
  versions <- vapply(pkgs,function(p) if(requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else NA_character_,character(1))
  wt(data.frame(package=pkgs,version=versions),'R_package_versions.tsv')
  py <- unname(paths['python3']);if(!nzchar(py)) py <- unname(paths['python'])
  if(nzchar(py)) {
    pycode <- file.path(out,'python_environment_check.py')
    writeLines(c('import sys, importlib.util','print("python_version:",sys.version)',
      'print("executable:",sys.executable)',
      'for p in ["numpy","scipy","pandas","bitarray","pybedtools"]:',
      '    print(p, "AVAILABLE" if importlib.util.find_spec(p) is not None else "MISSING")'),pycode)
    v <- tryCatch(system2(py,shQuote(pycode),stdout=TRUE,stderr=TRUE),error=function(e) conditionMessage(e))
    writeLines(v,file.path(out,'python_environment.txt'))
  }
  dirs <- c(file.path(base,'ldsc'),path.expand('~/mr_tools/ldsc'),path.expand('~/ldsc'),
    file.path(base,'reference'),file.path(base,'references'),path.expand('~/mr_tools'))
  candidates <- unique(unlist(lapply(dirs[dir.exists(dirs)],function(d)
    list.files(d,pattern='ldsc\\.py$|munge_sumstats\\.py$|w_hm3|hapmap|ldscore|eur_w_ld',
      full.names=TRUE,recursive=FALSE,ignore.case=TRUE))))
  wt(data.frame(path=as.character(candidates)),'known_location_reference_inventory.tsv')
  # 仅检查指定目录，不递归搜索整个主目录；空结果不代表服务器全局没有资源。
  info <- file.info(registry$file)
  wt(data.frame(id=registry$id,file=registry$file,bytes=info$size,mtime=as.character(info$mtime),
    tbi_exists=file.exists(paste0(registry$file,'.tbi'))),'source_file_inventory.tsv')
  meta <- list();audits <- list();run_status <- list()
  extract_attr <- function(s,key) {
    pat <- paste0('(?:<|,)',key,'=([^,>]+)')
    m <- regexec(pat,s,perl=TRUE);a <- regmatches(s,m)[[1L]]
    if(length(a)>=2L) a[2L] else NA_character_
  }
  # 单个文件完整扫描；最终输出计数是逐步筛选漏斗，不把重叠排除原因相加。
  scan_one <- function(f,id) {
    con <- gzfile(f,'rt');on.exit(close(con),add=TRUE)
    header <- character()
    repeat {
      h <- readLines(con,n=1L,warn=FALSE)
      if(!length(h)) stop('EOF before #CHROM header')
      header <- c(header,h)
      if(startsWith(h,'#CHROM')) break
      if(length(header)>10000L) stop('Unreasonable header length')
    }
    writeLines(header,file.path(out,'headers',paste0(id,'_header.txt')))
    cols <- strsplit(tail(header,1L),'\t',fixed=TRUE)[[1L]]
    if(length(cols)!=10L||cols[10L]!=id) stop('Unexpected GWAS-VCF sample columns')
    samp <- header[startsWith(header,'##SAMPLE=')]
    if(length(samp)!=1L||extract_attr(samp,'ID')!=id) stop('Unexpected SAMPLE header')
    type <- extract_attr(samp,'StudyType');N <- suppressWarnings(as.numeric(extract_attr(samp,'TotalControls')))
    expected_records <- suppressWarnings(as.numeric(extract_attr(samp,'HarmonisedVariants')))
    if(is.na(type)||type!='Continuous'||!is.finite(N)||N<=0) stop('Invalid study type or N')
    if(!any(grepl('assembly=HG19/GRCh37',header,fixed=TRUE))) stop('GRCh37 not confirmed in VCF header')
    es_header <- header[startsWith(header,'##FORMAT=<ID=ES,')]
    if(length(es_header)!=1L||!grepl('alternative allele',es_header,fixed=TRUE)) stop('ALT effect convention not confirmed')
    fmts <- list();chr_counts <- list();examples <- list()
    totals <- c(records=0,autosomal=0,autosomal_SNV=0,valid_effect=0,valid_LP=0,
      valid_AF=0,filter_accepted=0,single_rsid=0,nonpalindromic=0,MAF_ge_001=0)
    missing <- c(ES=0,SE=0,LP=0,AF=0,SS_present=0,SS_valid=0,SI_present=0,SI_valid=0)
    z2_sum <- 0;z2_n <- 0;ss_range <- c(Inf,-Inf);si_range <- c(Inf,-Inf)
    block <- 0L
    repeat {
      lines <- readLines(con,n=50000L,warn=FALSE)
      if(!length(lines)) break
      lines <- lines[nzchar(lines)]
      if(!length(lines)) next
      if(any(startsWith(lines,'#'))) stop('Unexpected header in variant body')
      d <- as.data.frame(data.table::fread(text=paste(lines,collapse='\n'),sep='\t',header=FALSE,
        colClasses='character',na.strings='.',quote='',showProgress=FALSE))
      if(ncol(d)!=10L) stop('VCF records do not have ten columns')
      names(d) <- c('chr','pos','rsid','REF','ALT','QUAL','FILTER','INFO','FORMAT','VALUE')
      if(anyNA(d$FORMAT)) stop('Missing FORMAT field')
      block <- block+1L
      tb <- table(d$FORMAT);fmts[[block]] <- data.frame(FORMAT=names(tb),records=as.numeric(tb))
      cb <- table(d$chr);chr_counts[[block]] <- data.frame(chr=names(cb),records=as.numeric(cb))
      val <- matrix(NA_real_,nrow(d),6L,dimnames=list(NULL,c('ES','SE','LP','AF','SS','SI')))
      for(fmt in unique(d$FORMAT)) {
        ix <- which(d$FORMAT==fmt);keys <- strsplit(fmt,':',fixed=TRUE)[[1L]]
        if(anyDuplicated(keys)||!all(c('ES','SE','LP','AF') %in% keys)) stop('Unsupported FORMAT: ',fmt)
        parts <- data.table::tstrsplit(d$VALUE[ix],':',fixed=TRUE,fill=NA_character_)
        if(length(parts)!=length(keys)) stop('FORMAT/value length mismatch')
        for(nm in intersect(colnames(val),keys)) val[ix,nm] <- suppressWarnings(as.numeric(parts[[match(nm,keys)]]))
        for(nm in c('SS','SI')) if(nm %in% keys) missing[paste0(nm,'_present')] <- missing[paste0(nm,'_present')]+length(ix)
      }
      for(nm in c('ES','SE','LP','AF')) missing[nm] <- missing[nm]+sum(!is.finite(val[,nm]))
      ss <- val[,'SS'];si <- val[,'SI']
      vss <- is.finite(ss)&ss>0;vsi <- is.finite(si)
      missing['SS_valid'] <- missing['SS_valid']+sum(vss);missing['SI_valid'] <- missing['SI_valid']+sum(vsi)
      if(any(vss)) ss_range <- range(c(ss_range,ss[vss])[is.finite(c(ss_range,ss[vss]))])
      if(any(vsi)) si_range <- range(c(si_range,si[vsi])[is.finite(c(si_range,si[vsi]))])
      pos <- suppressWarnings(as.numeric(d$pos))
      a <- d$chr %in% as.character(1:22) & is.finite(pos) & pos>0 & pos==floor(pos)
      snv <- a & d$REF %in% c('A','C','G','T') & d$ALT %in% c('A','C','G','T') & d$REF!=d$ALT
      z <- val[,'ES']/val[,'SE']
      eff <- snv & is.finite(val[,'ES']) & is.finite(val[,'SE']) & val[,'SE']>0 & is.finite(z)
      lp <- eff & is.finite(val[,'LP']) & val[,'LP']>=0
      af <- lp & is.finite(val[,'AF']) & val[,'AF']>0 & val[,'AF']<1
      filt <- af & (is.na(d$FILTER)|d$FILTER=='PASS')
      rs <- filt & !is.na(d$rsid) & grepl('^rs[0-9]+$',d$rsid)
      pal <- paste0(d$REF,d$ALT) %in% c('AT','TA','CG','GC')
      np <- rs & !pal
      common <- np & pmin(val[,'AF'],1-val[,'AF'])>=0.01
      masks <- list(records=rep(TRUE,nrow(d)),autosomal=a,autosomal_SNV=snv,valid_effect=eff,
        valid_LP=lp,valid_AF=af,filter_accepted=filt,single_rsid=rs,nonpalindromic=np,MAF_ge_001=common)
      for(nm in names(masks)) {masks[[nm]][is.na(masks[[nm]])] <- FALSE;totals[nm] <- totals[nm]+sum(masks[[nm]])}
      ii <- which(masks$MAF_ge_001)
      if(length(ii)) {
        z2_sum <- z2_sum+sum(z[ii]^2);z2_n <- z2_n+length(ii)
        if(length(examples)<1L) {
          j <- head(ii,10L)
          examples[[1L]] <- data.frame(SNP=d$rsid[j],CHR=d$chr[j],BP=pos[j],A1=d$ALT[j],A2=d$REF[j],
            BETA=val[j,'ES'],SE=val[j,'SE'],Z=z[j],LP=val[j,'LP'],EAF=val[j,'AF'],
            N_study=N,N_variant=val[j,'SS'])
        }
      }
      if(block %% 10L==0L) {cat(id,': ',format(totals['records'],big.mark=','),' records read\n',sep='');flush.console()}
    }
    fmt <- aggregate(records~FORMAT,do.call(rbind,fmts),sum)
    chr <- aggregate(records~chr,do.call(rbind,chr_counts),sum)
    wt(fmt,paste0('per_study/',id,'_FORMAT_counts.tsv'))
    wt(chr,paste0('per_study/',id,'_chromosome_counts.tsv'))
    wt(data.frame(stage=names(totals),n_records=unname(totals)),paste0('per_study/',id,'_QC_funnel.tsv'))
    if(length(examples)) wt(examples[[1L]],paste0('per_study/',id,'_signed_effect_examples.tsv'))
    row <- data.frame(id=id,study_type=type,build='GRCh37',N_study=N,
      N_source='VCF SAMPLE TotalControls for continuous trait; not case-control controls',
      header_HarmonisedVariants=expected_records,as.list(totals),
      count_matches_header=is.finite(expected_records)&&totals['records']==expected_records,
      mean_Z2_preliminary=if(z2_n) z2_sum/z2_n else NA_real_,
      SS_field_records=missing['SS_present'],SS_valid_records=missing['SS_valid'],
      SS_min=if(is.finite(ss_range[1])) ss_range[1] else NA_real_,SS_max=if(is.finite(ss_range[2])) ss_range[2] else NA_real_,
      SI_field_records=missing['SI_present'],SI_valid_records=missing['SI_valid'],
      SI_min=if(is.finite(si_range[1])) si_range[1] else NA_real_,SI_max=if(is.finite(si_range[2])) si_range[2] else NA_real_,
      effect_allele='ALT',other_allele='REF',effect_scale='NOT_CONFIRMED_FROM_VCF',
      duplicate_SNP_check='DEFERRED_TO_HAPMAP3_HARMONISATION',row.names=NULL,check.names=FALSE)
    wt(data.frame(field=names(missing),count=unname(missing)),paste0('per_study/',id,'_field_availability.tsv'))
    row
  }
  for(i in seq_len(nrow(registry))) {
    id <- registry$id[i];cat('\n===== ',id,' / ',registry$trait[i],' =====\n',sep='')
    # 保存Catalog原始响应便于核对表型、来源和祖源；不猜测效应单位。
    cm <- data.frame(id=id,catalog_trait=NA_character_,initial_sample_size=NA_character_,status='NOT_QUERIED')
    if(requireNamespace('httr',quietly=TRUE)&&requireNamespace('jsonlite',quietly=TRUE)) {
      tryCatch({
        u <- paste0('https://www.ebi.ac.uk/gwas/rest/api/studies/',registry$accession[i])
        h <- httr::GET(u,httr::timeout(25));httr::stop_for_status(h)
        txt <- httr::content(h,as='text',encoding='UTF-8');obj <- jsonlite::fromJSON(txt,simplifyVector=FALSE)
        writeLines(txt,file.path(out,'catalog_metadata',paste0(id,'.json')))
        if(!is.null(obj$diseaseTrait$trait)) cm$catalog_trait <- as.character(obj$diseaseTrait$trait)
        if(!is.null(obj$initialSampleSize)) cm$initial_sample_size <- as.character(obj$initialSampleSize)
        cm$status <- if(is.na(cm$catalog_trait)) 'RESPONSE_NEEDS_REVIEW' else 'RETRIEVED_REVIEW_REQUIRED'
      },error=function(e) {cm$status <<- 'LOOKUP_FAILED';logs <<- c(logs,paste(id,'Catalog:',conditionMessage(e)))})
    }
    meta[[id]] <- cm;wt(do.call(rbind,meta),'catalog_metadata_summary.tsv')
    tryCatch(withCallingHandlers({
      row <- scan_one(registry$file[i],id);audits[[id]] <- row
      run_status[[id]] <- data.frame(id=id,status=if(row$count_matches_header) 'SCANNED_COUNT_MATCH' else 'SCANNED_COUNT_MISMATCH',message='')
      wt(row,paste0('per_study/',id,'_summary.tsv'))
    },warning=function(w) {logs <<- c(logs,paste(id,conditionMessage(w),sep=': '));invokeRestart('muffleWarning')}),
      error=function(e) {run_status[[id]] <<- data.frame(id=id,status='ERROR',message=conditionMessage(e));logs <<- c(logs,paste(id,'ERROR',conditionMessage(e)))})
    wt(do.call(rbind,run_status),'run_status.tsv')
    if(length(audits)) wt(do.call(rbind,audits),'genomewide_GWAS_QC_summary.tsv')
  }
  if(any(vapply(run_status,function(x) x$status!='SCANNED_COUNT_MATCH',logical(1))))
    stop('Some studies failed or counts differ. Inspect run_status.tsv; successful files have saved summaries.')
  writeLines(c(
    'This step audits inputs; it does not calculate SNP heritability, genetic correlation, MR, or new chr8 results.',
    'QC funnel counts are sequential; no association-P threshold or LD clumping is applied.',
    'MAF>=0.01 is a preliminary count only; final SNP selection requires validated HapMap3 and European LD-score resources.',
    'Rows passing this audit are not yet final LDSC inputs: duplicate rsIDs and reference-allele compatibility remain to be checked.',
    'Mean Z squared is a preliminary diagnostic, not an LD-score intercept, heritability, or genetic correlation.',
    'VCF study-level sample sizes do not establish variant-specific effective N or exact overlap between traits.',
    'SI is the GWAS-VCF summary-data imputation field; it must not automatically be treated as original genotype imputation INFO.',
    'Phenotype transformations, GWAS covariates, ancestry, and original genotype QC require source-method/supplement review.',
    'No sdY=1 or common effect scale is inferred from these files.',
    'With overlapping UK Biobank samples, later bivariate LDSC intercept must not be fixed to zero by default.',
    'Known-location software inventory is limited; absence there does not prove software is absent on the server.',
    'After review, prepare harmonised genome-wide inputs, assess single-trait SNP heritability precision, then estimate cross-trait correlations.'),
    file.path(out,'interpretation_and_remaining_checks.txt'))
  s <- file.path(base,'18_audit_genomewide_GWAS.R');if(file.exists(s)) file.copy(s,file.path(out,basename(s)))
  cat('\n===== Genome-wide audit summary =====\n')
  res <- do.call(rbind,audits)
  print(res[,c('id','N_study','records','count_matches_header','MAF_ge_001','SS_field_records','effect_scale')],row.names=FALSE)
  done <- TRUE
  cat('\nCOMPLETED_REVIEW_REQUIRED_BEFORE_LDSC\nOutput:',out,'\nUpload this entire result folder as ZIP.\n')
})()
