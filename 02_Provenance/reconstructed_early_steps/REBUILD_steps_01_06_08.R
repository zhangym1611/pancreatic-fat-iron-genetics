# Reconstructed 2026-10-07 from archived inputs, results and logs.
# NOT the original historical source. New outputs only; no JWT/network needed.
(function() {
  # ---------------------- Configuration ----------------------
  base <- path.expand('~/胰腺')
  plink <- path.expand('~/mr_tools/plink')
  ref <- file.path(base, '00_ld_reference', 'EUR')
  strict_historical_versions <- TRUE
  seed <- 20261003L  # recovered from the original main MR RDS
  ids <- c(pancreas_fat='ebi-a-GCST90016675',
           pancreas_iron='ebi-a-GCST90016676',
           pancreas_volume='ebi-a-GCST90016669')
  Ns <- c(pancreas_fat=25617L, pancreas_iron=25617L, pancreas_volume=31758L)
  region_start <- 129047378L; region_end <- 130047378L
  scripts <- vapply(sys.frames(), function(e) {
    x <- e$ofile; if (is.null(x)) '' else as.character(x)[1L]
  }, character(1))
  scripts <- scripts[nzchar(scripts)]
  if (!length(scripts)) stop('Run this saved file with source(), not by pasting its contents.')
  script_file <- normalizePath(tail(scripts,1L), mustWork=TRUE)
  package_dir <- dirname(script_file)
  evidence <- file.path(package_dir,'reference_results')
  if (!dir.exists(evidence)) stop('Keep reference_results next to this script; unzip the entire package.')
  pkgs <- c('data.table','Rsamtools','GenomicRanges','IRanges','TwoSampleMR','coloc')
  ok <- vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)
  if (!all(ok)) stop('Missing packages: ',paste(pkgs[!ok],collapse=', '))
  historical <- c(TwoSampleMR='0.7.9',coloc='6.0.1')
  actual <- vapply(names(historical),function(p) as.character(utils::packageVersion(p)),character(1))
  if(strict_historical_versions && any(actual!=historical))
    stop('Historical package versions required: TwoSampleMR 0.7.9 and coloc 6.0.1. Installed: ',paste(names(actual),actual,collapse='; '),'. No software was changed.')
  vcfs <- setNames(file.path(base,paste0(ids,'.vcf.gz')),names(ids))
  required <- c(vcfs,paste0(vcfs,'.tbi'),plink,paste0(ref,c('.bed','.bim','.fam')))
  if(!all(file.exists(required))) stop('Missing inputs:\n',paste(required[!file.exists(required)],collapse='\n'))
  if(file.access(plink,1)!=0L) stop('PLINK is not executable.')
  pv <- system2(plink,'--version',stdout=TRUE,stderr=TRUE)
  if(!any(grepl('PLINK v1.90b6.21',pv,fixed=TRUE))) stop('Expected historical PLINK v1.90b6.21; found: ',paste(pv,collapse=' '))
  parent <- file.path(base,'reconstructed_early_steps')
  dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile('rebuild_',tmpdir=parent);stopifnot(dir.create(out))
  report <- file.path(out,'reports');dir.create(report)
  sink(file.path(report,'run_log.txt'),split=TRUE)
  on.exit(sink(),add=TRUE)
  complete <- FALSE; warning_messages <- character(); checks <- list()
  wt <- function(d,f) data.table::fwrite(d,f,sep='\t',na='NA')
  on.exit({
    writeLines(capture.output(sessionInfo()),file.path(report,'sessionInfo.txt'))
    writeLines(unique(warning_messages),file.path(report,'warnings.txt'))
    if(length(checks)) wt(data.table::rbindlist(checks,fill=TRUE),file.path(report,'historical_comparison.tsv'))
    writeLines(if(complete) 'COMPLETED_RECONSTRUCTION_REVIEW_COMPARISONS' else 'INCOMPLETE_CHECK_LOG',file.path(report,'completion_status.txt'))
    # Small review ZIP; complete working results stay in out.
    zb <- Sys.which('zip')
    if(nzchar(zb)) {
      wd <- getwd()
      tryCatch({
        setwd(out)
        folders <- c('reports','01_candidate_snps','02_ld_clumping','03_outcome_extraction','04_harmonisation','05_MR_results','06_colocalisation','08_pancreas_volume_coloc')
        folders <- folders[dir.exists(folders)]
        rc <- system2(zb,c('-q','-r',shQuote(file.path(out,'reconstruction_review.zip')),vapply(folders,shQuote,character(1))))
        if(rc!=0L) cat('Review ZIP failed; output files remain available.\n')
      },finally=setwd(wd))
    }
  },add=TRUE)
  file.copy(script_file,file.path(report,basename(script_file)))
  saveRDS(list(ids=ids,N=Ns,build='GRCh37',region=c(region_start,region_end),
    seed=seed,source='RECONSTRUCTED',plink_version=pv,ref=ref,vcfs=vcfs,
    package_versions=actual),file.path(report,'configuration.rds'))
  cat('Output:',out,'\nThis is reconstruction, not independent replication.\n')
  # Checks never change analysis values to match the old results.
  check <- function(label,status,detail='') {
    checks[[length(checks)+1L]] <<- data.frame(check=label,status=status,detail=as.character(detail))
    wt(data.table::rbindlist(checks,fill=TRUE),file.path(report,'historical_comparison.tsv'))
    cat('CHECK',label,status,detail,'\n')
  }
  mkdir <- function(n) {p<-file.path(out,n);dir.create(p);p}
  baseline <- function(suffix) {
    fs <- list.files(evidence,recursive=TRUE,full.names=TRUE)
    fs <- fs[endsWith(fs,suffix)]
    if(length(fs)!=1L) stop('Reference result missing or ambiguous: ',suffix)
    as.data.frame(data.table::fread(fs))
  }
  compare <- function(label,new,suffix,keys,cols,stochastic=character()) {
    old<-baseline(suffix);new<-as.data.frame(new)
    if(!all(c(keys,cols)%in%names(new)) || !all(c(keys,cols)%in%names(old))) stop('Comparison schema mismatch: ',label)
    keyfun<-function(d) do.call(paste,c(lapply(d[keys],as.character),sep='|'))
    kn<-keyfun(new);ko<-keyfun(old)
    if(anyDuplicated(kn)||anyDuplicated(ko)) stop('Nonunique comparison keys: ',label)
    if(!setequal(kn,ko)) {
      check(label,'DIFFERENT_VARIANT_OR_ROW_SET',paste('new',nrow(new),'old',nrow(old)))
      return(invisible(FALSE))
    }
    new<-new[match(ko,kn),,drop=FALSE]
    for(cc in cols) {
      x<-new[[cc]];y<-old[[cc]]
      if(is.numeric(x)&&is.numeric(y)) {
        same<- (is.na(x)&is.na(y)) | (!is.na(x)&!is.na(y)&abs(x-y)<=1e-10+1e-6*abs(y))
        same[is.na(same)]<-FALSE
        delta<-if(any(is.finite(x)&is.finite(y))) max(abs(x-y),na.rm=TRUE) else NA_real_
      } else {
        same<-(is.na(x)&is.na(y))|(!is.na(x)&!is.na(y)&as.character(x)==as.character(y));same[is.na(same)]<-FALSE;delta<-NA_real_
      }
      check(paste(label,cc,sep=':'),if(all(same)) 'MATCH' else if(cc%in%stochastic) 'BOOTSTRAP_DIFFERENCE_REVIEW' else 'DIFFERENCE_REVIEW',paste('different',sum(!same),'max_abs_diff',delta))
    }
  }
  # ---------------------- VCF helpers ----------------------
  metadata <- lapply(names(ids),function(tr) {
    co<-gzfile(vcfs[[tr]],'rt');on.exit(close(co))
    h<-character()
    repeat {s<-readLines(co,n=1L,warn=FALSE);if(!length(s))stop('VCF header incomplete');h<-c(h,s);if(startsWith(s,'#CHROM'))break}
    ss<-h[startsWith(h,'##SAMPLE=')]
    if(length(ss)!=1L || !grepl(ids[[tr]],ss,fixed=TRUE) || !grepl('StudyType=Continuous',ss,fixed=TRUE)) stop('Study metadata mismatch: ',tr)
    N<-as.numeric(sub('.*TotalControls=([0-9.]+).*','\\1',ss))
    if(!is.finite(N)||N!=Ns[[tr]]||!any(grepl('assembly=HG19/GRCh37',h,fixed=TRUE))) stop('N/build mismatch: ',tr)
    writeLines(h,file.path(report,paste0(ids[[tr]],'_header.txt')))
    list(N=N,header=h)
  });names(metadata)<-names(ids)
  parse_vcf <- function(lines) {
    lines<-lines[nzchar(lines)&!startsWith(lines,'#')]
    if(!length(lines)) return(data.table::data.table())
    d<-data.table::fread(text=paste(lines,collapse='\n'),sep='\t',header=FALSE,
      colClasses='character',na.strings='.',quote='',showProgress=FALSE)
    if(ncol(d)!=10L)stop('VCF must have exactly one study column.')
    data.table::setnames(d,c('chr','pos','rsid','REF','ALT','QUAL','FILTER','INFO','FORMAT','VALUE'))
    d$pos<-suppressWarnings(as.integer(d$pos))
    for(nm in c('ES','SE','LP','AF'))d[[nm]]<-NA_real_
    if(anyNA(d$FORMAT)||anyNA(d$VALUE))stop('Missing FORMAT or VALUE.')
    for(fmt in unique(d$FORMAT)) {
      ii<-which(d$FORMAT==fmt);kk<-strsplit(fmt,':',fixed=TRUE)[[1L]]
      if(anyDuplicated(kk)||!all(c('ES','SE','LP','AF')%in%kk))stop('Unsupported FORMAT: ',fmt)
      vals<-data.table::tstrsplit(d$VALUE[ii],':',fixed=TRUE,fill=NA_character_)
      if(length(vals)!=length(kk))stop('Unexpected FORMAT length.')
      for(nm in c('ES','SE','LP','AF'))d[[nm]][ii]<-suppressWarnings(as.numeric(vals[[match(nm,kk)]]))
    }
    d$key<-paste(d$chr,d$pos,d$REF,d$ALT,sep=':');d
  }
  as_effect <- function(d) data.frame(SNP=d$rsid,chr=d$chr,pos=d$pos,
    effect_allele=d$ALT,other_allele=d$REF,beta=d$ES,se=d$SE,eaf=d$AF,
    LP=d$LP,pval=10^(-d$LP),F_approx=(d$ES/d$SE)^2,stringsAsFactors=FALSE)
  valid_effect <- function(d) is.finite(d$ES)&is.finite(d$SE)&d$SE>0&
    is.finite(d$LP)&d$LP>=0&(is.na(d$FILTER)|d$FILTER=='PASS')
  read_query <- function(tr,gr) {
    tb<-Rsamtools::TabixFile(vcfs[[tr]],index=paste0(vcfs[[tr]],'.tbi'))
    open(tb);on.exit(close(tb));parse_vcf(unique(unlist(Rsamtools::scanTabix(tb,param=gr),use.names=FALSE)))
  }
  # ---------------------- 01: candidate instruments ----------------------
  step01 <- function() {
    dest<-mkdir('01_candidate_snps');ans<-list()
    for(tr in c('pancreas_fat','pancreas_iron')) {
      cat('\nSTEP 01:',tr,'\n');con<-gzfile(vcfs[[tr]],'rt')
      res<-tryCatch({
        kept<-list();forms<-list();nr<-0;ns<-0;block<-0L
        repeat {
          ln<-readLines(con,n=100000L,warn=FALSE);if(!length(ln))break
          d<-parse_vcf(ln);if(!nrow(d))next
          block<-block+1L;nr<-nr+nrow(d)
          forms[[length(forms)+1L]]<-data.frame(FORMAT=names(table(d$FORMAT)),records=as.integer(table(d$FORMAT)))
          sig<-is.finite(d$LP)&d$LP > -log10(5e-8);ns<-ns+sum(sig)
          good<-sig&d$chr%in%as.character(1:22)&d$REF%in%c('A','C','G','T')&d$ALT%in%c('A','C','G','T')&d$REF!=d$ALT&valid_effect(d)
          if(any(good))kept[[length(kept)+1L]]<-as_effect(d[which(good)])
          if(block%%10L==0L){cat('Read:',nr,'\n');flush.console()}
        }
        if(!length(kept))stop('No candidate SNPs: ',tr)
        cand<-as.data.frame(data.table::rbindlist(kept))
        if(anyNA(cand$SNP)||anyDuplicated(cand$SNP)||anyNA(cand$pos))stop('Candidate IDs/positions need manual review; no silent deduplication.')
        cand<-cand[order(cand$pval,cand$SNP),,drop=FALSE]
        ft<-data.table::rbindlist(forms);ft<-ft[,list(records=sum(records)),by='FORMAT']
        list(candidates=cand,formats=ft,n_read=nr,n_significant=ns)
      },finally=close(con))
      wt(res$candidates,file.path(dest,paste0(ids[[tr]],'_candidates_before_clumping.tsv')))
      saveRDS(res,file.path(dest,paste0(ids[[tr]],'_extraction_result.rds')))
      compare(paste('01',tr),res$candidates,paste0(ids[[tr]],'_candidates_before_clumping.tsv'),'SNP',c('chr','pos','effect_allele','other_allele','beta','se','eaf','LP','pval'))
      check(paste('01 rows',tr),if(res$n_read==9263091) 'MATCH' else 'DIFFERENCE_REVIEW',res$n_read)
      ans[[tr]]<-res$candidates
    };ans
  }
  # ---------------------- 02: original PLINK clumping ----------------------
  step02 <- function(candidates) {
    dest<-mkdir('02_ld_clumping')
    bim<-data.table::fread(paste0(ref,'.bim'),header=FALSE,select=2L,showProgress=FALSE)[[1L]]
    nr<-nrow(data.table::fread(paste0(ref,'.fam'),header=FALSE,showProgress=FALSE))
    if(nr!=503L)stop('Expected 503 reference individuals.')
    saveRDS(list(prefix=ref,files=file.info(paste0(ref,c('.bed','.bim','.fam'))),reference_N=nr,variants=length(bim)),file.path(dest,'reference_manifest.rds'))
    ans<-list()
    for(tr in names(candidates)) {
      d<-candidates[[tr]];td<-file.path(dest,tr);dir.create(td)
      missing<-d[!d$SNP%in%bim,,drop=FALSE];wt(missing,file.path(td,'SNPs_missing_from_EUR.tsv'))
      input<-d[d$SNP%in%bim,,drop=FALSE]
      if(anyDuplicated(bim[bim%in%input$SNP]))stop('Reference has duplicate instrument IDs.')
      inf<-file.path(td,'clump_input.tsv');wt(data.frame(SNP=input$SNP,P=input$pval),inf)
      pref<-file.path(td,'clump')
      args<-c('--bfile',ref,'--clump',inf,'--clump-field','P','--clump-snp-field','SNP',
        '--clump-kb','10000','--clump-p1','5e-8','--clump-p2','1','--clump-r2','0.001',
        '--memory','2048','--threads','1','--out',pref)
      writeLines(paste(shQuote(c(plink,args)),collapse=' '),file.path(td,'command.txt'))
      rc<-system2(plink,vapply(args,shQuote,character(1)),stdout=file.path(td,'plink_console.txt'),stderr=file.path(td,'plink_console.txt'))
      if(rc!=0L||!file.exists(paste0(pref,'.clumped')))stop('PLINK failed: ',td)
      cl<-data.table::fread(paste0(pref,'.clumped'),fill=TRUE,showProgress=FALSE)
      if(!'SNP'%in%names(cl)||!nrow(cl)||anyNA(match(cl$SNP,d$SNP)))stop('Invalid clumping output.')
      instruments<-d[match(cl$SNP,d$SNP),,drop=FALSE]
      wt(instruments,file.path(td,paste0(ids[[tr]],'_EUR_clumped.tsv')))
      saveRDS(list(instruments=instruments,missing_from_reference=missing,clumped=cl,plink_args=args,reference_prefix=ref),file.path(td,'local_clumping_result.rds'))
      compare(paste('02',tr),instruments,paste0(ids[[tr]],'_EUR_clumped.tsv'),'SNP',c('beta','se','eaf','pval'))
      ans[[tr]]<-instruments
    };ans
  }
  # ---------------------- 03: outcome extraction, both directions ----------------------
  directions<-list(IPFD_to_pancreas_iron=c('pancreas_fat','pancreas_iron'),Pancreas_iron_to_IPFD=c('pancreas_iron','pancreas_fat'))
  step03 <- function(instruments) {
    dest<-mkdir('03_outcome_extraction');ans<-list()
    for(nm in names(directions)) {
      tr<-directions[[nm]];ex<-instruments[[tr[1L]]]
      gr<-GenomicRanges::GRanges(ex$chr,IRanges::IRanges(ex$pos,ex$pos))
      d<-read_query(tr[2L],gr);ix<-integer(nrow(ex))
      for(i in seq_len(nrow(ex))) {
        hit<-which(d$rsid==ex$SNP[i]&d$chr==ex$chr[i]&d$pos==ex$pos[i]&
          ((d$REF==ex$other_allele[i]&d$ALT==ex$effect_allele[i])|
           (d$REF==ex$effect_allele[i]&d$ALT==ex$other_allele[i])))
        if(length(hit)!=1L)stop('Missing/ambiguous coordinate+allele match: ',ex$SNP[i])
        ix[i]<-hit
      }
      d<-d[ix];oy<-as_effect(d);oy$FILTER<-d$FILTER;oy$valid_effect<-valid_effect(d)
      wt(oy,file.path(dest,paste0(nm,'_outcome.tsv')))
      ans[[nm]]<-list(exposure=ex,outcome=oy,traits=tr)
      saveRDS(ans[[nm]],file.path(dest,paste0(nm,'_extracted.rds')))
      if(nm=='IPFD_to_pancreas_iron') compare('03 forward',oy,'IPFD_instruments_in_pancreas_iron.tsv','SNP',c('beta','se','eaf','LP','valid_effect'))
    };ans
  }
  # ---------------------- 04: harmonisation ----------------------
  phenotype<-c(pancreas_fat='Pancreas fat',pancreas_iron='Pancreas iron',pancreas_volume='Pancreas volume')
  fmt_mr<-function(d,type,tr) {
    d<-as.data.frame(d);d$phenotype<-phenotype[[tr]];d$id<-ids[[tr]]
    TwoSampleMR::format_data(d,type=type,phenotype_col='phenotype',snp_col='SNP',beta_col='beta',se_col='se',
      eaf_col='eaf',effect_allele_col='effect_allele',other_allele_col='other_allele',pval_col='pval',
      chr_col='chr',pos_col='pos',id_col='id')
  }
  step04<-function(extracted) {
    dest<-mkdir('04_harmonisation');ans<-list()
    for(nm in names(extracted)) {
      x<-extracted[[nm]];oy<-x$outcome[!is.na(x$outcome$valid_effect)&x$outcome$valid_effect,,drop=FALSE]
      ex<-fmt_mr(x$exposure,'exposure',x$traits[1L]);oy<-fmt_mr(oy,'outcome',x$traits[2L])
      h<-TwoSampleMR::harmonise_data(ex,oy,action=2)
      keep<-h[!is.na(h$mr_keep)&h$mr_keep,,drop=FALSE]
      if(!nrow(keep))stop('No harmonised instruments: ',nm)
      wt(h,file.path(dest,paste0(nm,'_all.tsv')));wt(keep,file.path(dest,paste0(nm,'_kept.tsv')))
      saveRDS(list(exposure_dat=ex,outcome_dat=oy,harmonised=h,mr_dat=keep,action=2),file.path(dest,paste0(nm,'_harmonised.rds')))
      if(nm=='IPFD_to_pancreas_iron') compare('04 forward',h,'IPFD_to_pancreas_iron_all.tsv','SNP',c('beta.exposure','beta.outcome','mr_keep','palindromic','ambiguous'))
      if(nm=='Pancreas_iron_to_IPFD') compare('04 reverse',keep,'harmonised_single_SNP.tsv','SNP',c('beta.exposure','beta.outcome','mr_keep'))
      ans[[nm]]<-keep
    };ans
  }
  # ---------------------- 05: MR + nonpalindromic + reverse ----------------------
  step05<-function(harm) {
    dest<-mkdir('05_MR_results');pars<-TwoSampleMR::default_parameters();saveRDS(pars,file.path(dest,'MR_parameters.rds'))
    RNGkind('Mersenne-Twister','Inversion','Rejection')
    methods<-c('mr_ivw','mr_weighted_median','mr_egger_regression')
    main<-harm$IPFD_to_pancreas_iron
    dats<-list(main=main,nonpalindromic=main[!main$palindromic,,drop=FALSE],reverse=harm$Pancreas_iron_to_IPFD)
    for(nm in names(dats)) {
      d<-dats[[nm]];td<-file.path(dest,nm);dir.create(td)
      set.seed(seed) # sensitivity seed reconstructed; main seed recovered
      ml<-if(nm=='reverse') 'mr_wald_ratio' else methods
      if(nm=='reverse'&&nrow(d)!=1L)stop('Reverse instruments changed; historical Wald analysis no longer applicable.')
      if(nm!='reverse'&&nrow(d)<3L)stop('Too few SNPs for specified MR methods.')
      fit<-TwoSampleMR::mr(d,parameters=pars,method_list=ml)
      if(!is.data.frame(fit))stop('Unexpected TwoSampleMR API return type.')
      fit$ci_lower<-fit$b-1.96*fit$se;fit$ci_upper<-fit$b+1.96*fit$se
      wt(fit,file.path(td,'MR_estimates.tsv'));wt(d,file.path(td,'harmonised.tsv'))
      het<-pleio<-loo<-NULL
      if(nm!='reverse') {
        het<-TwoSampleMR::mr_heterogeneity(d,method_list=c('mr_ivw','mr_egger_regression'))
        pleio<-TwoSampleMR::mr_pleiotropy_test(d)
        loo<-TwoSampleMR::mr_leaveoneout(d,method=TwoSampleMR::mr_ivw)
        wt(het,file.path(td,'heterogeneity.tsv'));wt(pleio,file.path(td,'Egger_intercept.tsv'));wt(loo,file.path(td,'leave_one_out.tsv'))
      }
      saveRDS(list(data=d,mr_results=fit,heterogeneity=het,pleiotropy=pleio,leave_one_out=loo,
        seed=seed,parameters=pars,session_info=sessionInfo()),file.path(td,'MR_results.rds'))
      oldpath<-switch(nm,main='IPFD_to_iron_8adc24e844673/01_MR_estimates.tsv',
        nonpalindromic='IPFD_to_iron_nonpalindromic/MR_estimates.tsv',reverse='Pancreas_iron_to_IPFD/Wald_ratio_result.tsv')
      old<-baseline(oldpath)
      for(i in seq_len(nrow(fit))) {
        j<-match(fit$method[i],old$method);if(is.na(j))stop('Unmatched MR method.')
        for(cc in c('b','se','pval')) {
          delta<-abs(fit[[cc]][i]-old[[cc]][j]);pass<-is.finite(delta)&&delta<=1e-10+1e-6*abs(old[[cc]][j])
          stochastic<-fit$method[i]=='Weighted median'&&cc%in%c('se','pval')
          check(paste('05',nm,fit$method[i],cc),if(pass)'MATCH' else if(stochastic)'BOOTSTRAP_DIFFERENCE_REVIEW' else 'DIFFERENCE_REVIEW',delta)
        }
      }
      if(nm=='main') {
        compare('05 heterogeneity',het,'IPFD_to_iron_8adc24e844673/02_heterogeneity.tsv','method',c('Q','Q_df','Q_pval'))
        compare('05 Egger intercept',pleio,'IPFD_to_iron_8adc24e844673/03_Egger_intercept.tsv',c('id.exposure','id.outcome'),c('egger_intercept','se','pval'))
        compare('05 leave-one-out',loo,'IPFD_to_iron_8adc24e844673/04_leave_one_out.tsv','SNP',c('b','se','p'))
      }
      grDevices::pdf(file.path(td,'MR_estimates.pdf'),width=8,height=4)
      yy<-seq_len(nrow(fit));graphics::par(mar=c(5,11,2,1))
      graphics::plot(fit$b,yy,xlim=range(c(fit$ci_lower,fit$ci_upper,0)),yaxt='n',ylab='',xlab='Beta (GWAS analysis scale)',pch=19)
      graphics::axis(2,at=yy,labels=fit$method,las=1);graphics::segments(fit$ci_lower,yy,fit$ci_upper,yy);graphics::abline(v=0,lty=2);grDevices::dev.off()
    }
  }
  # ---------------------- Dense regional ABF helpers ----------------------
  get_region<-function(tr,dest) {
    gr<-GenomicRanges::GRanges('8',IRanges::IRanges(region_start,region_end))
    d<-read_query(tr,gr)
    good<-d$chr=='8'&is.finite(d$pos)&d$pos>=region_start&d$pos<=region_end&
      grepl('^[ACGT]+$',d$REF)&grepl('^[ACGT]+$',d$ALT)&d$REF!=d$ALT&valid_effect(d)&is.finite(d$AF)&d$AF>0&d$AF<1
    good[is.na(good)]<-FALSE;wt(d[!good],file.path(dest,paste0(ids[[tr]],'_excluded.tsv')))
    d<-d[good];if(anyDuplicated(d$key))stop('Duplicate regional key: ',tr)
    wt(d,file.path(dest,paste0(ids[[tr]],'_region.tsv')))
    # Dense region includes valid indels and palindromic variants, no clumping/P filtering.
    compare(paste('region',tr),d,paste0(ids[[tr]],'_region.tsv'),'key',c('pos','REF','ALT','ES','SE','LP','AF'))
    as.data.frame(d)
  }
  coloc_ds<-function(d,tr) list(snp=d$key,position=d$pos,beta=d$ES,varbeta=d$SE^2,
    MAF=pmin(d$AF,1-d$AF),N=Ns[[tr]],type='quant')
  run_coloc<-function(a,b,tr1,tr2,label,dest) {
    common<-intersect(a$key,b$key);if(length(common)<100L)stop('Insufficient regional variants.')
    a<-a[match(common,a$key),,drop=FALSE];b<-b[match(common,b$key),,drop=FALSE]
    stopifnot(identical(a$key,b$key),identical(a$pos,b$pos))
    ds1<-coloc_ds(a,tr1);ds2<-coloc_ds(b,tr2)
    rows<-list();fits<-list()
    for(pr in c(1e-6,1e-5,5e-5)) {
      fit<-coloc::coloc.abf(ds1,ds2,p1=1e-4,p2=1e-4,p12=pr)
      s<-fit$summary
      rows[[length(rows)+1L]]<-data.frame(p1=1e-4,p2=1e-4,p12=pr,nsnps=unname(s['nsnps']),
        PP_H0=unname(s['PP.H0.abf']),PP_H1=unname(s['PP.H1.abf']),PP_H2=unname(s['PP.H2.abf']),PP_H3=unname(s['PP.H3.abf']),PP_H4=unname(s['PP.H4.abf']))
      fits[[as.character(pr)]]<-fit
      if(pr==1e-5) {
        v<-fit$results;v$rsid<-a$rsid[match(v$snp,a$key)];v<-v[order(v$SNP.PP.H4,decreasing=TRUE),];v$cumulative_PP_H4<-cumsum(v$SNP.PP.H4)
        wt(v,file.path(dest,paste0(label,'_variant_posteriors.tsv')))
      }
    }
    summary<-do.call(rbind,rows)
    wt(summary,file.path(dest,paste0(label,'_prior_sensitivity.tsv')))
    saveRDS(list(dataset1=ds1,dataset2=ds2,table1=a,table2=b,fits=fits),file.path(dest,paste0(label,'_coloc_results.rds')))
    list(summary=summary,fits=fits)
  }
  # ---------------------- 06: fat-iron ABF ----------------------
  step06<-function() {
    dest<-mkdir('06_colocalisation');a<-get_region('pancreas_fat',dest);b<-get_region('pancreas_iron',dest)
    fit<-run_coloc(a,b,'pancreas_fat','pancreas_iron','fat_iron',dest)
    compare('06 ABF',fit$summary,'coloc_prior_sensitivity.tsv','p12',c('nsnps','PP_H0','PP_H1','PP_H2','PP_H3','PP_H4'))
    v<-as.data.frame(data.table::fread(file.path(dest,'fat_iron_variant_posteriors.tsv')))
    compare('06 variant posterior',v,'coloc_variant_results.tsv','snp','SNP.PP.H4')
    list(pancreas_fat=a,pancreas_iron=b)
  }
  # ---------------------- 08: volume + regional plots ----------------------
  step08<-function(regional) {
    dest<-mkdir('08_pancreas_volume_coloc');regional$pancreas_volume<-get_region('pancreas_volume',dest)
    for(tr in names(regional))wt(regional[[tr]],file.path(dest,paste0(ids[[tr]],'_region.tsv')))
    common<-Reduce(intersect,lapply(regional,function(d)d$key))
    studies<-data.frame(trait=names(ids),id=unname(ids),expected_N=unname(Ns),file=unname(vcfs))
    saveRDS(list(studies=studies,metadata=metadata,datasets=regional,common_three=common,
      build='GRCh37',chromosome='8',start=region_start,end=region_end),file.path(dest,'three_trait_region_inputs.rds'))
    rows<-list();fits<-list()
    for(tr in c('pancreas_fat','pancreas_iron')) {
      label<-if(tr=='pancreas_fat')'fat_volume' else 'iron_volume'
      r<-run_coloc(regional[[tr]],regional$pancreas_volume,tr,'pancreas_volume',label,dest)
      s<-r$summary;s$comparison<-label;s$trait1<-tr;s$trait2<-'pancreas_volume';s$N1<-Ns[[tr]];s$N2<-Ns[['pancreas_volume']]
      rows[[label]]<-s;fits[[label]]<-r$fits
      v<-as.data.frame(data.table::fread(file.path(dest,paste0(label,'_variant_posteriors.tsv'))))
      compare(paste('08 posterior',label),v,paste0(label,'_variant_posteriors.tsv'),'snp','SNP.PP.H4')
    }
    s<-data.table::rbindlist(rows);wt(s,file.path(dest,'pairwise_coloc_summary.tsv'));saveRDS(fits,file.path(dest,'pairwise_coloc_results.rds'))
    compare('08 ABF',s,'pairwise_coloc_summary.tsv',c('comparison','p12'),c('nsnps','N1','N2','PP_H0','PP_H1','PP_H2','PP_H3','PP_H4'))
    grDevices::pdf(file.path(dest,'three_trait_region_plot.pdf'),width=9,height=9)
    graphics::par(mfrow=c(3,1),mar=c(4,4,2,1))
    for(tr in names(regional)) {
      d<-regional[[tr]];graphics::plot(d$pos/1e6,d$LP,pch=16,cex=.45,col='#32688A',
        xlim=c(region_start,region_end)/1e6,xlab='Chromosome 8 position (Mb, GRCh37)',ylab='-log10(P)',main=phenotype[[tr]])
      graphics::abline(h=-log10(5e-8),lty=2,col='grey50')
      graphics::abline(v=129547378/1e6,lty=3,col='#A54B3C')
    }
    grDevices::dev.off()
    cat('Step08 input for future use:',file.path(dest,'three_trait_region_inputs.rds'),'\n')
  }
  # ---------------------- Execute in dependency order ----------------------
  withCallingHandlers({
    cand<-step01();inst<-step02(cand);extracted<-step03(inst);harm<-step04(extracted)
    step05(harm);regional<-step06();step08(regional)
    complete<-TRUE
  },warning=function(w) {
    warning_messages<<-c(warning_messages,conditionMessage(w))
    cat('WARNING:',conditionMessage(w),'\n');invokeRestart('muffleWarning')
  })
  tab<-data.table::rbindlist(checks,fill=TRUE)
  cat('\n===== Reconstruction finished =====\n')
  print(table(tab$status));cat('Output:',out,'\nUpload:',file.path(out,'reconstruction_review.zip'),'\n')
  cat('Historical comparison differences must be reviewed before submission.\n')
  cat('No original results were overwritten. No causal conclusion is added by reconstruction.\n')
})()
