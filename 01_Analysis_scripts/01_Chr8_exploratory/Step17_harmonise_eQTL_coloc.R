# 第17步：区域 GWAS 与胰腺/胰岛 eQTL 的坐标协调、覆盖审计、探索性 ABF 共定位。
# 在远程 RStudio 中 source() 运行整个文件；无需逐段替换。
# 输入：原始 GRCh37 VCF 和第16步结果；输出：新建17_QTL_coloc子文件夹。
# 自动下载 UCSC Linux x86_64 liftOver 到个人 mr_tools 目录，无需管理员权限。
# 当前分析限定双等位 SNV（保留回文SNV，按基因组正链REF/ALT精确匹配）。
# indel 未进行跨版本规范化，明确排除并计数；不按P值或LD剪枝筛选。
# ABF每个性状每个区域至多一个因果变异；不是QTL SuSiE，也不是机制验证。
(function() {
  need <- c('data.table','Rsamtools','GenomicRanges','IRanges','coloc')
  available <- vapply(need,requireNamespace,logical(1),quietly=TRUE)
  if(!all(available)) stop('Missing R packages: ',paste(need[!available],collapse=', '))
  base <- path.expand('~/胰腺')
  qdir <- file.path(base,'16_regional_eQTL','chr8_eQTL_ff7a81b8b3096')
  ids <- c(pancreas_fat='ebi-a-GCST90016675',pancreas_iron='ebi-a-GCST90016676',
           pancreas_volume='ebi-a-GCST90016669')
  N_gwas <- c(pancreas_fat=25617,pancreas_iron=25617,pancreas_volume=31758)
  dids <- c('QTD000296','QTD000554','QTD000574')
  qfiles <- file.path(qdir,paste0(dids,'_region_annotated.tsv.gz'))
  meta_file <- file.path(qdir,'selected_dataset_metadata.tsv')
  vfiles <- file.path(base,paste0(ids,'.vcf.gz'))
  required_files <- c(qfiles,meta_file,vfiles,paste0(vfiles,'.tbi'))
  if(any(!file.exists(required_files))) stop('Missing input files:\n',paste(required_files[!file.exists(required_files)],collapse='\n'))
  parent <- file.path(base,'17_QTL_coloc')
  dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile('chr8_QTL_coloc_',tmpdir=parent)
  if(!dir.create(out)) stop('Cannot create output directory')
  for(nm in c('mapping','clean_QTL','matched_inputs','coloc_fits')) dir.create(file.path(out,nm))
  sink(file.path(out,'run_log.txt'),split=TRUE);on.exit(sink(),add=TRUE)
  old_options <- options(timeout=600);on.exit(options(old_options),add=TRUE)
  warnings_seen <- character(); completed <- FALSE
  on.exit({
    writeLines(if(length(warnings_seen)) warnings_seen else 'No captured warnings/errors.',file.path(out,'warnings_and_errors.txt'))
    writeLines(if(completed) 'COMPLETED' else 'INCOMPLETE: inspect run_log and error message.',file.path(out,'completion_status.txt'))
    writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  },add=TRUE)
  wt <- function(x,name) data.table::fwrite(x,file.path(out,name),sep='\t',na='NA')
  read_tsv <- function(f) as.data.frame(data.table::fread(f,showProgress=FALSE))
  # 使用连接读取gzip，避免依赖额外的R.utils包。
  read_gz_tsv <- function(f) {
    con <- gzfile(f,'rt');on.exit(close(con))
    as.data.frame(data.table::fread(text=paste(readLines(con,warn=FALSE),collapse='\n'),
                                  sep='\t',showProgress=FALSE))
  }
  cat('Output:',out,'\n')
  script_source <- file.path(base,'17_harmonise_eQTL_coloc.R')
  if(file.exists(script_source)) file.copy(script_source,file.path(out,'17_harmonise_eQTL_coloc.R'),overwrite=FALSE)
  meta <- read_tsv(meta_file)
  if(anyDuplicated(meta$dataset_id)) stop('Duplicate dataset metadata')
  meta <- meta[match(dids,meta$dataset_id),,drop=FALSE]
  if(anyNA(meta$dataset_id)||any(meta$quant_method!='ge')) stop('Invalid ge metadata')
  meta$sample_size <- as.numeric(meta$sample_size)
  if(any(!is.finite(meta$sample_size)|meta$sample_size<=0)) stop('Invalid QTL sample size')
  wt(meta,'selected_dataset_metadata.tsv')
  core <- data.frame(SNP=c('rs4733612','rs1561927','rs7824269','rs7819154','rs12676584'),
     pos37=c(129569999L,129568078L,129552344L,129558351L,129547378L),
     pos38=c(128557753L,128555832L,128540098L,128546105L,128535132L),
     REF=c('G','C','C','C','G'),ALT=c('A','T','T','A','A'))
  core$key37 <- with(core,paste('8',pos37,REF,ALT,sep=':'))
  core$key38 <- with(core,paste('8',pos38,REF,ALT,sep=':'))
  wt(core,'five_candidates_manifest.tsv')
  start37 <- 129047378L;end37 <- 130047378L
  warnings_seen <- c(warnings_seen,
    'Study-level N is used for sdY estimation; variant/gene-specific effective N is not confirmed.',
    'SNV-only analysis: indels are excluded, not assumed irrelevant.',
    'ABF assumes at most one causal variant per trait in the analysed interval; sample overlap is not modelled.')

  main <- function() {
    # 1. 重新读取小区域，按每条FORMAT解析，不依赖之前RDS内部命名。
    gwas_qc <- list();raw_gwas <- list()
    for(tr in names(ids)) {
      f <- file.path(base,paste0(ids[[tr]],'.vcf.gz'))
      con <- gzfile(f,'rt');h <- tryCatch(readLines(con,n=250L,warn=FALSE),finally=close(con))
      sh <- h[startsWith(h,'##SAMPLE=')]
      if(length(sh)!=1L || !grepl(ids[[tr]],sh,fixed=TRUE) ||
         !any(grepl('assembly=HG19/GRCh37',h,fixed=TRUE))) stop('Unexpected VCF metadata: ',tr)
      n_hdr <- as.numeric(sub('.*TotalControls=([0-9.]+).*','\\1',sh))
      if(!is.finite(n_hdr)||n_hdr!=N_gwas[[tr]]) stop('Sample size differs: ',tr)
      tb <- Rsamtools::TabixFile(f,index=paste0(f,'.tbi'))
      open(tb)
      lines <- tryCatch(Rsamtools::scanTabix(tb,param=GenomicRanges::GRanges('8',
        IRanges::IRanges(start37,end37)))[[1L]],finally=close(tb))
      if(!length(lines)) stop('Empty GWAS region: ',tr)
      d <- as.data.frame(data.table::fread(text=paste(lines,collapse='\n'),sep='\t',header=FALSE,
        colClasses='character',na.strings='.',quote='',showProgress=FALSE))
      if(ncol(d)!=10L) stop('Expected ten-column GWAS VCF')
      names(d) <- c('chr','pos','rsid','REF','ALT','QUAL','FILTER','INFO','FORMAT','VALUE')
      d$pos <- as.integer(d$pos)
      for(v in c('ES','SE','LP','AF')) d[[v]] <- NA_real_
      if(anyNA(d$FORMAT)) stop('Missing FORMAT')
      for(fmt in unique(d$FORMAT)) {
        ix <- which(d$FORMAT==fmt);ks <- strsplit(fmt,':',fixed=TRUE)[[1L]]
        if(anyDuplicated(ks)||!all(c('ES','SE','LP','AF') %in% ks)) stop('Unsupported FORMAT: ',fmt)
        vs <- strsplit(d$VALUE[ix],':',fixed=TRUE)
        if(any(lengths(vs)!=length(ks))) stop('Invalid sample field length')
        for(v in c('ES','SE','LP','AF')) d[[v]][ix] <- suppressWarnings(as.numeric(vapply(vs,`[`,character(1),match(v,ks))))
      }
      snv <- d$REF %in% c('A','C','G','T') & d$ALT %in% c('A','C','G','T') & d$REF!=d$ALT
      valid <- snv & d$chr=='8' & is.finite(d$pos) & d$pos>=start37 & d$pos<=end37 &
        is.finite(d$ES) & is.finite(d$SE) & d$SE>0 & is.finite(d$LP) & d$LP>=0 &
        is.finite(d$AF) & d$AF>0 & d$AF<1 & (is.na(d$FILTER)|d$FILTER=='PASS')
      valid[is.na(valid)] <- FALSE
      d$key37 <- with(d,paste(chr,pos,REF,ALT,sep=':'))
      wt(d[!valid,],paste0('mapping/',tr,'_excluded_before_mapping.tsv'))
      gwas_qc[[tr]] <- data.frame(trait=tr,raw_rows=nrow(d),non_SNV=sum(!snv,na.rm=TRUE),valid_SNV=sum(valid),N=n_hdr)
      d <- d[valid,,drop=FALSE]
      if(anyDuplicated(d$key37)) stop('Duplicate GWAS keys: ',tr)
      if(!all(core$key37 %in% d$key37)) stop('Core candidate missing in GWAS: ',tr)
      raw_gwas[[tr]] <- d
      cat(tr,': retained=',nrow(d),' valid SNVs; N=',n_hdr,'\n')
    }
    wt(do.call(rbind,gwas_qc),'GWAS_input_QC.tsv')

    # 2. 官方UCSC工具执行逐SNV liftOver；-multiple用于检测非唯一映射。
    tools_dir <- path.expand('~/mr_tools');dir.create(tools_dir,recursive=TRUE,showWarnings=FALSE)
    bin <- file.path(tools_dir,'liftOver');chain <- file.path(tools_dir,'hg19ToHg38.over.chain.gz')
    bin_url <- 'https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver'
    chain_url <- 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
    fetch <- function(u,p) {
      if(!file.exists(p)||file.info(p)$size==0) {
        tmp <- paste0(p,'.part');download.file(u,tmp,method='libcurl',mode='wb',quiet=FALSE)
        if(!file.exists(tmp)||file.info(tmp)$size<100L) stop('Download incomplete: ',u)
        if(!file.rename(tmp,p)) stop('Cannot save download: ',p)
      }
    }
    if(Sys.info()[['sysname']]!='Linux'||!Sys.info()[['machine']] %in% c('x86_64','amd64'))
      stop('This script downloads the Linux x86_64 UCSC executable; incompatible server architecture.')
    fetch(bin_url,bin);Sys.chmod(bin,'0700');fetch(chain_url,chain)
    probe <- suppressWarnings(system2(bin,stdout=TRUE,stderr=TRUE))
    if(!any(grepl('liftOver',probe,fixed=TRUE))||!any(grepl('usage:',probe,fixed=TRUE)))
      stop('UCSC liftOver cannot run on this server: ',paste(probe,collapse=' '))
    writeLines(probe,file.path(out,'mapping','liftOver_version_usage.txt'))
    u <- unique(do.call(rbind,lapply(raw_gwas,function(d) d[,c('key37','pos','REF','ALT')])))
    u <- u[order(u$pos,u$key37),];u$map_id <- sprintf('v%06d',seq_len(nrow(u)))
    bed <- data.frame(chr='chr8',start=u$pos-1L,end=u$pos,name=u$map_id,score=0,strand='+')
    oldbed <- file.path(out,'mapping','GWAS_GRCh37.bed')
    newbed <- file.path(out,'mapping','GWAS_GRCh38.bed')
    badbed <- file.path(out,'mapping','GWAS_unmapped.bed')
    data.table::fwrite(bed,oldbed,sep='\t',col.names=FALSE)
    msg <- system2(bin,args=c('-multiple','-noSerial',shQuote(oldbed),shQuote(chain),shQuote(newbed),shQuote(badbed)),stdout=TRUE,stderr=TRUE)
    writeLines(msg,file.path(out,'mapping','liftOver_log.txt'))
    rc <- attr(msg,'status');if((!is.null(rc)&&rc!=0L)||!file.exists(newbed)||file.info(newbed)$size==0) stop('liftOver failed')
    b <- as.data.frame(data.table::fread(newbed,header=FALSE))
    if(ncol(b)!=6L) stop('Unexpected liftOver BED structure')
    names(b) <- c('chr38','start38','end38','map_id','score','strand38')
    counts <- table(b$map_id)
    u$n_mappings <- as.integer(counts[match(u$map_id,names(counts))]);u$n_mappings[is.na(u$n_mappings)] <- 0L
    j <- match(u$map_id,b$map_id)
    u$chr38 <- b$chr38[j];u$pos38 <- b$start38[j]+1L;u$strand38 <- b$strand38[j]
    u$mapping_status <- ifelse(u$n_mappings==0,'UNMAPPED',ifelse(u$n_mappings>1,'MULTIPLE',
      ifelse(u$chr38!='chr8'|b$end38[j]-b$start38[j]!=1L,'NON_PRIMARY_OR_NON_UNIT','UNIQUE_CHR8')))
    comp <- function(x) chartr('ACGT','TGCA',x)
    u$REF38 <- ifelse(u$strand38=='-',comp(u$REF),u$REF)
    u$ALT38 <- ifelse(u$strand38=='-',comp(u$ALT),u$ALT)
    u$key38 <- paste('8',u$pos38,u$REF38,u$ALT38,sep=':')
    collide <- u$mapping_status=='UNIQUE_CHR8' & (duplicated(u$key38)|duplicated(u$key38,fromLast=TRUE))
    u$mapping_status[collide] <- 'TARGET_KEY_COLLISION'
    wt(u,'mapping/GWAS_liftOver_audit.tsv')
    test <- u[match(core$key37,u$key37),]
    if(anyNA(test$key37)||any(test$mapping_status!='UNIQUE_CHR8')||!identical(as.character(test$key38),as.character(core$key38)))
      stop('Five-candidate coordinate verification failed; do not proceed.')
    wt(data.frame(core,mapped_key38=test$key38,verified=TRUE),'mapping/five_candidate_coordinate_check.tsv')
    wt(data.frame(resource=c('UCSC liftOver','hg19ToHg38 chain'),url=c(bin_url,chain_url),
      local_path=c(bin,chain),md5=unname(tools::md5sum(c(bin,chain)))),'mapping/resource_manifest.tsv')
    good <- u[u$mapping_status=='UNIQUE_CHR8',]
    gwas <- lapply(raw_gwas,function(d) {
      a <- merge(d,good[,c('key37','pos38','REF38','ALT38','key38','strand38')],by='key37',sort=FALSE)
      # 负链映射仅互补碱基，不交换效应/非效应等位基因，因此beta和AF不变。
      a <- a[order(a$pos38,a$key38),];a$MAF <- pmin(a$AF,1-a$AF);a
    })
    for(tr in names(gwas)) wt(gwas[[tr]],paste0('mapping/',tr,'_GRCh38_SNV.tsv'))
    cat('Mapping passed: all five core candidates verified.\n')

    # 3. QTL去重与QC；同键统计量冲突即停止，不选择最小P值。
    qtls <- list();qc <- list();coverage <- list()
    for(i in seq_along(dids)) {
      id <- dids[i];d <- read_gz_tsv(qfiles[i])
      req <- c('molecular_trait_id','gene_id','chromosome','position','ref','alt','beta','se','pvalue','maf','rsid')
      if(!all(req %in% names(d))) stop('QTL required columns missing: ',id)
      for(nm in c('position','beta','se','pvalue','maf')) d[[nm]] <- as.numeric(d[[nm]])
      d$key38 <- paste(sub('^chr','',d$chromosome),d$position,d$ref,d$alt,sep=':')
      d$gene <- sub('\\..*$','',d$gene_id)
      d$pair_key <- paste(d$molecular_trait_id,d$key38,sep='|')
      dup <- duplicated(d$pair_key)|duplicated(d$pair_key,fromLast=TRUE)
      wt(d[dup,],paste0('clean_QTL/',id,'_duplicate_rows.tsv'))
      groups <- split(seq_len(nrow(d)),d$pair_key)
      stat_cols <- intersect(c('beta','se','pvalue','maf','gene_id','an','ac','r2','median_tpm'),names(d))
      conflicts <- vapply(groups,function(ix) nrow(unique(d[ix,stat_cols,drop=FALSE]))>1L,logical(1))
      if(any(conflicts)) stop('Conflicting duplicate QTL statistics: ',id,'; inspect duplicate_rows.tsv')
      alias <- vapply(groups,function(ix) paste(sort(unique(na.omit(as.character(d$rsid[ix])))),collapse=';'),character(1))
      raw_n <- nrow(d);d <- d[!duplicated(d$pair_key),];d$rsid_aliases <- alias[d$pair_key]
      unique_n <- nrow(d)
      snv <- d$ref %in% c('A','C','G','T') & d$alt %in% c('A','C','G','T') & d$ref!=d$alt
      valid <- snv & sub('^chr','',d$chromosome)=='8' & is.finite(d$position) &
        is.finite(d$beta) & is.finite(d$se) & d$se>0 & is.finite(d$pvalue) & d$pvalue>=0 & d$pvalue<=1 &
        is.finite(d$maf) & d$maf>0 & d$maf<=0.5 & !is.na(d$molecular_trait_id) & !is.na(d$gene)
      valid[is.na(valid)] <- FALSE
      wt(d[!valid,],paste0('clean_QTL/',id,'_excluded_unique_rows.tsv'))
      qc[[id]] <- data.frame(dataset_id=id,raw_rows=raw_n,unique_pairs=unique_n,
        collapsed_extra_rows=raw_n-unique_n,excluded_unique_rows=sum(!valid),retained_SNV_pairs=sum(valid))
      d <- d[valid,];qtls[[id]] <- d
      wt(d,paste0('clean_QTL/',id,'_clean_SNV.tsv'))
      for(g in unique(d$gene)) {
        a <- d[d$gene==g,]
        coverage[[paste(id,g)]] <- data.frame(dataset_id=id,gene_id=g,records=nrow(a),
          core_variants=sum(core$key38 %in% a$key38),minimum_nominal_p=min(a$pvalue),
          min_pos38=min(a$position),max_pos38=max(a$position))
      }
    }
    wt(do.call(rbind,qc),'QTL_deduplication_QC.tsv')
    wt(do.call(rbind,coverage),'all_gene_coverage_after_QC.tsv')

    # 4. 预先指定三个候选基因。无记录/无核心位点覆盖的组合如实记为不可评估。
    genes <- c(MYC='ENSG00000136997',PVT1='ENSG00000249859',LINC00824='ENSG00000254275')
    availability <- list();pair_qc <- list();summaries <- list();status_rows <- list()
    priors <- c(1e-6,1e-5,5e-5)
    pdffile <- file.path(out,'regional_GWAS_eQTL_comparison.pdf')
    grDevices::pdf(pdffile,width=10,height=11,onefile=TRUE)
    pdf_device <- grDevices::dev.cur();on.exit({if(pdf_device %in% grDevices::dev.list()) grDevices::dev.off(pdf_device)},add=TRUE)
    panel <- function(pos,lp,keys,label,xlim) {
      plot(pos/1e6,lp,pch=16,cex=0.45,col='#667788',xlim=xlim,
        xlab='GRCh38 chromosome 8 (Mb)',ylab='-log10(P)',main=label)
      h <- keys %in% core$key38;points(pos[h]/1e6,lp[h],pch=16,col='#C23B22',cex=0.9)
      abline(v=core$pos38/1e6,col='#E4C5BE',lty=3)
    }
    for(id in dids) for(symbol in names(genes)) {
      gene <- genes[[symbol]];qd <- qtls[[id]];qd <- qd[qd$gene==gene,,drop=FALSE]
      nc <- sum(core$key38 %in% qd$key38)
      state <- if(!nrow(qd)) 'NO_REGIONAL_RECORD' else if(nc!=5L) 'INCOMPLETE_CORE_COVERAGE' else 'EVALUABLE'
      availability[[paste(id,symbol)]] <- data.frame(dataset_id=id,gene_symbol=symbol,gene_id=gene,
        retained_SNV_pairs=nrow(qd),core_variants=nc,status=state)
      if(state!='EVALUABLE') next
      if(length(unique(qd$molecular_trait_id))!=1L) stop('More than one molecular trait for gene: ',id,' ',symbol)
      if(anyDuplicated(qd$key38)) stop('QTL key duplication remained')
      cat('\n===== ',id,' / ',symbol,' =====\n',sep='')
      par(mfrow=c(4,1),mar=c(3.8,4.2,2.2,1),oma=c(0,0,3,0))
      xr <- range(unlist(lapply(gwas,function(d)d$pos38)))/1e6
      for(tr in names(gwas)) {
        a <- gwas[[tr]];panel(a$pos38,a$LP,a$key38,tr,xr)
      }
      within <- qd$position>=xr[1]*1e6 & qd$position<=xr[2]*1e6
      panel(qd$position[within],-log10(pmax(qd$pvalue[within],1e-300)),qd$key38[within],paste(id,symbol,'expression'),xr)
      mtext(paste('Regional association comparison:',id,symbol,'| red = five core variants'),outer=TRUE,cex=1)
      for(tr in names(gwas)) {
        a <- gwas[[tr]];ix <- match(a$key38,qd$key38);found <- !is.na(ix)
        x <- a[found,,drop=FALSE];y <- qd[ix[found],,drop=FALSE]
        stopifnot(identical(x$key38,y$key38))
        nshared <- nrow(x);ncore <- sum(core$key38 %in% x$key38)
        label <- paste(id,symbol,tr,sep='__')
        # 保存每一条GWAS位点的匹配状态，REF/ALT反向时也不静默翻转。
        audit <- data.frame(SNP=a$rsid,key38=a$key38,status=ifelse(found,'EXACT_REF_ALT_MATCH',
          ifelse(a$pos38 %in% qd$position,'POSITION_PRESENT_ALLELE_MISMATCH','NO_QTL_RECORD_AT_POSITION')))
        wt(audit,paste0('matched_inputs/',label,'_matching_audit.tsv'))
        m <- data.frame(key38=x$key38,pos38=x$pos38,GWAS_rsid=x$rsid,QTL_rsid_aliases=y$rsid_aliases,
          REF=x$REF38,ALT=x$ALT38,GWAS_beta=x$ES,GWAS_se=x$SE,GWAS_LP=x$LP,GWAS_MAF=x$MAF,
          QTL_beta=y$beta,QTL_se=y$se,QTL_p=y$pvalue,QTL_MAF=y$maf,
          core_candidate=x$key38 %in% core$key38)
        wt(m,paste0('matched_inputs/',label,'.tsv'))
        eligible <- nshared>=100L && ncore==5L
        pair_qc[[label]] <- data.frame(dataset_id=id,gene_symbol=symbol,GWAS_trait=tr,
          mapped_GWAS_SNV=nrow(a),QTL_gene_SNV=nrow(qd),shared_SNV=nshared,
          fraction_GWAS_matched=nshared/nrow(a),core_variants=ncore,
          min_QTL_p=if(nshared) min(y$pvalue) else NA_real_,
          min_GWAS_p=if(nshared) 10^(-max(x$LP)) else NA_real_,
          max_MAF_difference=if(nshared) max(abs(x$MAF-y$maf)) else NA_real_,
          status=if(eligible) 'EXPLORATORY_ABF' else 'SKIPPED_INSUFFICIENT_COVERAGE')
        if(!eligible) next
        nq <- meta$sample_size[match(id,meta$dataset_id)]
        for(scale_mode in c('sdY_estimated','sdY1_sensitivity')) for(p12 in priors) {
          run <- paste(label,scale_mode,format(p12,scientific=TRUE),sep='__')
          d1 <- list(snp=x$key38,position=x$pos38,beta=x$ES,varbeta=x$SE^2,MAF=x$MAF,N=unname(N_gwas[[tr]]),type='quant')
          d2 <- list(snp=y$key38,position=y$position,beta=y$beta,varbeta=y$se^2,MAF=y$maf,N=nq,type='quant')
          if(scale_mode=='sdY1_sensitivity') {d1$sdY <- 1;d2$sdY <- 1}
          cat('Running:',run,' n=',nshared,'\n');flush.console()
          tryCatch({
            fit <- withCallingHandlers(coloc::coloc.abf(d1,d2,p1=1e-4,p2=1e-4,p12=p12),warning=function(w) {
              warnings_seen <<- c(warnings_seen,paste(run,conditionMessage(w),sep=': '));invokeRestart('muffleWarning')
            })
            s <- fit$summary
            row <- data.frame(dataset_id=id,gene_symbol=symbol,GWAS_trait=tr,scale_mode=scale_mode,
              GWAS_N=N_gwas[[tr]],QTL_N=nq,nsnps=unname(s['nsnps']),p1=1e-4,p2=1e-4,p12=p12,
              PP_H0=unname(s['PP.H0.abf']),PP_H1=unname(s['PP.H1.abf']),PP_H2=unname(s['PP.H2.abf']),
              PP_H3=unname(s['PP.H3.abf']),PP_H4=unname(s['PP.H4.abf']),minimum_QTL_p=min(y$pvalue))
            row$PP_QTL_association <- row$PP_H2+row$PP_H3+row$PP_H4
            row$H4_fraction_H3_H4 <- if((row$PP_H3+row$PP_H4)>0) row$PP_H4/(row$PP_H3+row$PP_H4) else NA_real_
            summaries[[run]] <- row
            saveRDS(list(input_GWAS=d1,input_QTL=d2,fit=fit,scale_mode=scale_mode),file.path(out,'coloc_fits',paste0(run,'.rds')))
            if(scale_mode=='sdY_estimated'&&p12==1e-5) {
              vr <- fit$results;vr <- vr[order(vr$SNP.PP.H4,decreasing=TRUE),]
              vr$GWAS_rsid <- x$rsid[match(vr$snp,x$key38)]
              wt(vr,paste0('coloc_fits/',label,'_variant_posterior_conditional_on_H4.tsv'))
            }
            status_rows[[run]] <- data.frame(run=run,status='COMPLETED',message='')
          },error=function(e) {
            warnings_seen <<- c(warnings_seen,paste(run,'ERROR',conditionMessage(e)))
            status_rows[[run]] <<- data.frame(run=run,status='ERROR',message=conditionMessage(e))
          })
        }
      }
    }
    grDevices::dev.off(pdf_device)
    wt(do.call(rbind,availability),'priority_gene_availability.tsv')
    if(length(pair_qc)) wt(do.call(rbind,pair_qc),'pair_coverage_and_signal_QC.tsv')
    if(length(status_rows)) wt(do.call(rbind,status_rows),'run_status.tsv')
    if(length(summaries)) {
      res <- do.call(rbind,summaries);rownames(res) <- NULL
      wt(res,'QTL_coloc_all_settings.tsv')
      primary <- res[res$scale_mode=='sdY_estimated' & res$p12==1e-5,]
      wt(primary,'QTL_coloc_primary_summary.tsv')
      cat('\n===== Primary exploratory results =====\n');print(primary,row.names=FALSE)
    } else stop('No completed coloc runs; inspect availability and errors.')
    saveRDS(list(input_region_GRCh37=c(start37,end37),QTL_input=qdir,GWAS_ids=ids,
      candidates=core,priors=priors,minimum_shared_SNV=100L,required_core=5L,
      method='SNV-only coloc.abf; single causal variant assumption',
      scale_modes=c('sdY_estimated','sdY1_sensitivity'),sample_size_source='study-level metadata'),
      file.path(out,'analysis_config.rds'))
    wt(data.frame(file=c(qfiles,meta_file),md5=unname(tools::md5sum(c(qfiles,meta_file)))),'QTL_input_checksums.tsv')
    writeLines(c(
      'Exploratory regional eQTL colocalisation audit, not proof of molecular mediation.',
      'Dataset1 = imaging GWAS; dataset2 = gene-expression QTL.',
      'H0 neither associated; H1 GWAS only; H2 QTL only; H3 both with distinct causal variants; H4 shared causal variant.',
      'Low H4 with high H1 means insufficient QTL association evidence, not demonstrated different causal variants.',
      'H4/(H3+H4) is conditional: never interpret it without H0-H4 and regional QTL signal strength.',
      'Variant posterior SNP.PP.H4 is conditional on H4; do not nominate causal genes from it when H4 is low.',
      'sdY=1 is an explicit sensitivity assumption, not a confirmed phenotype scale.',
      'sdY estimation uses study-level N and MAF; effective per-gene/per-variant N and residual phenotype variance remain uncertain.',
      'SNVs only; all available valid SNVs are retained without P-value or LD pruning. Indels require separate normalization.',
      '100 common SNVs plus five candidates is an operational minimum, not a guarantee of dense or sufficient coverage.',
      'QTL cis-window limits may truncate regional coverage. Review pair_coverage_and_signal_QC and matching audit files.',
      'No allele swaps based on allele frequency. Only exact post-liftOver genomic REF/ALT matches enter the analysis.',
      'No QTL LD-based multiple-signal fitting or sample-overlap correction is performed.',
      'Different datasets or scale/region checks are not automatically independent replication.'),file.path(out,'interpretation_notes.txt'))
    if(any(vapply(status_rows,function(x)x$status=='ERROR',logical(1))))
      stop('Some coloc runs failed; partial outputs saved. Inspect run_status.tsv.')
    completed <<- TRUE
    cat('\nCOMPLETED. Output:',out,'\nUpload the complete output folder as a ZIP.\n')
  }
  tryCatch(withCallingHandlers(main(),warning=function(w) {
    warnings_seen <<- c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')
  }),error=function(e) {
    warnings_seen <<- c(warnings_seen,paste('FATAL:',conditionMessage(e)))
    stop(conditionMessage(e),call.=FALSE)
  })
})()
