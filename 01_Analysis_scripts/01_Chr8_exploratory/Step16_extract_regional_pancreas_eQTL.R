# 第16步：三个基因表达QTL数据集的区域原始记录提取和覆盖审计。
# 不按P值筛选，不运行共定位，不翻转等位基因，不删除原始重复记录。
(function() {
  pkgs <- c('data.table','Rsamtools','GenomicRanges','IRanges')
  ok <- vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)
  if(!all(ok)) stop('Missing packages: ',paste(pkgs[!ok],collapse=', '))
  base <- path.expand('~/胰腺')
  meta_file <- file.path(base,'15_QTL_availability','pancreas_QTL_ff7a86fd2dece',
    'eQTL_Catalogue_pancreas_islet_datasets.tsv')
  if(!file.exists(meta_file)) stop('Step-15 metadata not found: ',meta_file)
  meta <- as.data.frame(data.table::fread(meta_file))
  dataset_ids <- c('QTD000296','QTD000554','QTD000574')
  if(anyDuplicated(meta$dataset_id)) stop('Duplicated dataset IDs in metadata')
  meta <- meta[match(dataset_ids,meta$dataset_id),,drop=FALSE]
  if(anyNA(meta$dataset_id)||any(meta$quant_method!='ge')) stop('Missing or incorrect ge datasets')
  # 这是与第15步相同的GRCh38探索性提取窗口，不声称是旧GRCh37窗口的精确liftOver。
  region_start <- 128035132L;region_end <- 129057753L
  candidates <- data.frame(SNP=c('rs4733612','rs1561927','rs7824269','rs7819154','rs12676584'),
    chromosome='8',position=c(128557753L,128555832L,128540098L,128546105L,128535132L),
    ref=c('G','C','C','C','G'),alt=c('A','T','T','A','A'),stringsAsFactors=FALSE)
  candidates$key <- with(candidates,paste(chromosome,position,ref,alt,sep=':'))
  parent <- file.path(base,'16_regional_eQTL')
  dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile('chr8_eQTL_',tmpdir=parent)
  if(!dir.create(out)) stop('Cannot create output folder')
  dir.create(file.path(out,'indices'))
  sink(file.path(out,'run_log.txt'),split=TRUE);on.exit(sink(),add=TRUE)
  oldopt <- options(timeout=180);on.exit(options(oldopt),add=TRUE)
  wt <- function(x,n) data.table::fwrite(x,file.path(out,n),sep='\t',na='NA')
  wt(meta,'selected_dataset_metadata.tsv');wt(candidates,'candidate_manifest_GRCh38.tsv')
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  saveRDS(list(build='GRCh38',chromosome='8',start=region_start,end=region_end,
    purpose='Exploratory regional extraction, no p-value filtering',metadata_input=meta_file,
    candidates=candidates,download_mode='Remote tabix with local index'),file.path(out,'analysis_config.rds'))
  logs <- character(); statuses <- list(); coverage <- list(); all_genes <- list(); targets <- list()
  cat('Output:',out,'\nGRCh38 chr8:',region_start,'-',region_end,'\n',sep='')
  # 读取远程gzip文件的第一行，关闭连接；不下载整个大文件。
  get_header <- function(u) {
    con <- gzcon(url(u,open='rb',method='libcurl'))
    on.exit(close(con),add=TRUE)
    first <- readLines(con,n=1L,warn=FALSE)
    if(length(first)!=1L || !grepl('molecular_trait_id',first,fixed=TRUE)) stop('Unexpected remote header')
    h <- strsplit(sub('^#','',sub('\r$','',first)),'\t',fixed=TRUE)[[1]]
    if(anyDuplicated(h)) stop('Duplicated column names')
    h
  }
  # 每个数据集只做一次区域查询，避免逐SNP连续访问远程文件。
  scan_region <- function(u,idx) {
    tb <- Rsamtools::TabixFile(file=u,index=idx)
    open(tb);on.exit(close(tb),add=TRUE)
    seqs <- Rsamtools::seqnamesTabix(tb)
    chromosome <- intersect(c('8','chr8'),seqs)
    if(length(chromosome)!=1L) stop('Cannot identify unique chr8 sequence in index')
    param <- GenomicRanges::GRanges(chromosome,
      IRanges::IRanges(start=region_start,end=region_end))
    Rsamtools::scanTabix(tb,param=param)[[1L]]
  }
  for(i in seq_len(nrow(meta))) {
    m <- meta[i,];id <- m$dataset_id
    u <- sub('^ftp://','https://',m$ftp_path)
    idx <- file.path(out,'indices',paste0(id,'.all.tsv.gz.tbi'))
    cat('\n===== ',id,' / ',m$study_label,' =====\n',sep='');flush.console()
    status <- data.frame(dataset_id=id,study=m$study_label,metadata_N=m$sample_size,
      status='ERROR',raw_rows=NA_integer_,unique_variants=NA_integer_,
      molecular_traits=NA_integer_,gene_ids=NA_integer_,valid_effect_rows=NA_integer_,
      duplicate_pair_rows=NA_integer_,message='',stringsAsFactors=FALSE)
    tryCatch(withCallingHandlers({
      Sys.sleep(3)
      download.file(paste0(u,'.tbi'),idx,method='libcurl',mode='wb',quiet=TRUE)
      if(!file.exists(idx)||file.info(idx)$size==0) stop('Empty index download')
      Sys.sleep(3)
      header <- get_header(u)
      writeLines(paste(header,collapse='\t'),file.path(out,paste0(id,'_source_header.tsv')))
      Sys.sleep(3)
      records <- scan_region(u,idx)
      # 首先保存原始提取记录，后续任何QC失败也不丢失。
      con <- gzfile(file.path(out,paste0(id,'_region_raw.tsv.gz')),'wt')
      tryCatch(writeLines(c(paste(header,collapse='\t'),records),con),finally=close(con))
      status$raw_rows <- length(records)
      if(!length(records)) {
        status$status <- 'QUERY_RETURNED_ZERO_ROWS_REVIEW_REQUIRED'
      } else {
        d <- data.table::fread(text=paste(records,collapse='\n'),sep='\t',header=FALSE,
          colClasses='character',na.strings=c('NA','.'),showProgress=FALSE)
        if(ncol(d)!=length(header)) stop('Header / data column count mismatch')
        data.table::setnames(d,header)
        required <- c('molecular_trait_id','chromosome','position','ref','alt','variant','gene_id','pvalue','beta','se','maf')
        if(!all(required %in% names(d))) stop('Required fields missing from source')
        chr <- sub('^chr','',d$chromosome)
        pos <- suppressWarnings(as.integer(d$position))
        if(anyNA(pos)||any(chr!='8')||any(pos<region_start|pos>region_end)) stop('Returned coordinates inconsistent with query')
        if(anyNA(d$ref)||anyNA(d$alt)||anyNA(d$molecular_trait_id)||anyNA(d$gene_id)) stop('Missing alleles or trait/gene identifiers')
        d[,key_GRCh38 := paste(chr,pos,ref,alt,sep=':')]
        d[,gene_id_unversioned := sub('\\..*$','',gene_id)]
        beta_num <- suppressWarnings(as.numeric(d$beta));se_num <- suppressWarnings(as.numeric(d$se))
        p_num <- suppressWarnings(as.numeric(d$pvalue));maf_num <- suppressWarnings(as.numeric(d$maf))
        d[,valid_effect := is.finite(beta_num)&is.finite(se_num)&se_num>0]
        d[,valid_pvalue := is.finite(p_num)&p_num>=0&p_num<=1]
        d[,pvalue_zero := is.finite(p_num)&p_num==0]
        d[,valid_MAF := is.finite(maf_num)&maf_num>0&maf_num<=0.5]
        d[,duplicate_trait_variant := duplicated(paste(molecular_trait_id,key_GRCh38,sep='|')) |
             duplicated(paste(molecular_trait_id,key_GRCh38,sep='|'),fromLast=TRUE)]
        d[,dataset_id := id]
        # 注释表保留所有记录；上述标记尚未用于过滤。
        data.table::fwrite(d,file.path(out,paste0(id,'_region_annotated.tsv.gz')),sep='\t',compress='gzip',na='NA')
        status$unique_variants <- data.table::uniqueN(d$key_GRCh38)
        status$molecular_traits <- data.table::uniqueN(d$molecular_trait_id)
        status$gene_ids <- data.table::uniqueN(d$gene_id_unversioned)
        status$valid_effect_rows <- sum(d$valid_effect)
        status$duplicate_pair_rows <- sum(d$duplicate_trait_variant)
        gs <- d[,{
          pv <- suppressWarnings(as.numeric(pvalue));validp <- is.finite(pv)&pv>=0&pv<=1
          list(records=.N,unique_variants=data.table::uniqueN(key_GRCh38),
            molecular_traits=data.table::uniqueN(molecular_trait_id),
            valid_effect_rows=sum(valid_effect),invalid_MAF_rows=sum(!valid_MAF),
            duplicate_flagged_rows=sum(duplicate_trait_variant),
            minimum_nominal_p=if(any(validp)) min(pv[validp]) else NA_real_,
            core_variants_covered=sum(candidates$key %in% key_GRCh38))
        },by=gene_id_unversioned]
        gs[,dataset_id := id];wt(gs,paste0(id,'_gene_coverage.tsv'));all_genes[[id]] <- gs
        hits <- d[d$key_GRCh38 %in% candidates$key,]
        hits[,candidate_rsid := candidates$SNP[match(key_GRCh38,candidates$key)]]
        wt(hits,paste0(id,'_five_candidate_associations.tsv'))
        cv <- lapply(seq_len(nrow(candidates)),function(j) {
          x <- candidates[j,];sel <- d$key_GRCh38==x$key
          data.frame(dataset_id=id,SNP=x$SNP,key_GRCh38=x$key,exact_allele_records=sum(sel),
            genes_with_exact_allele=data.table::uniqueN(d$gene_id_unversioned[sel]),
            same_position_other_allele_records=sum(pos==x$position & !sel),
            status=if(any(sel)) 'EXACT_REF_ALT_PRESENT' else 'NO_EXACT_RECORD_IN_REGION')
        })
        coverage[[id]] <- data.table::rbindlist(cv)
        focus <- c(LINC00824='ENSG00000254275',PVT1='ENSG00000249859',MYC='ENSG00000136997')
        targets[[id]] <- data.frame(dataset_id=id,gene_symbol=names(focus),gene_id=unname(focus),
          regional_record_count=vapply(focus,function(g) sum(d$gene_id_unversioned==g),integer(1)),
          status=vapply(focus,function(g) if(any(d$gene_id_unversioned==g)) 'REGIONAL_ASSOCIATIONS_AVAILABLE' else 'NO_REGIONAL_RECORD_NOT_PROOF_UNTESTED',character(1)))
        status$status <- 'EXTRACTED_AND_AUDITED'
      }
    },warning=function(w) {
      logs <<- c(logs,paste(id,'WARNING',conditionMessage(w)));invokeRestart('muffleWarning')
    }),error=function(e) {
      status$message <<- conditionMessage(e)
      logs <<- c(logs,paste(id,'ERROR',conditionMessage(e)))
      cat('Dataset failed:',conditionMessage(e),'\n')
    })
    statuses[[id]] <- status
    wt(data.table::rbindlist(statuses),'extraction_status.tsv')
    if(length(coverage)) wt(data.table::rbindlist(coverage),'five_candidate_coverage.tsv')
    if(length(all_genes)) wt(data.table::rbindlist(all_genes,fill=TRUE),'all_gene_coverage.tsv')
    if(length(targets)) wt(data.table::rbindlist(targets),'priority_gene_availability.tsv')
    writeLines(if(length(logs)) logs else 'No captured warnings/errors.',file.path(out,'warnings_and_errors.txt'))
    print(status,row.names=FALSE)
  }
  writeLines(c('Coordinates: GRCh38. ALT is the QTL effect allele according to eQTL Catalogue documentation.',
    'No P-value filter was applied. Minimum P-values are descriptive, not corrected significance tests.',
    'All original records are retained, including rsID-related duplicates and invalid-effect flags.',
    'NO_REGIONAL_RECORD does not establish that a gene was unmeasured or unexpressed.',
    'Failed queries and zero-row queries are not negative biological evidence.',
    'The extraction window is exploratory; GWAS coordinate conversion and pair-level deduplication remain necessary.',
    'Exact REF/ALT matches are reported; no strand or effect flipping has been performed.',
    'This step extracts gene-expression QTLs only, not splicing QTLs.',
    'GTEx Catalogue reanalysis is not independent replication of GTEx Portal data.',
    'No colocalisation was run. Do not interpret candidate coverage as colocalisation evidence.'),file.path(out,'interpretation_notes.txt'))
  cat('\n===== Final extraction status =====\n');print(as.data.frame(data.table::rbindlist(statuses)),row.names=FALSE)
  cat('\nOutput:',out,'\nIf any dataset failed, retain this folder and send the error log; do not repeatedly retry remote queries.\n')
})()
