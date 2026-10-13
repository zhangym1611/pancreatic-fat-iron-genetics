# 第15步：胰腺QTL可用性核查。固定GTEx v8 / GENCODE v26 / GRCh38。
# 显著关联接口的空结果不等于未检测、不表达或无调控作用。
(function() {
  userlib <- path.expand('~/R_pancreas_library')
  if(dir.exists(userlib)) .libPaths(c(userlib,.libPaths()))
  pkgs <- c('httr','jsonlite','data.table')
  missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) {
    dir.create(userlib,recursive=TRUE,showWarnings=FALSE)
    .libPaths(c(userlib,.libPaths()))
    install.packages(missing,lib=userlib,repos='https://cloud.r-project.org')
  }
  if(!all(vapply(pkgs,requireNamespace,logical(1),quietly=TRUE))) stop('Packages unavailable')
  base <- path.expand('~/胰腺')
  input <- file.path(base,'14_functional_annotation','chr8_annotation_ff7a87bd402b5')
  if(!file.exists(file.path(input,'regional_gene_distances.tsv'))) stop('Step-14 input not found: ',input)
  parent <- file.path(base,'15_QTL_availability')
  dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile('pancreas_QTL_',tmpdir=parent);dir.create(out);dir.create(file.path(out,'raw'))
  sink(file.path(out,'run_log.txt'),split=TRUE);on.exit(sink(),add=TRUE)
  cat('Output:',out,'\nGTEx v8; GENCODE v26; GRCh38; tissue Pancreas\n')
  logs <- list();issues <- character()
  wt <- function(d,n) data.table::fwrite(d,file.path(out,n),sep='\t',na='NA')
  scalar <- function(x) if(is.null(x)||!length(x)) NA_character_ else paste(unlist(x),collapse=';')
  request <- function(url,q=list(),label,json=TRUE) {
    status <- NA_integer_;msg <- '';value <- NULL
    tryCatch({
      for(attempt in 1:3) {
        r <- httr::GET(url,query=q,httr::timeout(45),httr::user_agent('PancreasQTL-audit/1.0'))
        status <- httr::status_code(r);body <- httr::content(r,as='text',encoding='UTF-8')
        writeLines(body,file.path(out,'raw',paste0(label,if(json) '.json' else '.tsv')),useBytes=TRUE)
        if(status==200 || !(status %in% c(429,500,502,503,504)) || attempt==3) break
        Sys.sleep(3)
      }
      if(status!=200) stop('HTTP ',status)
      value <- if(json) jsonlite::fromJSON(body,simplifyVector=FALSE) else body
    },error=function(e) {msg <<- conditionMessage(e);value <<- NULL})
    if(nzchar(msg)) issues <<- c(issues,paste(label,msg))
    logs[[label]] <<- data.frame(label=label,url=url,query=jsonlite::toJSON(q,auto_unbox=TRUE),
      HTTP_status=status,status=if(nzchar(msg)) 'FAILED' else 'retrieved',message=msg,
      UTC=format(Sys.time(),tz='UTC',usetz=TRUE))
    wt(data.table::rbindlist(logs),'request_log.tsv');Sys.sleep(0.3)
    value
  }
  # 将嵌套字段保留为JSON字符串，不丢弃API返回的字段。
  flatten <- function(rows) {
    if(!length(rows)) return(data.table::data.table())
    data.table::rbindlist(lapply(rows,function(r) as.data.frame(lapply(r,function(v) {
      if(is.null(v)||!length(v)) NA_character_ else if(is.atomic(v)&&length(v)==1) as.character(v)
      else as.character(jsonlite::toJSON(v,auto_unbox=TRUE))
    }),stringsAsFactors=FALSE)),fill=TRUE)
  }
  api <- 'https://gtexportal.org/api/v2/'
  paged <- function(endpoint,q,label) {
    rows <- list();complete <- FALSE
    for(page in 0:49) {
      q$page <- page;q$itemsPerPage <- 1000L
      z <- request(paste0(api,endpoint),q,paste0(label,'_p',page))
      if(is.null(z)) break
      if(!'data' %in% names(z)) {issues <<- c(issues,paste(label,'unexpected API schema'));break}
      rows <- c(rows,z$data)
      np <- z$paging_info$numberOfPages
      if(!is.null(np) && page+1>=np) {complete <- TRUE;break}
      if(is.null(np) && length(z$data)<1000) {complete <- TRUE;break}
    }
    list(rows=rows,complete=complete)
  }
  # 区域显著eQTL扫描：用于发现其他靶基因；不把空结果当作全区域无QTL。
  region <- request(paste0(api,'association/singleTissueEqtlByLocation'),
    list(tissueSiteDetailId='Pancreas',start=128035132L,end=129057753L,
         chromosome='chr8',datasetId='gtex_v8'),'GTEx_region_significant_eQTL')
  regional_rows <- if(!is.null(region)) region$singleTissueEqtl else NULL
  if(length(regional_rows)) wt(flatten(regional_rows),'GTEx_region_significant_eQTL.tsv')
  region_schema <- !is.null(region) && 'singleTissueEqtl' %in% names(region)
  wt(data.frame(start_GRCh38=128035132,end_GRCh38=129057753,
    n_returned=if(region_schema) length(regional_rows) else NA_integer_,
    status=if(region_schema) 'SIGNIFICANT_ASSOCIATION_QUERY_COMPLETED' else 'REQUEST_OR_SCHEMA_FAILED'),
    'regional_query_status.tsv')
  # 基因候选来自第14步；另纳入远端背景基因MYC，不能据此认定它是靶基因。
  old <- data.table::fread(file.path(input,'regional_gene_distances.tsv'))
  genes <- unique(as.character(old$gene_id))
  genes <- unique(c('ENSG00000254275',genes,'ENSG00000136997'))
  if(length(regional_rows)) {
    extra <- unlist(lapply(regional_rows,function(x) x$gencodeId))
    genes <- unique(c(genes,sub('\\..*$','',extra)))
  }
  genes <- genes[!is.na(genes)&grepl('^ENSG[0-9]+$',genes)]
  wt(data.frame(gene_id=genes),'candidate_gene_ids.tsv')
  summary <- list()
  for(gene in genes) {
    cat('\nChecking:',gene,'\n');flush.console()
    row <- data.frame(gene_id=gene,gene_symbol=NA_character_,gencodeId=NA_character_,
      reference_status='REQUEST_FAILED',median_TPM=NA_real_,expression_status='NOT_QUERIED',
      significant_eQTL_records=NA_integer_,eQTL_status='NOT_QUERIED',
      significant_sQTL_records=NA_integer_,sQTL_status='NOT_QUERIED',
      expression_filter_status='NOT_ESTABLISHED_BY_THIS_AUDIT')
    ref <- paged('reference/gene',list(geneId=gene,gencodeVersion='v26',genomeBuild='GRCh38/hg38'),paste0(gene,'_reference'))
    exact <- Filter(function(x) identical(sub('\\..*$','',scalar(x$gencodeId)),gene),ref$rows)
    if(ref$complete && length(exact)==1L) {
      v <- exact[[1]];gid <- scalar(v$gencodeId)
      row$reference_status <- 'FOUND';row$gene_symbol <- scalar(v$geneSymbol);row$gencodeId <- gid
      q <- list(gencodeId=gid,tissueSiteDetailId='Pancreas',datasetId='gtex_v8')
      ex <- paged('expression/medianGeneExpression',q,paste0(gene,'_expression'))
      if(length(ex$rows)) wt(flatten(ex$rows),paste0(gene,'_expression.tsv'))
      row$expression_status <- if(!ex$complete) 'REQUEST_FAILED_OR_INCOMPLETE' else if(!length(ex$rows)) 'NO_EXPRESSION_RECORD_RETURNED' else 'RECORD_RETURNED'
      if(length(ex$rows)==1L) row$median_TPM <- suppressWarnings(as.numeric(scalar(ex$rows[[1]]$median)))
      for(kind in c('eQTL','sQTL')) {
        endpoint <- if(kind=='eQTL') 'association/singleTissueEqtl' else 'association/singleTissueSqtl'
        a <- paged(endpoint,q,paste0(gene,'_',kind))
        if(length(a$rows)) wt(flatten(a$rows),paste0(gene,'_significant_',kind,'.tsv'))
        row[[paste0('significant_',kind,'_records')]] <- if(a$complete) length(a$rows) else NA_integer_
        row[[paste0(kind,'_status')]] <- if(!a$complete) 'REQUEST_FAILED_OR_INCOMPLETE' else if(!length(a$rows)) 'NO_SIGNIFICANT_RECORD_RETURNED' else 'SIGNIFICANT_RECORDS_RETURNED'
      }
    } else if(ref$complete) row$reference_status <- if(!length(exact)) 'NOT_IN_REFERENCE_RESPONSE' else 'AMBIGUOUS_REFERENCE_RESPONSE'
    summary[[gene]] <- row;wt(data.table::rbindlist(summary),'GTEx_pancreas_gene_availability.tsv')
  }
  # eQTL Catalogue旧REST API已停用，使用官方全量下载清单；本步不下载大体积汇总统计。
  metadata_url <- 'https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv'
  body <- request(metadata_url,label='eQTL_Catalogue_download_metadata',json=FALSE)
  if(!is.null(body)) tryCatch({
    m <- data.table::fread(text=body,sep='\t',header=TRUE)
    needed <- c('study_id','dataset_id','tissue_label','sample_group','quant_method','ftp_path')
    if(!all(needed %in% names(m))) stop('Catalogue metadata schema changed')
    keep <- grepl('pancrea|islet',paste(m$tissue_label,m$sample_group),ignore.case=TRUE)
    selected <- m[which(keep)]
    selected$summary_https <- sub('^ftp://','https://',selected$ftp_path)
    selected$index_https <- paste0(selected$summary_https,'.tbi')
    selected$file_access_status <- 'LISTED_IN_OFFICIAL_METADATA_NOT_DOWNLOADED'
    selected$LINC00824_tested_status <- 'NOT_YET_CHECKED_IN_DATASET'
    wt(selected,'eQTL_Catalogue_pancreas_islet_datasets.tsv')
  },error=function(e) {issues <<- c(issues,paste('Catalogue parsing:',conditionMessage(e)))})
  writeLines(if(length(issues)) issues else 'No captured errors.',file.path(out,'warnings_and_errors.txt'))
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  saveRDS(list(GTEx='gtex_v8',annotation='GENCODE v26',build='GRCh38',tissue='Pancreas',
    input_step14=input,genes=genes,metadata_url=metadata_url),file.path(out,'analysis_config.rds'))
  writeLines(c('This is an availability audit, not QTL colocalisation.',
    'Median TPM=0 does not prove absence of expression in every sample or cell type.',
    'Expression records do not establish that a gene passed the QTL expression filter.',
    'Empty significant-QTL results do not distinguish untested genes from tested non-significant genes.',
    'A listed Catalogue dataset does not establish that LINC00824 was quantified in it.',
    'Full regional summary statistics, including non-significant variants, are required before colocalisation.',
    'The regional GTEx scan is significant-eQTL discovery only; no complete regional sQTL scan was performed.',
    'GTEx Catalogue reanalysis is not independent replication of GTEx Portal results.',
    'GTEx gene-level association tables include variants outside the imaging interval; allele/build matching is still required.'),
    file.path(out,'interpretation_notes.txt'))
  cat('\n===== GTEx pancreas availability =====\n');print(as.data.frame(data.table::rbindlist(summary)),row.names=FALSE)
  cat('\nOutput:',out,'\nCaptured issues:',length(issues),'\n')
})()
