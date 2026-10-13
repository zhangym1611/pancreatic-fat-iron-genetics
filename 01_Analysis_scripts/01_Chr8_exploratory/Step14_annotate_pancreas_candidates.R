# 第14步：候选变异注释与证据检索。需要服务器可访问Ensembl和Europe PMC。
# 自动输出的是注释与检索记录；不是已完成的QTL共定位或系统文献综述。
(function() {
  packages <- c('httr','jsonlite','data.table')
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if (length(missing)) {
    userlib <- path.expand('~/R_pancreas_library')
    dir.create(userlib,recursive=TRUE,showWarnings=FALSE)
    .libPaths(c(userlib,.libPaths()))
    install.packages(missing,lib=userlib,repos='https://cloud.r-project.org')
  }
  if (!all(vapply(packages,requireNamespace,logical(1),quietly=TRUE))) stop('Required packages unavailable')
  base <- path.expand('~/胰腺')
  parent <- file.path(base,'14_functional_annotation')
  dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile('chr8_annotation_',tmpdir=parent)
  if (!dir.create(out)) stop('Cannot create output directory')
  dir.create(file.path(out,'raw_json'))
  sink(file.path(out,'run_log.txt'),split=TRUE); on.exit(sink(),add=TRUE)
  cat('Output:',out,'\n')
  candidates <- data.frame(SNP=c('rs4733612','rs1561927','rs7824269','rs7819154','rs12676584'),
    chr='8',pos_GRCh37=c(129569999L,129568078L,129552344L,129558351L,129547378L),
    REF=c('G','C','C','C','G'),ALT=c('A','T','T','A','A'),stringsAsFactors=FALSE)
  write_tab <- function(x,n) data.table::fwrite(x,file.path(out,n),sep='\t',na='NA')
  write_tab(candidates,'candidate_manifest.tsv')
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  requests <- list(); errors <- character()
  scalar <- function(x,default=NA_character_) if(is.null(x)||!length(x)) default else paste(unlist(x),collapse=';')
  fetch <- function(url,query=list(),label) {
    cat('Request:',label,'\n'); flush.console()
    status <- NA_integer_; message <- ''; payload <- NULL
    tryCatch({
      for (attempt in 1:3) {
        resp <- httr::GET(url,query=query,httr::accept_json(),httr::timeout(45),
                          httr::user_agent('PancreasCandidateAnnotation/1.0'))
        status <- httr::status_code(resp)
        body <- httr::content(resp,as='text',encoding='UTF-8')
        writeLines(body,file.path(out,'raw_json',paste0(label,'.json')),useBytes=TRUE)
        if (status==200L) break
        if (!(status %in% c(429L,500L,502L,503L,504L)) || attempt==3L) break
        Sys.sleep(3)
      }
      if(status!=200L) stop('HTTP ',status)
      payload <- jsonlite::fromJSON(body,simplifyVector=FALSE)
      if (!is.null(payload$error)) stop(scalar(payload$error))
    },error=function(e) {message <<- conditionMessage(e); payload <<- NULL})
    if(nzchar(message)) errors <<- c(errors,paste(label,message))
    requests[[label]] <<- data.frame(label=label,url=url,
      query=jsonlite::toJSON(query,auto_unbox=TRUE),HTTP_status=status,
      status=if(nzchar(message)) 'FAILED' else 'retrieved',message=message,
      retrieved_UTC=format(Sys.time(),tz='UTC',usetz=TRUE))
    write_tab(data.table::rbindlist(requests),'request_log.tsv')
    Sys.sleep(0.3)
    payload
  }
  # GRCh37端点与现有GWAS坐标一致，不进行隐式坐标转换。
  host <- 'https://grch37.rest.ensembl.org'
  release <- fetch(paste0(host,'/info/data'),list('content-type'='application/json'),'ensembl_GRCh37_release')
  summaries <- list(); consequences <- list(); mappings <- list()
  for(i in seq_len(nrow(candidates))) {
    d <- candidates[i,]; rs <- d$SNP
    v <- fetch(paste0(host,'/vep/human/id/',rs),
      list('content-type'='application/json',canonical=1,symbol=1,hgvs=1,regulatory=1),paste0(rs,'_VEP_GRCh37'))
    accepted <- list()
    if(!is.null(v)) for (entry in v) {
      alleles <- strsplit(scalar(entry$allele_string,''),'/',fixed=TRUE)[[1]]
      matched <- identical(scalar(entry$assembly_name),'GRCh37') &&
        identical(scalar(entry$seq_region_name),'8') &&
        isTRUE(as.numeric(scalar(entry$start))==d$pos_GRCh37) &&
        isTRUE(as.numeric(scalar(entry$end))==d$pos_GRCh37) &&
        isTRUE(as.numeric(scalar(entry$strand))==1) &&
        length(alleles)>=2L && identical(alleles[1],d$REF) && d$ALT %in% alleles[-1]
      if(matched) accepted[[length(accepted)+1L]] <- entry
    }
    state <- if(is.null(v)) 'REQUEST_FAILED' else if(!length(accepted)) 'COORDINATE_OR_ALLELE_REVIEW' else 'MATCHED_GRCh37'
    genes <- character(); terms <- character()
    for (entry in accepted) {
      for (kind in c('transcript_consequences','regulatory_feature_consequences','motif_feature_consequences','intergenic_consequences')) {
        for (a in entry[[kind]]) {
          if(!identical(scalar(a$variant_allele),d$ALT)) next
          genes <- c(genes,scalar(a$gene_symbol,'')); terms <- c(terms,unlist(a$consequence_terms))
          consequences[[length(consequences)+1L]] <- data.frame(SNP=rs,assembly='GRCh37',
            effect_allele=d$ALT,annotation_type=kind,gene_id=scalar(a$gene_id),
            gene_symbol=scalar(a$gene_symbol),transcript_id=scalar(a$transcript_id),
            consequence=scalar(a$consequence_terms),impact=scalar(a$impact),
            biotype=scalar(a$biotype),canonical=scalar(a$canonical),
            hgvsc=scalar(a$hgvsc),hgvsp=scalar(a$hgvsp),
            regulatory_feature_id=scalar(a$regulatory_feature_id),
            motif_feature_id=scalar(a$motif_feature_id))
        }
      }
    }
    summaries[[rs]] <- data.frame(SNP=rs,pos_GRCh37=d$pos_GRCh37,REF=d$REF,ALT=d$ALT,
      VEP_status=state,annotated_gene_symbols=paste(sort(unique(genes[nzchar(genes)])),collapse=';'),
      consequences=paste(sort(unique(terms)),collapse=';'))
    # GRCh38仅记录rsID映射，供后续GTEx等资源核查；不自动拿来匹配QTL。
    m <- fetch(paste0('https://rest.ensembl.org/variation/human/',rs),
       list('content-type'='application/json'),paste0(rs,'_mapping_GRCh38'))
    if(!is.null(m)) for(a in m$mappings) mappings[[length(mappings)+1L]] <- data.frame(
      SNP=rs,assembly=scalar(a$assembly_name),chr=scalar(a$seq_region_name),
      start=scalar(a$start),end=scalar(a$end),strand=scalar(a$strand),alleles=scalar(a$allele_string))
  }
  write_tab(data.table::rbindlist(summaries),'variant_annotation_summary.tsv')
  if(length(consequences)) write_tab(data.table::rbindlist(consequences,fill=TRUE),'all_ALT_consequences.tsv')
  if(length(mappings)) write_tab(data.table::rbindlist(mappings),'rsID_coordinate_mappings.tsv')
  # 整个既定共定位区域内的基因；位置接近不等于靶基因。
  g <- fetch(paste0(host,'/overlap/region/human/8:129047378-130047378'),
    list(feature='gene','content-type'='application/json'),'regional_genes_GRCh37')
  gene_rows <- list()
  if(!is.null(g)) for(a in g) {
    if(!identical(scalar(a$assembly_name),'GRCh37') || !identical(scalar(a$seq_region_name),'8')) next
    for(i in seq_len(nrow(candidates))) {
      start <- as.numeric(scalar(a$start)); end <- as.numeric(scalar(a$end)); strand <- as.numeric(scalar(a$strand))
      tss <- if(isTRUE(strand==1)) start else end
      pos <- candidates$pos_GRCh37[i]
      gene_rows[[length(gene_rows)+1L]] <- data.frame(SNP=candidates$SNP[i],gene_id=scalar(a$id),
        gene_symbol=scalar(a$external_name),biotype=scalar(a$biotype),start=start,end=end,
        strand=strand,TSS=tss,distance_to_TSS_bp=abs(pos-tss),within_gene=pos>=start & pos<=end)
    }
  }
  if(length(gene_rows)) {
    tab <- data.table::rbindlist(gene_rows);tab <- tab[order(tab$SNP,tab$distance_to_TSS_bp),]
    write_tab(tab,'regional_gene_distances.tsv')
  }
  # 文献检索为首轮发现，记录总命中数；最多保存每个查询前100条，不宣称全面覆盖。
  literature <- list(); search_log <- list()
  queries <- setNames(paste0('"',candidates$SNP,'"'),candidates$SNP)
  queries[['pancreas_fat_iron_volume']] <- '(pancreas OR pancreatic) AND (fat OR iron OR volume) AND (colocalisation OR colocalization OR "fine mapping")'
  for(label in names(queries)) {
    q <- queries[[label]]
    z <- fetch('https://www.ebi.ac.uk/europepmc/webservices/rest/search',
      list(query=q,format='json',pageSize=100,resultType='core'),paste0('literature_',label))
    search_log[[label]] <- data.frame(query_label=label,query=q,
      hit_count=if(is.null(z)) NA_character_ else scalar(z$hitCount),
      status=if(is.null(z)) 'REQUEST_FAILED' else 'retrieved_first_100_at_most')
    if(!is.null(z)) for(a in z$resultList$result) literature[[length(literature)+1L]] <- data.frame(
      query_label=label,id=scalar(a$id),source=scalar(a$source),PMID=scalar(a$pmid),
      PMCID=scalar(a$pmcid),DOI=scalar(a$doi),title=scalar(a$title),year=scalar(a$pubYear),
      journal=scalar(a$journalInfo$journal$title))
  }
  write_tab(data.table::rbindlist(search_log),'literature_search_log.tsv')
  if(length(literature)) write_tab(data.table::rbindlist(literature),'literature_hits.tsv')
  links <- data.frame(SNP=candidates$SNP,
    Ensembl_GRCh37=paste0('https://grch37.ensembl.org/Homo_sapiens/Variation/Explore?v=',candidates$SNP),
    GWAS_Catalog=paste0('https://www.ebi.ac.uk/gwas/variants/',candidates$SNP),
    dbSNP=paste0('https://www.ncbi.nlm.nih.gov/snp/',candidates$SNP),
    PubMed=paste0('https://pubmed.ncbi.nlm.nih.gov/?term=',candidates$SNP),
    GTEx='https://gtexportal.org/home/',eQTL_Catalogue='https://www.ebi.ac.uk/eqtl/',
    QTL_evidence_status='NOT_YET_REVIEWED',target_tissue='Pancreas; islet datasets where available')
  write_tab(links,'evidence_review_links.tsv')
  writeLines(if(length(errors)) errors else 'No captured request errors.',file.path(out,'warnings_and_errors.txt'))
  writeLines(c('This is annotation and initial literature discovery, not causal-gene identification.',
    'No eQTL/sQTL association or QTL colocalisation has been performed in this script.',
    'Zero search hits, empty annotation and failed requests do not establish absence of evidence.',
    'Nearby genes and VEP transcript assignments are not confirmed causal targets.',
    'GRCh38 mappings require allele/build verification before QTL matching.',
    'Original GWAS supplementary tables and newer studies require full-text review.',
    'Candidate list is the union of retained credible sets, not an exhaustive list of possible causal variants.'),
    file.path(out,'interpretation_notes.txt'))
  cat('\n===== Annotation summary =====\n');print(as.data.frame(data.table::rbindlist(summaries)),row.names=FALSE)
  cat('\nOutput:',out,'\nCaptured request failures:',length(errors),'\n')
  cat('QTL evidence remains NOT_YET_REVIEWED. Keep all raw responses and request logs.\n')
})()
