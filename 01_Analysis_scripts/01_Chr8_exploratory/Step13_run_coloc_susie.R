# 第13步：复用第12步拟合对象进行SuSiE共定位；不重新拟合、不修改输入。
(function() {
  required <- c('coloc', 'susieR', 'data.table')
  available <- vapply(required, requireNamespace, logical(1), quietly=TRUE)
  if (!all(available)) stop('Missing packages: ', paste(required[!available], collapse=', '))
  base <- path.expand('~/胰腺')
  scenarios <- c('baseline_2729','exclude_rs117601696','reference_MAF_ge_005')
  traits <- c('pancreas_fat','pancreas_iron','pancreas_volume')
  fit_names <- unlist(lapply(scenarios, function(s) paste0(s,'__',traits,'_susie_fit.rds')))
  # 唯一完整目录时自动读取；多个完整目录时由用户选择配置文件，避免混用批次。
  configs <- list.files(file.path(base,'12_SuSiE_sensitivity'),
    pattern='^analysis_config\\.rds$', recursive=TRUE, full.names=TRUE)
  configs <- configs[vapply(dirname(configs), function(d)
    all(file.exists(file.path(d,fit_names))), logical(1))]
  if (!length(configs)) stop('No complete step-12 folder found under ~/胰腺/12_SuSiE_sensitivity')
  if (length(configs)==1L) {
    config_file <- configs[1]
  } else {
    cat('Multiple complete runs found. Select analysis_config.rds from the intended run.\n')
    print(configs)
    if (!interactive()) stop('Run interactively to select the intended input batch.')
    config_file <- file.choose()
    if (!normalizePath(config_file,mustWork=TRUE) %in% normalizePath(configs,mustWork=TRUE))
      stop('Selected file is not a configuration from a complete step-12 run.')
  }
  input_dir <- dirname(config_file)
  cfg <- readRDS(config_file)
  if (!all(scenarios %in% names(cfg$scenarios))) stop('Missing scenario identifiers')
  fits <- list(); audit <- list(); members <- list(); components <- list()
  for (s in scenarios) {
    ids <- as.character(cfg$scenarios[[s]])
    if (!length(ids) || anyNA(ids) || anyDuplicated(ids)) stop('Invalid configured SNPs: ',s)
    for (t in traits) {
      label <- paste(s,t,sep='__')
      f <- readRDS(file.path(input_dir,paste0(label,'_susie_fit.rds')))
      if (!inherits(f,'susie') || !isTRUE(f$converged)) stop('Invalid or nonconverged fit: ',label)
      for (field in c('alpha','lbf_variable')) {
        m <- f[[field]]
        if (!is.matrix(m) || !identical(colnames(m),ids) || any(!is.finite(m)))
          stop('Invalid matrix or SNP order: ',label,' / ',field)
      }
      if (!identical(dim(f$alpha),dim(f$lbf_variable))) stop('Matrix dimension mismatch: ',label)
      if (is.null(f$sets$cs) || !length(f$sets$cs)) stop('No retained credible set: ',label)
      if (length(f$sets$cs_index)!=length(f$sets$cs)) stop('CS indices missing: ',label)
      # 原对象已经带有SNP列名；再次核验后为PIP附上同一顺序的名称。
      if (length(f$pip)!=length(ids)) stop('PIP length mismatch: ',label)
      names(f$pip) <- ids
      fits[[label]] <- f
      audit[[label]] <- data.frame(scenario=s,trait=t,nsnps=length(ids),
        converged=f$converged,n_CS=length(f$sets$cs),snp_order='PASS')
      for (j in seq_along(f$sets$cs)) {
        ix <- f$sets$cs[[j]]; component <- f$sets$cs_index[j]
        if (anyNA(ix) || any(ix<1 | ix>length(ids))) stop('Invalid CS membership: ',label)
        members[[paste(label,j)]] <- data.frame(scenario=s,trait=t,
          CS=names(f$sets$cs)[j],component=component,SNP=ids[ix],
          alpha=f$alpha[component,ix],PIP=f$pip[ix])
      }
      components[[label]] <- data.frame(scenario=s,trait=t,
        component=seq_len(nrow(f$alpha)),prior_variance=f$V,logBF=f$lbf,
        retained_CS=seq_len(nrow(f$alpha)) %in% f$sets$cs_index)
    }
  }
  parent <- file.path(base,'13_coloc_susie')
  dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile('chr8_coloc_susie_',tmpdir=parent)
  if (!dir.create(out)) stop('Cannot create output directory')
  sink(file.path(out,'run_log.txt'),split=TRUE)
  on.exit(sink(),add=TRUE)
  cat('Input:',input_dir,'\nOutput:',out,'\n')
  cat('Reusing existing fits; no SuSiE refitting.\n')
  write_tab <- function(d,name) data.table::fwrite(d,file.path(out,name),sep='\t')
  write_tab(data.table::rbindlist(audit),'input_audit.tsv')
  write_tab(data.table::rbindlist(members),'credible_set_members.tsv')
  write_tab(data.table::rbindlist(components),'all_component_audit.tsv')
  input_files <- c(config_file,file.path(input_dir,fit_names))
  write_tab(data.frame(file=basename(input_files),md5=unname(tools::md5sum(input_files))),
            'input_checksums.tsv')
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  priors <- c(1e-6,1e-5,5e-5)
  saveRDS(list(input_dir=input_dir,source_config=cfg,p1=1e-4,p2=1e-4,p12=priors,
    back_calculate_lbf=FALSE,analysis='Pairwise coloc.susie of retained credible sets'),
    file.path(out,'step13_config.rds'))
  pairs <- list(fat_iron=traits[c(1,2)],fat_volume=traits[c(1,3)],iron_volume=traits[c(2,3)])
  results <- list(); statuses <- list(); notes <- character()
  pp_cols <- paste0('PP.H',0:4,'.abf')
  for (s in scenarios) for (pair in names(pairs)) for (p12 in priors) {
    ts <- pairs[[pair]]; label <- paste(s,pair,format(p12,scientific=TRUE),sep='__')
    cat('\nRunning:',label,'\n'); flush.console()
    ans <- tryCatch(withCallingHandlers({
      x <- coloc::coloc.susie(
        dataset1=fits[[paste(s,ts[1],sep='__')]],
        dataset2=fits[[paste(s,ts[2],sep='__')]],
        back_calculate_lbf=FALSE,p1=1e-4,p2=1e-4,p12=p12)
      saveRDS(x,file.path(out,paste0(label,'_coloc_susie.rds')))
      sm <- as.data.frame(x$summary)
      if (!nrow(sm) || !all(pp_cols %in% names(sm))) stop('No valid signal-pair posterior table returned')
      pp <- as.matrix(sm[,pp_cols,drop=FALSE])
      if (any(!is.finite(pp)) || any(pp<0 | pp>1) || any(abs(rowSums(pp)-1)>1e-6))
        stop('Invalid posterior probabilities')
      sm <- cbind(data.frame(scenario=s,pair=pair,trait1=ts[1],trait2=ts[2],
        p1=1e-4,p2=1e-4,p12=p12,stringsAsFactors=FALSE),sm)
      if (!is.null(x$results)) write_tab(as.data.frame(x$results),paste0(label,'_variant_posteriors.tsv'))
      sm
    },warning=function(w) {
      notes <<- c(notes,paste(label,'WARNING',conditionMessage(w)))
      invokeRestart('muffleWarning')
    }),error=function(e) {
      notes <<- c(notes,paste(label,'ERROR',conditionMessage(e)))
      NULL
    })
    statuses[[label]] <- data.frame(scenario=s,pair=pair,p12=p12,
      status=if(is.null(ans)) 'ERROR' else 'completed',
      signal_pair_rows=if(is.null(ans)) 0L else nrow(ans))
    if (!is.null(ans)) results[[label]] <- ans
    write_tab(data.table::rbindlist(statuses),'run_status.tsv')
    if (length(results)) write_tab(data.table::rbindlist(results,fill=TRUE),'coloc_susie_summary.tsv')
    writeLines(if(length(notes)) notes else 'No captured warnings/errors.',
               file.path(out,'warnings_and_errors.txt'))
  }
  if (length(results)) {
    all_results <- as.data.frame(data.table::rbindlist(results,fill=TRUE))
    main <- all_results[all_results$p12==1e-5,,drop=FALSE]
    write_tab(main,'coloc_susie_main_prior.tsv')
    cat('\n===== Main prior p12=1e-5 =====\n'); print(main,row.names=FALSE)
  }
  cat('\n===== Run status =====\n'); print(as.data.frame(data.table::rbindlist(statuses)),row.names=FALSE)
  cat('\nOutput:',out,'\n')
  cat('Interpret only retained signal comparisons. Unresolved secondary components remain in all_component_audit.tsv.\n')
  cat('Pairwise colocalisation does not establish three-trait sharing, mediation or causal direction.\n')
  if (any(vapply(statuses,function(x) x$status=='ERROR',logical(1))))
    cat('Some comparisons failed; inspect warnings_and_errors.txt before interpretation.\n')
})()
