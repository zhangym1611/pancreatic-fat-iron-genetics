# Exploratory SuSiE-RSS and matched-variant ABF sensitivity analyses.
# Run the entire file in the user's remote RStudio. No automatic deletions.
(function() {
  pkgs <- c('data.table', 'susieR', 'coloc')
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) stop('Missing packages: ', paste(pkgs[!ok], collapse=', '))
  base <- path.expand('~/胰腺')
  input <- file.path(base, '11_GWAS_LD_diagnostics',
    'chr8_diagnostic_c35bc64d238b', 'LD_GWAS_after_frequency_QC.rds')
  if (!file.exists(input)) stop('Input not found: ', input)
  obj <- readRDS(input)
  traits <- c('pancreas_fat', 'pancreas_iron', 'pancreas_volume')
  keys <- c('rs4733612','rs1561927','rs7824269','rs7819154','rs12676584')
  R <- obj$LD
  if (is.null(R) || ncol(R)!=2729L) stop('Expected the audited 2729-variant LD matrix')
  if (!is.matrix(R) || nrow(R)!=ncol(R) || is.null(colnames(R)) ||
      !identical(rownames(R),colnames(R)) || anyDuplicated(colnames(R))) stop('Invalid LD identifiers')
  if (any(!is.finite(R)) || max(abs(R-t(R)))>1e-8 || max(abs(diag(R)-1))>1e-8) stop('Invalid LD values')
  if ('rs200708125' %in% colnames(R)) stop('Frequency-discrepant SNP still present')
  if (!all(traits %in% names(obj$aligned))) stop('Missing trait inputs')
  for (t in traits) {
    d <- obj$aligned[[t]]
    if (!identical(d$snp,colnames(R))) stop('GWAS-LD order mismatch: ',t)
    if (any(!is.finite(d$beta)) || any(!is.finite(d$varbeta)) || any(d$varbeta<=0) ||
        any(!is.finite(d$MAF)) || any(d$MAF<=0 | d$MAF>0.5) ||
        length(d$N)!=1L || !is.finite(d$N) || d$N<=2) stop('Invalid GWAS values: ',t)
  }
  meta <- as.data.frame(obj$retained_variants)
  if (!all(c('rsid','ref_AF_A1') %in% names(meta)) || anyDuplicated(meta$rsid))
    stop('Missing or duplicated retained-variant metadata')
  meta <- meta[match(colnames(R),meta$rsid),,drop=FALSE]
  if (!identical(meta$rsid,colnames(R)) || any(!is.finite(meta$ref_AF_A1)) ||
      any(meta$ref_AF_A1<=0 | meta$ref_AF_A1>=1)) stop('Reference frequency mismatch')
  ref_maf <- pmin(meta$ref_AF_A1,1-meta$ref_AF_A1)
  masks <- list(
    baseline_2729 = rep(TRUE,ncol(R)),
    exclude_rs117601696 = colnames(R)!='rs117601696',
    reference_MAF_ge_005 = ref_maf>=0.05
  )
  if (sum(colnames(R)=='rs117601696')!=1L) stop('Sensitivity SNP not uniquely present')
  for (nm in names(masks)) {
    if (sum(masks[[nm]])<100L || !all(keys %in% colnames(R)[masks[[nm]]])) stop('Unexpected subset: ',nm)
  }
  parent <- file.path(base,'12_SuSiE_sensitivity')
  dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile('chr8_susie_',tmpdir=parent)
  if (!dir.create(out)) stop('Cannot create output folder')
  sink(file.path(out,'run_log.txt'),split=TRUE)
  on.exit(sink(),add=TRUE)
  cat('Output:',out,'\n')
  cat('susieR:',as.character(packageVersion('susieR')),' coloc:',as.character(packageVersion('coloc')),'\n')
  cat('Exploratory external-LD analysis; L=10 is a cap, not the detected signal count.\n')
  cat('No LD ridge adjustment, automatic allele flipping or residual-based SNP deletion.\n')
  config <- list(input=input, scenarios=lapply(masks,function(x) colnames(R)[x]),
    L=10L,max_iter=1000L,tol=1e-3,coverage=0.95,min_abs_corr=0.5,
    estimate_residual_variance=FALSE,refine=FALSE,reference_N=obj$reference_n,
    ABF_p1=1e-4,ABF_p2=1e-4,ABF_p12=c(1e-6,1e-5,5e-5),
    ABF_sdY='estimated from trait MAF, varbeta and study-level N')
  saveRDS(config,file.path(out,'analysis_config.rds'))
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  log_messages <- character(); fit_rows <- list(); abf_rows <- list()
  write_logs <- function() writeLines(if(length(log_messages)) log_messages else 'No captured warnings/errors.',file.path(out,'warnings_and_errors.txt'))
  capture_run <- function(expr,label) {
    ans <- tryCatch(withCallingHandlers(expr,warning=function(w) {
      log_messages <<- c(log_messages,paste(label,'WARNING',conditionMessage(w)))
      invokeRestart('muffleWarning')
    }),error=function(e) {
      log_messages <<- c(log_messages,paste(label,'ERROR',conditionMessage(e)))
      structure(list(message=conditionMessage(e)),class='audit_error')
    })
    write_logs()
    ans
  }
  write_table <- function(x,name) data.table::fwrite(x,file.path(out,name),sep='\t')
  pairs <- list(fat_iron=traits[c(1,2)],fat_volume=traits[c(1,3)],iron_volume=traits[c(2,3)])
  subset_data <- function(d,k) list(snp=d$snp[k],position=d$position[k],beta=d$beta[k],
    varbeta=d$varbeta[k],MAF=d$MAF[k],N=d$N,type='quant')
  for (scenario in names(masks)) {
    keep <- masks[[scenario]]
    Rs <- R[keep,keep,drop=FALSE]
    attr(Rs,'eigen') <- NULL
    ds <- lapply(obj$aligned[traits],subset_data,k=keep)
    prefix <- paste0(scenario,'__')
    cat('\n===== ',scenario,' n=',ncol(Rs),' =====\n',sep='')
    write_table(data.frame(SNP=colnames(Rs),reference_MAF=ref_maf[keep]),paste0(prefix,'variants.tsv'))
    # Recompute ABF on exactly the same variants as each SuSiE scenario.
    for (pn in names(pairs)) for (p12 in config$ABF_p12) {
      ts <- pairs[[pn]]; label <- paste(scenario,pn,p12,sep='__')
      fit <- capture_run(coloc::coloc.abf(ds[[ts[1]]],ds[[ts[2]]],
                       p1=1e-4,p2=1e-4,p12=p12),label)
      saveRDS(fit,file.path(out,paste0(label,'_abf.rds')))
      vals <- rep(NA_real_,5)
      if (!inherits(fit,'audit_error')) vals <- unname(fit$summary[paste0('PP.H',0:4,'.abf')])
      abf_rows[[label]] <- data.frame(scenario=scenario,pair=pn,nsnps=ncol(Rs),p12=p12,
        PP_H0=vals[1],PP_H1=vals[2],PP_H2=vals[3],PP_H3=vals[4],PP_H4=vals[5],
        status=if(inherits(fit,'audit_error')) 'ERROR' else 'completed')
      write_table(data.table::rbindlist(abf_rows),'ABF_matched_variants_summary.tsv')
    }
    for (trait in traits) {
      d <- ds[[trait]]; label <- paste(scenario,trait,sep='__')
      z <- d$beta/sqrt(d$varbeta); names(z) <- d$snp
      cat('\nFitting SuSiE:',label,'\n'); flush.console()
      set.seed(20261003)
      fit <- capture_run(susieR::susie_rss(z=z,R=Rs,n=d$N,L=10,
        estimate_residual_variance=FALSE,estimate_prior_variance=TRUE,
        coverage=0.95,min_abs_corr=0.5,n_purity=ncol(Rs),
        max_iter=1000,tol=1e-3,refine=FALSE,verbose=FALSE),label)
      saveRDS(fit,file.path(out,paste0(label,'_susie_fit.rds')))
      row <- data.frame(scenario=scenario,trait=trait,nsnps=length(z),N=d$N,
        status='ERROR',converged=FALSE,niter=NA_integer_,n_CS=NA_integer_,
        max_PIP=NA_real_,top_SNP=NA_character_,elbo_drop_count=NA_integer_)
      if (!inherits(fit,'audit_error')) {
        post <- capture_run({
          if (length(fit$pip)!=length(z) || any(!is.finite(fit$pip))) stop('Unexpected PIP output')
          sets <- fit$sets$cs
          ncs <- if(is.null(sets)) 0L else length(sets)
          pip <- data.frame(SNP=d$snp,position=d$position,PIP=as.numeric(fit$pip),
                            beta_A1=d$beta,se=sqrt(d$varbeta))
          write_table(pip[order(-pip$PIP),],paste0(label,'_PIP.tsv'))
          write_table(pip[match(keys,pip$SNP),],paste0(label,'_five_candidates.tsv'))
          csrows <- list(); members <- list()
          if(ncs>0L) for(j in seq_along(sets)) {
            ix <- sets[[j]]
            component <- fit$sets$cs_index[j]
            if(length(component)!=1L || is.na(component)) stop('Missing CS component index')
            cr <- abs(Rs[ix,ix,drop=FALSE])
            rvalues <- if(length(ix)==1L) 1 else cr[upper.tri(cr)]
            minr <- min(rvalues)
            csrows[[j]] <- data.frame(CS=names(sets)[j],component=component,size=length(ix),
              posterior_mass=sum(fit$alpha[component,ix]),min_abs_r=minr,
              mean_abs_r=mean(rvalues),median_abs_r=median(rvalues),
              effect_logBF=fit$lbf[component],
              lead_SNP=d$snp[ix[which.max(fit$alpha[component,ix])]])
            members[[j]] <- data.frame(CS=names(sets)[j],component=component,
              SNP=d$snp[ix],position=d$position[ix],alpha=fit$alpha[component,ix],PIP=fit$pip[ix])
          }
          cs_tab <- if(ncs) data.table::rbindlist(csrows) else data.frame(
            CS=character(),component=integer(),size=integer(),posterior_mass=numeric(),
            min_abs_r=numeric(),mean_abs_r=numeric(),median_abs_r=numeric(),
            effect_logBF=numeric(),lead_SNP=character())
          member_tab <- if(ncs) data.table::rbindlist(members) else data.frame(
            CS=character(),component=integer(),SNP=character(),position=integer(),alpha=numeric(),PIP=numeric())
          write_table(cs_tab,paste0(label,'_credible_sets.tsv'))
          write_table(member_tab,paste0(label,'_CS_members.tsv'))
          write_table(data.frame(iteration=seq_along(fit$elbo),ELBO=fit$elbo),paste0(label,'_ELBO.tsv'))
          list(ncs=ncs,top=d$snp[which.max(fit$pip)],maxpip=max(fit$pip))
        },paste0(label,'_postprocess'))
        row$converged <- isTRUE(fit$converged)
        row$niter <- if(length(fit$niter)==1L) fit$niter else NA_integer_
        row$elbo_drop_count <- sum(diff(fit$elbo)< -1e-6,na.rm=TRUE)
        if (!inherits(post,'audit_error')) {
          row$n_CS <- post$ncs; row$max_PIP <- post$maxpip; row$top_SNP <- post$top
          row$status <- if(row$converged) 'completed_review_required' else 'NOT_CONVERGED'
        } else row$status <- 'POSTPROCESS_ERROR_fit_saved'
      }
      fit_rows[[label]] <- row
      write_table(data.table::rbindlist(fit_rows),'SuSiE_fit_summary.tsv')
      print(row,row.names=FALSE)
    }
  }
  cat('\n===== SuSiE summary =====\n')
  print(as.data.frame(data.table::rbindlist(fit_rows)),row.names=FALSE)
  cat('\n===== ABF main prior p12=1e-5 =====\n')
  at <- as.data.frame(data.table::rbindlist(abf_rows))
  print(at[at$p12==1e-5,],row.names=FALSE)
  write_logs()
  cat('\nOutput:',out,'\n')
  cat('All scenarios are retained. Convergence is not proof of valid LD or causal signals.\n')
  cat('No coloc.susie inference performed before reviewing fit quality and credible sets.\n')
})()
