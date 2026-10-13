# 第19步：用户目录安装独立LDSC环境、验证EUR LD scores和HapMap3。
# 仅软件及资源准备；模拟数据测试绝不是正式遗传力结果；不读取用户GWAS。
(function() {
  if(Sys.info()[['sysname']]!='Linux'||!Sys.info()[['machine']] %in% c('x86_64','amd64'))
    stop('This installer is for Linux x86_64, as in your server audit.')
  if(!requireNamespace('data.table',quietly=TRUE)) stop('Missing R package: data.table')
  base <- path.expand('~/胰腺')
  parent <- file.path(base,'19_LDSC_setup');dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile('LDSC_setup_',tmpdir=parent);if(!dir.create(out)) stop('Cannot create output folder')
  stack <- path.expand('~/mr_tools/ldsc_stack19')
  downloads <- file.path(stack,'downloads');dir.create(downloads,recursive=TRUE,showWarnings=FALSE)
  refroot <- file.path(base,'ldsc_reference','EUR_tutorial_zenodo14993076')
  dir.create(refroot,recursive=TRUE,showWarnings=FALSE)
  complete <- FALSE;warnings_seen <- character()
  sink(file.path(out,'run_log.txt'),split=TRUE);on.exit(sink(),add=TRUE)
  oldopt <- options(timeout=1200);on.exit(options(oldopt),add=TRUE)
  on.exit({
    writeLines(if(complete) 'COMPLETED_SOFTWARE_AND_REFERENCES_READY' else 'INCOMPLETE_CHECK_LOGS',file.path(out,'completion_status.txt'))
    writeLines(if(length(warnings_seen)) warnings_seen else 'No captured R warnings/errors.',file.path(out,'warnings_and_errors.txt'))
    writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  },add=TRUE)
  wt <- function(x,n) data.table::fwrite(x,file.path(out,n),sep='\t',na='NA')
  cat('Output:',out,'\nSoftware:',stack,'\nReferences:',refroot,'\n')
  steps <- list();resources <- list()
  mark <- function(stage,status) {
    steps[[stage]] <<- data.frame(stage=stage,status=status,time=as.character(Sys.time()))
    wt(do.call(rbind,steps),'setup_status.tsv')
  }
  fetch <- function(url,dest,label,md5=NULL) {
    valid <- file.exists(dest)&&file.info(dest)$size>0
    if(valid&&!is.null(md5)) valid <- identical(unname(tools::md5sum(dest)),md5)
    if(!valid) {
      for(attempt in 1:3) {
        cat('Downloading:',label,' attempt',attempt,'\n');flush.console()
        part <- paste0(dest,'.part')
        ok <- tryCatch({
          rc <- download.file(url,part,mode='wb',method='libcurl',quiet=FALSE)
          if(rc!=0L||!file.exists(part)||file.info(part)$size<100) stop('Empty/incomplete download')
          if(!is.null(md5)&&unname(tools::md5sum(part))!=md5) stop('Published MD5 verification failed')
          if(!file.rename(part,dest)) stop('Cannot finalize download')
          TRUE
        },error=function(e) {cat(conditionMessage(e),'\n');FALSE})
        if(ok) break
        if(attempt==3) stop('Download failed: ',label,'; see run_log.txt. Do not disable certificate checks.')
        Sys.sleep(3)
      }
    }
    resources[[label]] <<- data.frame(resource=label,url=url,path=dest,bytes=file.info(dest)$size,
      md5=unname(tools::md5sum(dest)),published_MD5=if(is.null(md5)) NA_character_ else md5)
    wt(do.call(rbind,resources),'download_manifest.tsv')
  }
  run <- function(command,args,label) {
    log <- file.path(out,paste0(label,'.txt'))
    cat('Running:',label,'; progress is recorded in',basename(log),'\n');flush.console()
    rc <- system2(command,args=vapply(args,shQuote,character(1)),stdout=log,stderr=log,
      env=c('PYTHONNOUSERSITE=1','PYTHONPATH=','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1'))
    if(rc!=0L) {
      if(file.exists(log)) cat(tail(readLines(log,warn=FALSE),30L),sep='\n')
      stop('Command failed: ',label,' (exit ',rc,'). The full log is ',log)
    }
    invisible(log)
  }
  safe_untar <- function(archive,dest) {
    members <- utils::untar(archive,list=TRUE)
    if(!length(members)||any(grepl('(^/|(^|/)\\.\\.(/|$))',members,perl=TRUE))) stop('Unsafe or empty archive')
    rc <- utils::untar(archive,exdir=dest)
    if(rc!=0L) stop('Archive extraction failed: ',archive)
  }
  work <- function() {
    # 1. 固定micromamba版本。无需sudo，不更改系统Python或shell配置。
    mark('software_environment','RUNNING')
    ma <- file.path(downloads,'micromamba_2.3.2_linux64.tar.bz2')
    fetch('https://micro.mamba.pm/api/micromamba/linux-64/2.3.2',ma,'micromamba_2.3.2')
    bin <- file.path(stack,'bin','micromamba')
    if(!file.exists(bin)) {
      rc <- utils::untar(ma,files='bin/micromamba',exdir=stack)
      if(rc!=0L||!file.exists(bin)) stop('Cannot extract micromamba binary')
      Sys.chmod(bin,'0700')
    }
    vlog <- run(bin,'--version','micromamba_version')
    if(!any(grepl('2.3.2',readLines(vlog,warn=FALSE),fixed=TRUE))) stop('Unexpected micromamba version')
    envdir <- file.path(stack,'env_py39');py <- file.path(envdir,'bin','python')
    specs <- c('python=3.9','numpy=1.21.5','pandas=1.3.3','scipy=1.7.3','bitarray=2.9.2')
    writeLines(specs,file.path(out,'requested_environment_specs.txt'))
    # 数值核心版本沿用CBIIT environment3.yml；不使用该分支遗留的Python2 requirements.txt。
    common_args <- c('--no-rc','--root-prefix',file.path(stack,'mamba_root'))
    if(!file.exists(py)) {
      if(dir.exists(envdir)&&length(list.files(envdir,all.files=TRUE,no..=TRUE))>0L &&
         !file.exists(file.path(envdir,'conda-meta','history')))
        stop('An incomplete environment exists at ',envdir,'. Send the log; do not delete unrelated environments.')
      run(bin,c(common_args,'create','--yes','--prefix',envdir,'--override-channels','--channel','conda-forge',specs),
          'environment_install')
    }
    if(!file.exists(py)) stop('Environment python not found')
    run(bin,c(common_args,'list','--prefix',envdir,'--explicit'),'environment_explicit_lock')
    check <- paste('import sys,numpy,pandas,scipy,bitarray;',
      'assert sys.version_info[:2]==(3,9);',
      'assert numpy.__version__=="1.21.5";',
      'assert pandas.__version__=="1.3.3";',
      'assert scipy.__version__=="1.7.3";',
      'assert bitarray.__version__=="2.9.2";',
      'print(sys.version);print("Pinned dependency checks passed")')
    run(py,c('-c',check),'dependency_check')
    mark('software_environment','PASS')
    # 2. 固定源码commit，不运行git pull，不修改源码。
    mark('LDSC_source','RUNNING')
    commit <- '6c673952cee74bd5c57aef1555a03b1c015399a0'
    src_archive <- file.path(downloads,paste0('CBIIT_ldsc_',commit,'.tar.gz'))
    fetch(paste0('https://codeload.github.com/CBIIT/ldsc/tar.gz/',commit),src_archive,'CBIIT_LDSC_fixed_commit')
    src_parent <- file.path(stack,'source');dir.create(src_parent,showWarnings=FALSE)
    repo <- file.path(src_parent,paste0('ldsc-',commit))
    if(!file.exists(file.path(repo,'ldsc.py'))) safe_untar(src_archive,src_parent)
    if(!all(file.exists(file.path(repo,c('ldsc.py','munge_sumstats.py','environment3.yml'))))) stop('Source archive incomplete')
    run(py,c(file.path(repo,'ldsc.py'),'--help'),'ldsc_help')
    run(py,c(file.path(repo,'munge_sumstats.py'),'--help'),'munge_help')
    wt(data.frame(repository='https://github.com/CBIIT/ldsc',branch_origin='ldsc39',commit=commit,
      python=py,ldsc=file.path(repo,'ldsc.py'),munge=file.path(repo,'munge_sumstats.py')),'software_manifest.tsv')
    sourcefiles <- list.files(repo,pattern='\\.py$',recursive=TRUE,full.names=TRUE)
    wt(data.frame(path=sourcefiles,md5=unname(tools::md5sum(sourcefiles))),'source_code_checksums.tsv')
    mark('LDSC_source','PASS')
    # 3. 旧Broad下载链接失效；以下为注明原资源来源的Zenodo存档，校验发布MD5。
    mark('reference_download','RUNNING')
    ar <- file.path(downloads,'eur_w_ld_chr_zenodo14993076.tgz')
    fetch('https://zenodo.org/records/14993076/files/eur_w_ld_chr.tgz?download=1',ar,
      'EUR_LD_scores_Zenodo_14993076','76c1890c8cf22d99d05c6707cc8441b4')
    refdir <- file.path(refroot,'eur_w_ld_chr')
    expected <- unlist(lapply(1:22,function(ch) file.path(refdir,paste0(ch,c('.l2.ldscore.gz','.l2.M_5_50')))))
    if(!all(file.exists(expected))) safe_untar(ar,refroot)
    if(!all(file.exists(expected))) stop('Reference archive layout differs; inspect extracted files before proceeding')
    hg <- file.path(downloads,'w_hm3_zenodo7773502.snplist.gz')
    fetch('https://zenodo.org/records/7773502/files/w_hm3.snplist.gz?download=1',hg,
      'HapMap3_Zenodo_7773502','153ecc2bcfa740afafe656e6a384d769')
    hm <- file.path(refroot,'w_hm3.snplist')
    # 始终从校验过的gzip重新生成专用参考文本。
    ci <- gzfile(hg,'rb');co <- file(hm,'wb')
    tryCatch({repeat {b <- readBin(ci,'raw',n=1048576L);if(!length(b)) break;writeBin(b,co)}},
      finally={close(ci);close(co)})
    mark('reference_download','PASS')
    # 4. 全部参考文件结构及校验值核验，并进行小规模模拟数据软件测试。
    mark('reference_QC_and_software_test','RUNNING')

validator_lines <- c(
"# Generated by step 19. Software test only; no user GWAS is analysed.",
"import sys, os, json, gzip, hashlib, subprocess",
"from pathlib import Path",
"import numpy as np",
"import pandas as pd",
"from scipy.stats import norm",
"out, repo, ref, hmfile = map(Path, sys.argv[1:5])",
"def sha(p):",
"    h=hashlib.sha256()",
"    with open(p,'rb') as f:",
"        for b in iter(lambda:f.read(1048576),b''): h.update(b)",
"    return h.hexdigest()",
"def run(args, filename):",
"    p=subprocess.run([str(x) for x in args],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)",
"    (out/filename).write_text(p.stdout,encoding='utf-8')",
"    if p.returncode or 'Traceback (most recent call last)' in p.stdout:",
"        raise RuntimeError('Software test failed: '+filename+'\\n'+p.stdout[-2000:])",
"    return p.stdout",
"hm=pd.read_csv(hmfile,sep=r'\\s+',dtype=str)",
"assert list(hm.columns)==['SNP','A1','A2'], 'Unexpected HapMap3 columns'",
"assert len(hm)==1217311, 'Unexpected HapMap3 SNP count'",
"assert hm.SNP.is_unique and hm.notna().all().all(), 'Duplicate or missing HapMap3 entries'",
"assert hm.A1.isin(list('ACGT')).all() and hm.A2.isin(list('ACGT')).all()",
"assert (hm.A1!=hm.A2).all()",
"checks=[]; fingerprints=[]; samples=[]; seen=set()",
"rng=np.random.default_rng(190019)",
"for chrom in range(1,23):",
"    p=ref/(str(chrom)+'.l2.ldscore.gz')",
"    d=pd.read_csv(p,sep=r'\\s+')",
"    assert {'CHR','SNP','BP','L2'}.issubset(d.columns), str(p)",
"    assert len(d)>1000 and d.SNP.is_unique and not seen.intersection(d.SNP), str(p)",
"    seen.update(d.SNP)",
"    assert (d.CHR==chrom).all() and (d.BP>0).all() and d.BP.is_monotonic_increasing",
"    assert np.isfinite(d.L2).all(), 'Nonfinite LD scores'",
"    mm=ref/(str(chrom)+'.l2.M');mf=ref/(str(chrom)+'.l2.M_5_50')",
"    m5=np.atleast_1d(np.loadtxt(mf))",
"    assert m5.size==1 and np.isfinite(m5[0]) and m5[0]>0",
"    merged=d.merge(hm,on='SNP',how='inner',validate='one_to_one')",
"    checks.append(dict(CHR=chrom,n_LD_SNP=len(d),n_HapMap3_overlap=len(merged),",
"                       min_LD=float(d.L2.min()),max_LD=float(d.L2.max()),M_5_50=float(m5[0])))",
"    for f in [p,mf]: fingerprints.append(dict(path=str(f),sha256=sha(f),bytes=f.stat().st_size))",
"    assert len(merged)>1000",
"    pal=(merged.A1+merged.A2).isin(['AT','TA','CG','GC'])",
"    merged=merged.loc[~pal]",
"    idx=rng.choice(len(merged),size=min(1500,len(merged)),replace=False)",
"    samples.append(merged.iloc[np.sort(idx)].copy())",
"    print('Validated reference chromosome',chrom,flush=True)",
"pd.DataFrame(checks).to_csv(out/'reference_QC.tsv',sep='\\t',index=False)",
"fingerprints.append(dict(path=str(hmfile),sha256=sha(hmfile),bytes=hmfile.stat().st_size))",
"pd.DataFrame(fingerprints).to_csv(out/'reference_file_checksums.tsv',sep='\\t',index=False)",
"(out/'HapMap3_QC.json').write_text(json.dumps(dict(rows=len(hm),unique_SNP=int(hm.SNP.nunique()),",
"    columns=list(hm.columns),reference_overlap=len(seen.intersection(set(hm.SNP)))),indent=2))",
"# Independent synthetic test statistics: these are not patient or study data.",
"simdir=out/'software_test_only';simdir.mkdir(exist_ok=True)",
"s=pd.concat(samples,ignore_index=True)",
"s['Z']=rng.normal(size=len(s))*np.sqrt(1+0.01*np.maximum(s.L2.to_numpy(),0))",
"s['P']=2*norm.sf(np.abs(s.Z));s['N']=30000",
"raw=simdir/'SYNTHETIC_ONLY.tsv';prefix=simdir/'SYNTHETIC_ONLY'",
"s[['SNP','A1','A2','Z','P','N']].to_csv(raw,sep='\\t',index=False)",
"run([sys.executable,repo/'munge_sumstats.py','--sumstats',raw,'--out',prefix,",
"     '--merge-alleles',hmfile,'--signed-sumstats','Z,0'], 'software_test_munge_console.txt')",
"sumfile=Path(str(prefix)+'.sumstats.gz')",
"assert sumfile.is_file()",
"v=pd.read_csv(sumfile,sep=r'\\s+').dropna(subset=['Z','N'])",
"assert len(v)>1000 and v.SNP.is_unique",
"hprefix=simdir/'SYNTHETIC_ONLY_h2'",
"run([sys.executable,repo/'ldsc.py','--h2',sumfile,'--ref-ld-chr',str(ref)+os.sep,",
"     '--w-ld-chr',str(ref)+os.sep,'--out',hprefix], 'software_test_h2_console.txt')",
"log=Path(str(hprefix)+'.log').read_text()",
"assert 'Total Observed scale h2' in log and 'Traceback' not in log",
"# NaN/Inf estimates or standard errors are not accepted as a working smoke test.",
"import re",
"match=re.search(r'Total Observed scale h2:\\s*([-+0-9.eE]+)\\s*\\(([-+0-9.eE]+)\\)',log)",
"assert match and np.isfinite([float(match[1]),float(match[2])]).all(), 'Nonfinite software test estimate'",
"(out/'software_test_status.tsv').write_text('test\\tstatus\\nHapMap3_structure\\tPASS\\nEUR_22_chromosomes\\tPASS\\nmunge_synthetic\\tPASS\\nLDSC_h2_synthetic\\tPASS\\n')",
"# Avoid confusing simulated estimates with scientific results.",
"(simdir/'NOT_RESEARCH_RESULTS.txt').write_text('All data and estimates here are synthetic SOFTWARE TESTS ONLY. Never use them in a manuscript or as GWAS results.\\n')",
"versions={'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__}",
"import scipy,bitarray",
"versions.update(scipy=scipy.__version__,bitarray=bitarray.__version__)",
"(out/'python_package_versions.json').write_text(json.dumps(versions,indent=2))",
"print('Reference validation and synthetic software tests PASSED.',flush=True)"
)
    validator <- file.path(out,'validate_reference_and_test.py')
    writeLines(validator_lines,validator,useBytes=TRUE)
    run(py,c(validator,out,repo,refdir,hm),'reference_validation_console')
    if(!file.exists(file.path(out,'software_test_status.tsv'))) stop('Software test summary missing')
    mark('reference_QC_and_software_test','PASS')
    config <- list(python=py,ldsc=file.path(repo,'ldsc.py'),munge_sumstats=file.path(repo,'munge_sumstats.py'),
      source_commit=commit,reference_ld_chr=paste0(refdir,'/'),weight_ld_chr=paste0(refdir,'/'),
      hapmap3_snplist=hm,reference_population='1000 Genomes European',
      reference_source='Zenodo 14993076: archived original tutorial EUR scores',
      hapmap3_source='Zenodo 7773502: archived Alkes Group HapMap3 list',
      analysis_scope='Unpartitioned LDSC preparation only; no study GWAS analysed')
    saveRDS(config,file.path(out,'LDSC_paths.rds'))
    wt(data.frame(item=names(config),value=unlist(config,use.names=FALSE)),'LDSC_paths.tsv')
    writeLines(c('Step 19 prepares software and reference data only.',
      'The installed implementation is CBIIT LDSC, fixed commit 6c673952cee74bd5c57aef1555a03b1c015399a0; do not silently update it.',
      'Python 3.9 and core numerical package versions follow the compatible environment3.yml; legacy requirements.txt is not used.',
      'Reference archive and HapMap3 compressed file were checked against the MD5 values published in their Zenodo records.',
      'The current numeric regression test is synthetic SOFTWARE VALIDATION ONLY, never a scientific h2 result.',
      'For unpartitioned LDSC, these tutorial LD scores are used for both reference and regression weights.',
      'Do not reuse the previous 503-person local LD matrix as genome-wide LD scores.',
      'Next: harmonise GWAS with HapMap3 and inspect duplicates, sample size assumptions, and allele compatibility.',
      'Later: estimate single-trait h2 and its precision before interpreting rg; retain the bivariate intercept by default for overlapping UKB samples.'),
      file.path(out,'README_next_step.txt'))
    script <- file.path(base,'19_setup_LDSC_resources.R')
    if(file.exists(script)) file.copy(script,file.path(out,basename(script)))
    complete <<- TRUE
    cat('\nCOMPLETED_SOFTWARE_AND_REFERENCES_READY\nOutput:',out,'\n')
    cat('Upload this output folder only. No need to upload the installed environment or large reference files.\n')
  }
  tryCatch(withCallingHandlers(work(),warning=function(w) {
    warnings_seen <<- c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')
  }),error=function(e) {
    warnings_seen <<- c(warnings_seen,paste('ERROR:',conditionMessage(e)))
    stop(conditionMessage(e),call.=FALSE)
  })
})()
