# 第22步：8表型28组LDSC遗传相关，自由截距，固定28项BH-FDR。
(function() {
  base <- path.expand('~/胰腺')
  prior <- file.path(base,'19_LDSC_setup','LDSC_setup_1566e430f1187')
  inputs <- file.path(base,'20_LDSC_inputs','GWAS_inputs_16c90c74e8c0e1')
  h2dir <- file.path(base,'21_LDSC_heritability','h2_16c90c661ce88d')
  check_status <- function(path,expected) {
    if(!file.exists(path)||!identical(readLines(path,warn=FALSE),expected)) stop('完成标记不一致或文件缺失：',path)
  }
  check_status(file.path(prior,'completion_status.txt'),'COMPLETED_SOFTWARE_AND_REFERENCES_READY')
  check_status(file.path(inputs,'reports','completion_status.txt'),'COMPLETED_INPUT_PREPARATION_REVIEW_QC_BEFORE_LDSC')
  check_status(file.path(h2dir,'reports','completion_status.txt'),'COMPLETED_8_TRAIT_H2_REVIEW_ESTIMATES')
  if(!requireNamespace('data.table',quietly=TRUE)) stop('缺少data.table。')
  cfg <- readRDS(file.path(prior,'LDSC_paths.rds'))
  if(!file.exists(cfg$python)) stop('第19步Python不存在。')
  parent <- file.path(base,'22_LDSC_genetic_correlation');dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile('rg_',tmpdir=parent);if(!dir.create(out)) stop('无法建立结果目录。')
  reports <- file.path(out,'reports');dir.create(reports)
  writeLines('INCOMPLETE_CHECK_LOGS',file.path(reports,'completion_status.txt'))
  writeLines(capture.output(sessionInfo()),file.path(reports,'R_sessionInfo.txt'))
  writeLines(c(
    'All 28 unordered pairs of eight quantitative traits are attempted once; no selection on significance.',
    'Original CBIIT LDSC estimate_rg, fixed commit 6c673952cee74bd5c57aef1555a03b1c015399a0.',
    'Unconstrained h2 and genetic-covariance intercepts; allele checking enabled; 200 jackknife blocks; two-step cutoff 30.',
    'No --no-intercept, --intercept-gencov, --return-silly-things, --no-check-alleles, or additional chisq-max.',
    'Full-precision original object estimates are saved; pairwise h2 is re-estimated on each SNP intersection.',
    'Step21 low precision flag and pairwise h2 Z<4 flag are retained; neither automatically removes a trait.',
    'The original implementation suppresses its printed rg summary when abs(rg)>1.2. These values remain in *_internal columns only.',
    'Nonfinite estimates, nonpositive SE, nonpositive pair h2, original out-of-bounds warning and unexpected SNP counts require review.',
    'Absolute rg between 1 and 1.2 is retained with a flag, not truncated to 1. Wald intervals are untruncated.',
    'BH-FDR uses stats::p.adjust(p, method="BH", n=28). The family never shrinks to successful or significant tests.',
    'Missing/non-interpretable tests remain NA; with missing tests adjustment is provisional using a fixed family of 28.',
    'Pancreas-fat seven-pair subset uses the SAME 28-test FDR, not a separate seven-test correction.',
    'FDR<0.05 indicates statistical evidence under model assumptions, not causality or a mechanistic pathway.',
    'Study-level N and unavailable per-SNP INFO remain limitations. Free intercept does not resolve every possible bias.',
    'Heatmap colours saturate at +/-1; displayed estimates and tables are never truncated. Diagonal 1 is identity, not an extra fit.',
    'Official tutorial: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation'),
    file.path(reports,'methods_and_interpretation.txt'))
  script <- file.path(base,'22_run_LDSC_genetic_correlation.R')
  if(file.exists(script)) file.copy(script,file.path(reports,basename(script)))
  python_lines <- c(
"\"\"\"Step21: call the unchanged fixed LDSC estimator; save full-precision Hsq fields.\"\"\"",
"import sys, json, re, hashlib, traceback, time, importlib.util, gc, platform",
"from pathlib import Path",
"import numpy as np",
"import pandas as pd",
"",
"EXPECTED = {",
" 'ebi-a-GCST90016666':('liver_volume',32860),",
" 'ebi-a-GCST90016669':('pancreas_volume',31758),",
" 'ebi-a-GCST90016671':('visceral_adipose_tissue_volume',32860),",
" 'ebi-a-GCST90016672':('abdominal_subcutaneous_adipose_tissue_volume',32860),",
" 'ebi-a-GCST90016673':('liver_fat',32858),",
" 'ebi-a-GCST90016674':('liver_iron',32858),",
" 'ebi-a-GCST90016675':('pancreas_fat',25617),",
" 'ebi-a-GCST90016676':('pancreas_iron',25617)}",
"",
"def digest(path,kind='sha256'):",
"    h=hashlib.new(kind)",
"    with open(path,'rb') as f:",
"        for b in iter(lambda:f.read(2**20),b''):h.update(b)",
"    return h.hexdigest()",
"",
"def scalar(value):",
"    try:",
"        a=np.asarray(value,dtype=float)",
"        if a.size!=1:raise ValueError('Expected scalar')",
"        v=float(a.reshape(-1)[0])",
"        return v if np.isfinite(v) else None",
"    except (TypeError,ValueError):return None",
"",
"from itertools import combinations",
"class Log:",
"    def __init__(self,p):self.f=open(p,'w',encoding='utf-8')",
"    def log(self,s):self.f.write(str(s)+'\\n');self.f.flush();print(s,flush=True)",
"    def close(self):self.f.close()",
"",
"def write_table(rows,path):pd.DataFrame(rows).to_csv(path,sep='\\t',index=False,float_format='%.17g',na_rep='NA')",
"",
"def extract_rg(fit):",
"    row={}",
"    for name,attr in [('rg_internal','rg_ratio'),('rg_se_internal','rg_se'),('z_internal','z'),('p_internal','p')]:",
"        row[name]=scalar(getattr(fit,attr,None))",
"    flags=[]",
"    for k,obj in [('h2_1',fit.hsq1),('h2_2',fit.hsq2),('gencov',fit.gencov)]:",
"        row[k]=scalar(obj.tot);row[k+'_se']=scalar(obj.tot_se)",
"        row[k+'_intercept']=scalar(obj.intercept);row[k+'_intercept_se']=scalar(obj.intercept_se)",
"        if k!='gencov':",
"            v,se=row[k],row[k+'_se']",
"            row[k+'_z']=v/se if v is not None and se is not None and se>0 else None",
"            if row[k+'_z'] is None or row[k+'_z']<4:flags.append(k+'_LOW_PRECISION')",
"    valid=all(row[k] is not None for k in ['rg_internal','rg_se_internal','z_internal','p_internal'])",
"    valid=valid and row['rg_se_internal']>0 and 0<=row['p_internal']<=1",
"    if valid and not np.isclose(row['rg_internal']/row['rg_se_internal'],row['z_internal'],rtol=1e-8,atol=1e-10):",
"        valid=False;flags.append('Z_INCONSISTENCY')",
"    if getattr(fit,'_negative_hsq',False):valid=False;flags.append('NONPOSITIVE_PAIR_H2')",
"    # Preserve internal values but respect the original fixed implementation's warning threshold.",
"    if row['rg_internal'] is not None and abs(row['rg_internal'])>1.2:",
"        valid=False;flags.append('LDSC_RG_OUT_OF_BOUNDS')",
"    elif row['rg_internal'] is not None and abs(row['rg_internal'])>1:",
"        flags.append('ABS_RG_GT_1_REVIEW')",
"    row['rg']=row['rg_internal'] if valid else None",
"    row['rg_se']=row['rg_se_internal'] if valid else None",
"    row['z']=row['z_internal'] if valid else None",
"    row['p']=row['p_internal'] if valid else None",
"    row['ci_lower']=row['rg']-1.959963984540054*row['rg_se'] if valid else None",
"    row['ci_upper']=row['rg']+1.959963984540054*row['rg_se'] if valid else None",
"    row['estimate_status']='ESTIMABLE' if valid else 'NOT_INTERPRETABLE_REVIEW'",
"    row['pair_QC_flag']=';'.join(flags) if flags else 'NO_AUTOMATIC_FLAG'",
"    return row",
"",
"def main():",
"    prior=Path(sys.argv[1]);inputs=Path(sys.argv[2]);out=Path(sys.argv[3]);h2dir=Path(sys.argv[4])",
"    report=out/'reports';logs=out/'ldsc_logs';report.mkdir(exist_ok=True);logs.mkdir(exist_ok=True)",
"    status=report/'completion_status.txt';status.write_text('INCOMPLETE_CHECK_LOGS\\n')",
"    cfg=pd.read_csv(prior/'LDSC_paths.tsv',sep='\\t').set_index('item').value.to_dict()",
"    if cfg['source_commit']!='6c673952cee74bd5c57aef1555a03b1c015399a0':raise ValueError('Wrong LDSC commit')",
"    manifest=pd.read_csv(inputs/'reports'/'GWAS_LDSC_QC_summary.tsv',sep='\\t')",
"    if len(manifest)!=8 or not manifest.id.is_unique or set(manifest.id)!=set(EXPECTED):raise ValueError('Wrong study set')",
"    checks=[]",
"    print('Checking step19 resources and step20 input hashes...',flush=True)",
"    for filename,col,kind in [('reference_file_checksums.tsv','sha256','sha256'),('source_code_checksums.tsv','md5','md5')]:",
"        for row in pd.read_csv(prior/filename,sep='\\t').to_dict('records'):",
"            p=Path(row['path']);ok=p.is_file() and digest(p,kind)==row[col]",
"            checks.append(dict(kind=kind,path=str(p),match=ok))",
"            if not ok:raise ValueError('Resource checksum mismatch: '+str(p))",
"    for row in manifest.to_dict('records'):",
"        sid=row['id'];trait,n=EXPECTED[sid]",
"        if row['trait']!=trait or row['N']!=n or row['status']!='PASS_INPUT_QC':raise ValueError('Input QC/trait mismatch: '+sid)",
"        p=inputs/'sumstats'/(sid+'.sumstats.gz')",
"        if p.resolve()!=Path(row['sumstats']).resolve():raise ValueError('Input path differs from step20')",
"        ok=p.is_file() and digest(p)==row['sumstats_sha256']",
"        checks.append(dict(kind='input_sha256',path=str(p),match=ok))",
"        if not ok:raise ValueError('Input checksum mismatch: '+sid)",
"    write_table(checks,report/'input_and_resource_checks.tsv')",
"    manifest.to_csv(report/'step20_input_manifest.tsv',sep='\\t',index=False)",
"    sys.path.insert(0,str(Path(cfg['ldsc']).parent))",
"    spec=importlib.util.spec_from_file_location('fixed_ldsc21',cfg['ldsc'])",
"    ldsc=importlib.util.module_from_spec(spec);spec.loader.exec_module(ldsc)",
"    import scipy,bitarray",
"    versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,",
"                  scipy=scipy.__version__,bitarray=bitarray.__version__,LDSC=ldsc.__version__,commit=cfg['source_commit'])",
"    (report/'software_versions.json').write_text(json.dumps(versions,indent=2))",
"    if sys.version_info[:2]!=(3,9) or np.__version__!='1.21.5' or pd.__version__!='1.3.3' or scipy.__version__!='1.7.3' or bitarray.__version__!='2.9.2':",
"        raise ValueError('Environment differs from validated step19')",
"    h2=pd.read_csv(h2dir/'reports'/'LDSC_h2_summary.tsv',sep='\\t').set_index('id')",
"    if len(h2)!=8 or not h2.index.is_unique or set(h2.index)!=set(EXPECTED) or not (h2.run_status=='COMPLETED').all():",
"        raise ValueError('Step21 h2 table incomplete')",
"    h2.to_csv(report/'step21_h2_summary.tsv',sep='\\t')",
"    overlap=pd.read_csv(inputs/'reports'/'pairwise_SNP_overlap.tsv',sep='\\t')",
"    expected_overlap={frozenset([r['trait1'],r['trait2']]):int(r['n_common']) for r in overlap.to_dict('records')}",
"    if len(expected_overlap)!=28:raise ValueError('Expected 28 SNP overlap comparisons')",
"    results=[];commands=[]",
"    pairs=list(combinations(EXPECTED,2))",
"    write_table([dict(pair_id='pair_%02d'%(i+1),id1=a,id2=b,trait1=EXPECTED[a][0],trait2=EXPECTED[b][0])",
"                 for i,(a,b) in enumerate(pairs)],report/'planned_28_pairs.tsv')",
"    for i,(a,b) in enumerate(pairs,1):",
"        pair_id='pair_%02d'%i;ta,na=EXPECTED[a];tb,nb=EXPECTED[b]",
"        prefix=logs/(pair_id+'_');logpath=logs/(pair_id+'.log')",
"        argv=['--rg',str(inputs/'sumstats'/(a+'.sumstats.gz'))+','+str(inputs/'sumstats'/(b+'.sumstats.gz')),",
"              '--ref-ld-chr',cfg['reference_ld_chr'],'--w-ld-chr',cfg['weight_ld_chr'],'--out',str(prefix),",
"              '--n-blocks','200','--two-step','30','--print-delete-vals']",
"        args=ldsc.parser.parse_args(argv)",
"        if args.no_intercept or args.intercept_h2 is not None or args.intercept_gencov is not None or args.no_check_alleles or args.return_silly_things:",
"            raise ValueError('Unexpected intercept/allele/rg setting')",
"        commands.append(dict(pair_id=pair_id,argv=[cfg['python'],cfg['ldsc']]+argv,entrypoint='ldsc.sumstats.estimate_rg(args, log)'))",
"        (report/'commands.json').write_text(json.dumps(commands,indent=2,ensure_ascii=False))",
"        row=dict(pair_id=pair_id,id1=a,id2=b,trait1=ta,trait2=tb,N1=na,N2=nb,",
"                 step21_h2_z1=float(h2.loc[a,'h2_z']),step21_h2_z2=float(h2.loc[b,'h2_z']),",
"                 involves_liver_iron='liver_iron' in (ta,tb),involves_pancreas_fat='pancreas_fat' in (ta,tb),",
"                 run_status='FAILED',estimate_status='NOT_ESTIMATED',rg=None,rg_se=None,p=None,",
"                 ci_lower=None,ci_upper=None,pair_QC_flag='NOT_EVALUATED')",
"        row['step21_low_precision']=row['step21_h2_z1']<4 or row['step21_h2_z2']<4",
"        log=Log(logpath);start=time.time()",
"        try:",
"            log.log(ldsc.MASTHEAD);log.log('Pair %d/28: %s vs %s'%(i,ta,tb));log.log(json.dumps(argv,ensure_ascii=False))",
"            fits=ldsc.sumstats.estimate_rg(args,log)",
"            text=logpath.read_text()",
"            if len(fits)!=1 or fits[0] is None or 'ERROR computing rg' in text or 'Traceback (most recent call last)' in text:",
"                raise ValueError('Original LDSC reported an internal pair failure; see log')",
"            fit=fits[0];row.update(extract_rg(fit))",
"            n=re.findall(r'([0-9,]+) SNPs with valid alleles\\.',text)",
"            if not n:raise ValueError('Cannot confirm pair SNP count')",
"            row['n_regression']=int(n[-1].replace(',',''))",
"            row['n_expected_overlap']=expected_overlap[frozenset([ta,tb])]",
"            if row['n_regression']!=row['n_expected_overlap']:",
"                row['pair_QC_flag']+=';SNP_COUNT_CHANGED_REVIEW'",
"                row['estimate_status']='NOT_INTERPRETABLE_REVIEW'",
"                for k in ['rg','rg_se','z','p','ci_lower','ci_upper']:row[k]=None",
"            row['n_jackknife_blocks']=int(np.asarray(fit.gencov.tot_delete_values).shape[0])",
"            row['run_status']='COMPLETED';row['elapsed_seconds']=time.time()-start",
"            log.log('Estimate status: '+row['estimate_status']+'; QC: '+row['pair_QC_flag'])",
"            del fit,fits;gc.collect()",
"        except Exception:",
"            row['error']=traceback.format_exc();log.log(row['error'])",
"        finally:log.close()",
"        (report/(pair_id+'_full_precision.json')).write_text(json.dumps(row,indent=2,ensure_ascii=False,allow_nan=False))",
"        results.append(row);write_table(results,report/'LDSC_rg_raw_summary.tsv')",
"        print('FINISHED',i,'/28',ta,tb,row['run_status'],flush=True)",
"    messages=[]",
"    for p in sorted(logs.glob('*.log')):",
"        for line in p.read_text().splitlines():",
"            if any(k in line.lower() for k in ['warning','traceback','error','ratio <']):messages.append(dict(file=p.name,message=line))",
"    pd.DataFrame(messages,columns=['file','message']).to_csv(report/'LDSC_log_messages.tsv',sep='\\t',index=False)",
"    success=all(r['run_status']=='COMPLETED' for r in results)",
"    status.write_text('COMPLETED_28_PAIR_FITS_REVIEW_ESTIMATES\\n' if success else 'PARTIAL_FAILURE_CHECK_LOGS\\n')",
"    if not success:raise RuntimeError('Some pairs failed; all 28 attempts retained. Upload review ZIP.')",
"",
"if __name__=='__main__':",
"    try:main()",
"    except Exception:",
"        out=Path(sys.argv[3]);(out/'reports').mkdir(exist_ok=True)",
"        (out/'reports'/'python_error.txt').write_text(traceback.format_exc())",
"        raise"
)
  helper <- file.path(out,'run_rg.py');writeLines(python_lines,helper,useBytes=TRUE)
  log <- file.path(reports,'python_console.txt')
  cat('结果目录：',out,'\n将顺序计算28组。进度日志：',log,'\n');flush.console()
  err <- NULL
  rc <- tryCatch(system2(cfg$python,args=vapply(c('-u',helper,prior,inputs,out,h2dir),shQuote,character(1)),
    stdout=log,stderr=log,env=c('PYTHONNOUSERSITE=1','PYTHONPATH=','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1')),
    error=function(e) {err <<- conditionMessage(e);127L})
  if(!is.null(err)) writeLines(err,file.path(reports,'R_error.txt'))
  post_error <- NULL
  tryCatch({
    raw <- file.path(reports,'LDSC_rg_raw_summary.tsv')
    if(file.exists(raw)) {
      d <- as.data.frame(data.table::fread(raw,na.strings='NA'))
      for(k in c('p','rg','rg_se','ci_lower','ci_upper')) d[[k]] <- as.numeric(d[[k]])
      valid <- d$run_status=='COMPLETED' & d$estimate_status=='ESTIMABLE' & is.finite(d$p) & d$p>=0 & d$p<=1
      d$FDR_BH_28 <- NA_real_
      if(any(valid)) d$FDR_BH_28[valid] <- p.adjust(d$p[valid],method='BH',n=28L)
      d$FDR_significant <- ifelse(is.na(d$FDR_BH_28),NA,d$FDR_BH_28<0.05)
      complete <- nrow(d)==28L && all(valid)
      d$FDR_scope_status <- if(complete) 'ALL_28_ESTIMABLE' else 'PROVISIONAL_FIXED_28_REVIEW_MISSING'
      data.table::fwrite(d,file.path(reports,'LDSC_rg_summary_FDR.tsv'),sep='\t',na='NA')
      sub <- d[d$trait1=='pancreas_fat'|d$trait2=='pancreas_fat',,drop=FALSE]
      sub$other_trait <- ifelse(sub$trait1=='pancreas_fat',sub$trait2,sub$trait1)
      data.table::fwrite(sub,file.path(reports,'pancreas_fat_rg_FDR.tsv'),sep='\t',na='NA')
      ids <- c('liver_volume','pancreas_volume','visceral_adipose_tissue_volume',
        'abdominal_subcutaneous_adipose_tissue_volume','liver_fat','liver_iron','pancreas_fat','pancreas_iron')
      labels <- c('Liver volume','Pancreas volume','VAT volume','ASAT volume','Liver fat','Liver iron','Pancreas fat','Pancreas iron')
      mat <- matrix(NA_real_,8,8,dimnames=list(ids,ids));diag(mat) <- 1
      qmat <- matrix(NA_real_,8,8,dimnames=list(ids,ids))
      for(i in seq_len(nrow(d))) {
        a <- match(d$trait1[i],ids);b <- match(d$trait2[i],ids)
        mat[a,b] <- mat[b,a] <- d$rg[i];qmat[a,b] <- qmat[b,a] <- d$FDR_BH_28[i]
      }
      data.table::fwrite(data.frame(trait=ids,mat,check.names=FALSE),file.path(reports,'rg_matrix.tsv'),sep='\t',na='NA')
      data.table::fwrite(data.frame(trait=ids,qmat,check.names=FALSE),file.path(reports,'FDR_matrix.tsv'),sep='\t',na='NA')
      # 图为结果核查用途，显著性星号使用统一的28项FDR。
      tryCatch({
        local({
          pdf(file.path(reports,'rg_heatmap.pdf'),width=9,height=9);on.exit(dev.off())
          par(mar=c(9,9,4,2));plot.new();plot.window(xlim=c(.5,8.5),ylim=c(.5,8.5),asp=1)
          pal <- colorRampPalette(c('#2166AC','white','#B2182B'))(201)
          for(a in 1:8) for(b in 1:8) {
            v <- mat[a,b];color <- if(is.finite(v)) pal[1+round((max(-1,min(1,v))+1)*100)] else 'grey90'
            rect(b-.5,8-a+.5,b+.5,8-a+1.5,col=color,border='white')
            lab <- if(is.finite(v)) sprintf('%.2f',v) else 'NA'
            if(is.finite(qmat[a,b])&&qmat[a,b]<.05) lab <- paste0(lab,'*')
            text(b,9-a,lab,cex=.78,col=if(is.finite(v)&&abs(v)>.7) 'white' else 'black')
          }
          axis(1,at=1:8,labels=labels,las=2,cex.axis=.8);axis(2,at=8:1,labels=labels,las=1,cex.axis=.8)
          title('LDSC genetic correlations');mtext('* BH-FDR <0.05 across 28 tests; blue negative, red positive',side=3,line=.2,cex=.7)
        })
        local({
          s <- sub[is.finite(sub$rg)&is.finite(sub$rg_se),,drop=FALSE]
          pdf(file.path(reports,'pancreas_fat_rg_forest.pdf'),width=10,height=6);on.exit(dev.off())
          if(!nrow(s)) {plot.new();text(.5,.5,'No estimable pancreas-fat pairs; review logs.')} else {
            par(mar=c(5,14,4,2));y <- rev(seq_len(nrow(s)))
            lab <- labels[match(s$other_trait,ids)]
            low <- s$step21_h2_z1<4|s$step21_h2_z2<4
            lab[low] <- paste0(lab[low],' [low h2 precision]')
            xr <- range(c(-1,0,1,s$ci_lower,s$ci_upper));pad <- max(diff(xr)*.04,.05)
            plot(s$rg,y,xlim=xr+c(-pad,pad),ylim=c(.5,nrow(s)+.5),pch=19,yaxt='n',ylab='',
              xlab='Genetic correlation (95% Wald CI)',main='Pancreas fat versus other traits')
            segments(s$ci_lower,y,s$ci_upper,y);abline(v=0,lty=2,col='grey50')
            axis(2,at=y,labels=lab,las=1,cex.axis=.75)
            mtext('All seven displayed pairs retain the same 28-test FDR; see accompanying table.',side=3,line=.2,cex=.7)
          }
        })
      },error=function(e) writeLines(conditionMessage(e),file.path(reports,'plot_error.txt')))
      cat('\n===== 胰腺脂肪相关比较：统一28项FDR =====\n')
      print(sub[,c('other_trait','rg','rg_se','p','FDR_BH_28','estimate_status'),drop=FALSE],row.names=FALSE)
      writeLines(if(complete) 'ALL_28_ESTIMABLE_REVIEW_QC_AND_FDR' else 'REVIEW_INCOMPLETE_OR_NONINTERPRETABLE_PAIRS',
                 file.path(reports,'analysis_review_status.txt'))
    }
  },error=function(e) {post_error <<- conditionMessage(e);writeLines(post_error,file.path(reports,'R_postprocess_error.txt'))})
  pack <- paste('import sys,zipfile,pathlib;',
    'o=pathlib.Path(sys.argv[1]);z=zipfile.ZipFile(o/"22_rg_review.zip","w",zipfile.ZIP_DEFLATED);',
    '[z.write(p,str(p.relative_to(o))) for folder in ("reports","ldsc_logs") for p in sorted((o/folder).glob("*")) if p.is_file()];',
    'z.write(o/"run_rg.py","run_rg.py");z.close()')
  pack_result <- tryCatch(system2(cfg$python,args=vapply(c('-c',pack,out),shQuote,character(1)),stdout=TRUE,stderr=TRUE),
                          error=function(e) conditionMessage(e))
  cat('\n请上传：',file.path(out,'22_rg_review.zip'),'\n')
  if(!file.exists(file.path(out,'22_rg_review.zip'))) cat('打包未成功，请压缩本次结果文件夹上传。\n')
  if(rc!=0L||!is.null(post_error)) stop('部分计算或整理未完成，请上传核查包。')
  check_status(file.path(reports,'completion_status.txt'),'COMPLETED_28_PAIR_FITS_REVIEW_ESTIMATES')
  cat('\nCOMPLETED_28_PAIR_FITS_REVIEW_ESTIMATES\n请先核查精度和异常标记，再解释FDR结果。\n')
})()
