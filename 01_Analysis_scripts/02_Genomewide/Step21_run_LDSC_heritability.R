# 第21步：8种连续表型的单表型LDSC遗传力和精度评估。
# 固定第19步环境、第20步输入；不修改数据/源码；不计算rg。
(function() {
  base <- path.expand('~/胰腺')
  prior <- file.path(base,'19_LDSC_setup','LDSC_setup_1566e430f1187')
  inputs <- file.path(base,'20_LDSC_inputs','GWAS_inputs_16c90c74e8c0e1')
  check_status <- function(path,expected) {
    if(!file.exists(path)||!identical(readLines(path,warn=FALSE),expected))
      stop('完成标记不匹配或文件缺失：',path)
  }
  check_status(file.path(prior,'completion_status.txt'),'COMPLETED_SOFTWARE_AND_REFERENCES_READY')
  check_status(file.path(inputs,'reports','completion_status.txt'),'COMPLETED_INPUT_PREPARATION_REVIEW_QC_BEFORE_LDSC')
  if(!requireNamespace('data.table',quietly=TRUE)) stop('缺少data.table。')
  cfg <- readRDS(file.path(prior,'LDSC_paths.rds'))
  if(!file.exists(cfg$python)) stop('第19步Python路径不存在。')
  parent <- file.path(base,'21_LDSC_heritability')
  dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile('h2_',tmpdir=parent)
  if(!dir.create(out)) stop('无法建立输出目录。')
  reports <- file.path(out,'reports');dir.create(reports)
  writeLines('INCOMPLETE_CHECK_LOGS',file.path(reports,'completion_status.txt'))
  writeLines(capture.output(sessionInfo()),file.path(reports,'R_sessionInfo.txt'))
  writeLines(c(
    'Step21: unpartitioned single-trait LDSC for all 8 quantitative traits; no rg estimated.',
    'Fixed CBIIT LDSC commit: 6c673952cee74bd5c57aef1555a03b1c015399a0.',
    'Calls the original ldsc.sumstats.estimate_h2 function via original argument parser; does not modify source.',
    'Free intercept; two-step cutoff 30; 200 block-jackknife blocks; no extra chisq-max filter.',
    'Original tutorial EUR LD scores used for both reference and weights; default M_5_50 denominator.',
    'Inputs already use HapMap3, GWAS MAF>=0.01, nonpalindromic SNVs, duplicate removal and MHC exclusion.',
    'Observed-scale SNP h2 of the analysed quantitative phenotype, NOT total/family heritability or a causal estimate.',
    'Full-precision h2, SE, intercept and summary statistics are extracted from Hsq, not rounded text logs.',
    '95% intervals are untruncated normal/Wald intervals; estimates outside [0,1] are retained and flagged.',
    'h2/SE<4 is a pragmatic precision-review flag, not an automatic exclusion rule or a biological conclusion.',
    'Do not select traits on favourable estimates. Review all 8 outputs before planning rg.',
    'Study-level N is assumed at each SNP; per-SNP N and imputation quality remain unavailable.',
    'GWAS phenotype transformation/covariates and sample overlap still matter for interpretation.',
    'The intercept and ratio should be interpreted with uncertainty; they alone do not establish absence/presence of confounding.',
    'Official tutorial: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation'),
    file.path(reports,'methods_and_interpretation.txt'))
  script <- file.path(base,'21_run_LDSC_heritability.R')
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
"def summarize(fit):",
"    r={k:scalar(getattr(fit,a)) for k,a in [",
"        ('h2','tot'),('h2_se','tot_se'),('intercept','intercept'),('intercept_se','intercept_se'),",
"        ('mean_chi2','mean_chisq'),('lambda_GC','lambda_gc'),('ratio','ratio'),('ratio_se','ratio_se')]}",
"    flags=[]",
"    if r['h2'] is None or r['h2_se'] is None or r['h2_se']<=0:",
"        flags.append('INVALID_H2_OR_SE');r.update(h2_z=None,h2_ci_lower=None,h2_ci_upper=None)",
"    else:",
"        r.update(h2_z=r['h2']/r['h2_se'],h2_ci_lower=r['h2']-1.959963984540054*r['h2_se'],",
"                 h2_ci_upper=r['h2']+1.959963984540054*r['h2_se'])",
"        if r['h2']<=0:flags.append('NONPOSITIVE_H2')",
"        if r['h2_z']<4:flags.append('H2_Z_LT_4_REVIEW_PRECISION')",
"        if r['h2']>1:flags.append('H2_GT_1_REVIEW_MODEL')",
"    if r['intercept'] is None or r['intercept_se'] is None or r['intercept_se']<0:",
"        flags.append('INVALID_INTERCEPT_OR_SE')",
"    else:",
"        r.update(intercept_ci_lower=r['intercept']-1.959963984540054*r['intercept_se'],",
"                 intercept_ci_upper=r['intercept']+1.959963984540054*r['intercept_se'])",
"    r['precision_flag']=';'.join(flags) if flags else 'H2_Z_GE_4_REVIEW_OTHER_DIAGNOSTICS'",
"    r['ratio_note']='UNDEFINED_MEAN_CHI2_LE_1' if r['ratio'] is None else (",
"        'NEGATIVE_RATIO_RETAINED_NOT_TRUNCATED' if r['ratio']<0 else 'INTERPRET_WITH_MEAN_CHI2_AND_SE')",
"    return r",
"",
"class Log:",
"    def __init__(self,p):self.f=open(p,'w',encoding='utf-8')",
"    def log(self,s):self.f.write(str(s)+'\\n');self.f.flush();print(s,flush=True)",
"    def close(self):self.f.close()",
"",
"def write_table(rows,path):pd.DataFrame(rows).to_csv(path,sep='\\t',index=False,float_format='%.17g',na_rep='NA')",
"",
"def main():",
"    prior=Path(sys.argv[1]);inputs=Path(sys.argv[2]);out=Path(sys.argv[3])",
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
"    results=[];commands=[]",
"    for sid,(trait,n) in EXPECTED.items():",
"        item=manifest.loc[manifest.id==sid].iloc[0];prefix=logs/(sid+'_h2');logpath=Path(str(prefix)+'.log')",
"        argv=['--h2',str(inputs/'sumstats'/(sid+'.sumstats.gz')),'--ref-ld-chr',cfg['reference_ld_chr'],",
"              '--w-ld-chr',cfg['weight_ld_chr'],'--out',str(prefix),'--n-blocks','200','--two-step','30','--print-delete-vals']",
"        args=ldsc.parser.parse_args(argv)",
"        if args.no_intercept or args.intercept_h2 is not None or args.chisq_max is not None:raise ValueError('Unexpected estimator setting')",
"        commands.append(dict(id=sid,argv=[cfg['python'],cfg['ldsc']]+argv,entrypoint='ldsc.sumstats.estimate_h2(args, log)'))",
"        (report/'commands.json').write_text(json.dumps(commands,indent=2,ensure_ascii=False))",
"        (report/(sid+'_effective_arguments.json')).write_text(json.dumps(vars(args),indent=2,ensure_ascii=False))",
"        row=dict(id=sid,trait=trait,N=n,n_input=int(item.n_final),run_status='FAILED',precision_flag='NOT_EVALUATED')",
"        log=Log(logpath);start=time.time()",
"        try:",
"            log.log(ldsc.MASTHEAD);log.log('Single-trait h2; fixed original estimator; full-precision export.')",
"            log.log('Arguments: '+json.dumps(argv,ensure_ascii=False))",
"            # Same function called by ldsc.py --h2; source code remains unchanged.",
"            fit=ldsc.sumstats.estimate_h2(args,log)",
"            row.update(summarize(fit))",
"            np.savetxt(str(prefix)+'.intercept_delete',np.asarray(fit.intercept_delete_values),fmt='%.18e')",
"            text=logpath.read_text()",
"            match=re.search(r'After merging with regression SNP LD,\\s*([0-9,]+) SNPs remain',text)",
"            if not match:raise ValueError('Cannot verify regression SNP count from log')",
"            row['n_regression']=int(match[1].replace(',',''))",
"            row['n_jackknife_blocks']=int(np.asarray(fit.tot_delete_values).shape[0])",
"            if row['n_regression']!=row['n_input']:row['precision_flag']+=';SNP_COUNT_CHANGED_REVIEW'",
"            if row['n_regression']<200000:row['precision_flag']+=';LT_200K_SNPS'",
"            row['run_status']='COMPLETED'",
"            row['elapsed_seconds']=time.time()-start",
"            log.log('Full-precision h2: '+str(row['h2'])+'; SE: '+str(row['h2_se']))",
"            log.log('Fit completed. Precision flag: '+row['precision_flag'])",
"            del fit;gc.collect()",
"        except Exception:",
"            row['error']=traceback.format_exc();log.log(row['error'])",
"        finally:log.close()",
"        (report/(sid+'_full_precision.json')).write_text(json.dumps(row,indent=2,ensure_ascii=False,allow_nan=False))",
"        results.append(row);write_table(results,report/'LDSC_h2_summary.tsv')",
"        print('FINISHED',sid,row['run_status'],row['precision_flag'],flush=True)",
"    warnings=[]",
"    for p in sorted(logs.glob('*.log')):",
"        for line in p.read_text().splitlines():",
"            if any(k in line.lower() for k in ['warning','traceback','error','ratio <','ratio na']):warnings.append(dict(file=p.name,message=line))",
"    pd.DataFrame(warnings,columns=['file','message']).to_csv(report/'LDSC_log_messages.tsv',sep='\\t',index=False)",
"    success=all(r['run_status']=='COMPLETED' for r in results)",
"    status.write_text('COMPLETED_8_TRAIT_H2_REVIEW_ESTIMATES\\n' if success else 'PARTIAL_FAILURE_CHECK_LOGS\\n')",
"    if not success:raise RuntimeError('One or more traits failed; other results retained. Upload review ZIP.')",
"    print('COMPLETED_8_TRAIT_H2_REVIEW_ESTIMATES',flush=True)",
"",
"if __name__=='__main__':",
"    try:main()",
"    except Exception:",
"        out=Path(sys.argv[3]);(out/'reports').mkdir(exist_ok=True)",
"        (out/'reports'/'python_error.txt').write_text(traceback.format_exc())",
"        raise"
)
  helper <- file.path(out,'run_h2.py');writeLines(python_lines,helper,useBytes=TRUE)
  log <- file.path(reports,'python_console.txt')
  cat('结果目录：',out,'\n')
  cat('将顺序计算8个表型。详细进度：',log,'\n');flush.console()
  err <- NULL
  rc <- tryCatch(system2(cfg$python,
    args=vapply(c('-u',helper,prior,inputs,out),shQuote,character(1)),stdout=log,stderr=log,
    env=c('PYTHONNOUSERSITE=1','PYTHONPATH=','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1')),
    error=function(e) {err <<- conditionMessage(e);127L})
  if(!is.null(err)) writeLines(err,file.path(reports,'R_error.txt'))
  # 图仅展示实际完成且数值有效的结果，不将失败或负值改成零。
  summary_path <- file.path(reports,'LDSC_h2_summary.tsv')
  if(file.exists(summary_path)) {
    tab <- data.table::fread(summary_path,na.strings='NA')
    if(all(c('h2','h2_se','h2_ci_lower','h2_ci_upper') %in% names(tab))) {
      d <- as.data.frame(tab[run_status=='COMPLETED'])
      d <- d[is.finite(d$h2)&is.finite(d$h2_se)&d$h2_se>0,,drop=FALSE]
      if(nrow(d)>0L) {
        tryCatch({
          local({
            pdf(file.path(reports,'LDSC_h2_forest.pdf'),width=10,height=6)
            on.exit(dev.off())
            par(mar=c(5,17,3,2))
            y <- rev(seq_len(nrow(d)))
            limits <- range(c(0,d$h2_ci_lower,d$h2_ci_upper))
            pad <- max(diff(limits)*0.08,0.01)
            plot(d$h2,y,xlim=limits+c(-pad,pad),ylim=c(0.5,nrow(d)+0.5),
                 yaxt='n',ylab='',xlab='Observed-scale SNP heritability (95% Wald CI)',
                 pch=19,main='Single-trait LDSC: eight-trait assessment')
            axis(2,at=y,labels=gsub('_',' ',d$trait),las=1,cex.axis=0.75)
            segments(d$h2_ci_lower,y,d$h2_ci_upper,y)
            abline(v=0,lty=2,col='grey50')
          })
        },error=function(e) writeLines(conditionMessage(e),file.path(reports,'plot_error.txt')))
      }
    }
    cat('\n===== 第21步结果汇总 =====\n')
    cols <- intersect(c('trait','N','h2','h2_se','h2_z','intercept','intercept_se','run_status','precision_flag'),names(tab))
    print(tab[,cols,with=FALSE])
  }
  # 打包所有小型报告、原始日志、jackknife数值和执行脚本；不复制大型GWAS/参考文件。
  pack <- paste('import sys,zipfile,pathlib;',
    'o=pathlib.Path(sys.argv[1]);z=zipfile.ZipFile(o/"21_h2_review.zip","w",zipfile.ZIP_DEFLATED);',
    '[z.write(p,str(p.relative_to(o))) for folder in ("reports","ldsc_logs") for p in sorted((o/folder).glob("*")) if p.is_file()];',
    'z.write(o/"run_h2.py","run_h2.py");z.close()')
  pack_result <- tryCatch(system2(cfg$python,args=vapply(c('-c',pack,out),shQuote,character(1)),
                                  stdout=TRUE,stderr=TRUE),error=function(e) conditionMessage(e))
  if(!file.exists(file.path(out,'21_h2_review.zip'))) cat('自动打包未成功，请直接压缩本次结果文件夹上传。\n')
  cat('\n请上传：',file.path(out,'21_h2_review.zip'),'\n')
  if(rc!=0L) {
    if(file.exists(log)) cat(tail(readLines(log,warn=FALSE),20L),sep='\n')
    stop('有步骤或表型未完成，请上传核查包，不要直接进入rg分析。')
  }
  check_status(file.path(reports,'completion_status.txt'),'COMPLETED_8_TRAIT_H2_REVIEW_ESTIMATES')
  cat('\nCOMPLETED_8_TRAIT_H2_REVIEW_ESTIMATES\n')
  cat('计算完成不等于遗传力显著。请先核查估计值、标准误和精度标记。\n')
})()
