# 第23步：排除整条chr8后的8表型28组LDSC遗传相关敏感性分析。
# 目的：检验第22步全基因组rg是否依赖chr8强共享区域。
# 原则：不修改第20步原始LDSC输入；不改变LDSC版本、自由截距、200 jackknife区块、two-step=30或28项BH-FDR。
# 仅做两项变化：
#   1) 从8份LDSC sumstats中删除chr8 SNP；
#   2) 构建对应的“无chr8”LD-score参考副本：chr1-7/9-22复用原参考，chr8置空且M_5_50=0。
# 输出：28组no-chr8 rg、统一28项FDR、与第22步full-genome结果的描述性比较。
# 注意：full vs no-chr8高度相关，因此不做把两者当独立估计的正式差异Z检验。

(function() {
  base <- path.expand('~/胰腺')
  prior <- file.path(base,'19_LDSC_setup','LDSC_setup_1566e430f1187')
  inputs <- file.path(base,'20_LDSC_inputs','GWAS_inputs_16c90c74e8c0e1')
  h2dir <- file.path(base,'21_LDSC_heritability','h2_16c90c661ce88d')
  step22_parent <- file.path(base,'22_LDSC_genetic_correlation')

  check_status <- function(path,expected) {
    if(!file.exists(path) || !identical(readLines(path,warn=FALSE),expected))
      stop('完成标记不一致或文件缺失：',path)
  }
  check_status(file.path(prior,'completion_status.txt'),'COMPLETED_SOFTWARE_AND_REFERENCES_READY')
  check_status(file.path(inputs,'reports','completion_status.txt'),'COMPLETED_INPUT_PREPARATION_REVIEW_QC_BEFORE_LDSC')
  check_status(file.path(h2dir,'reports','completion_status.txt'),'COMPLETED_8_TRAIT_H2_REVIEW_ESTIMATES')
  if(!requireNamespace('data.table',quietly=TRUE)) stop('缺少data.table。')

  cfg <- readRDS(file.path(prior,'LDSC_paths.rds'))
  if(!identical(cfg$source_commit,'6c673952cee74bd5c57aef1555a03b1c015399a0'))
    stop('第19步LDSC源码版本不一致。')
  if(!file.exists(cfg$python) || !file.exists(cfg$ldsc)) stop('第19步Python/LDSC路径不存在。')

  # 自动寻找唯一一个已经完成且28组均可解释的第22步目录；若有多个则停止，避免误用旧批次。
  find_step22 <- function(parent) {
    if(!dir.exists(parent)) stop('第22步父目录不存在：',parent)
    ds <- list.dirs(parent,recursive=FALSE,full.names=TRUE)
    ok <- vapply(ds,function(d) {
      s1 <- file.path(d,'reports','completion_status.txt')
      s2 <- file.path(d,'reports','analysis_review_status.txt')
      f  <- file.path(d,'reports','LDSC_rg_summary_FDR.tsv')
      file.exists(s1) && identical(readLines(s1,warn=FALSE),'COMPLETED_28_PAIR_FITS_REVIEW_ESTIMATES') &&
        file.exists(s2) && identical(readLines(s2,warn=FALSE),'ALL_28_ESTIMABLE_REVIEW_QC_AND_FDR') &&
        file.exists(f)
    },logical(1))
    ds <- ds[ok]
    if(length(ds)!=1L) {
      stop('应有且仅有1个通过核查的第22步目录；当前找到',length(ds),'个：\n',paste(ds,collapse='\n'),
           '\n请先把旧批次移出该父目录，或在脚本顶部改为手工指定source22。')
    }
    ds[[1]]
  }
  source22 <- find_step22(step22_parent)
  full_file <- file.path(source22,'reports','LDSC_rg_summary_FDR.tsv')
  full <- as.data.frame(data.table::fread(full_file,na.strings='NA'))
  if(nrow(full)!=28L || any(full$run_status!='COMPLETED') || any(full$estimate_status!='ESTIMABLE'))
    stop('第22步汇总表不是28组全部完成且可解释。')

  parent <- file.path(base,'23_LDSC_rg_exclude_chr8')
  dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  out <- tempfile('rg_nochr8_',tmpdir=parent)
  if(!dir.create(out)) stop('无法创建第23步结果目录。')
  reports <- file.path(out,'reports');dir.create(reports)
  writeLines('INCOMPLETE_CHECK_LOGS',file.path(reports,'completion_status.txt'))
  writeLines(capture.output(sessionInfo()),file.path(reports,'R_sessionInfo.txt'))
  data.table::fwrite(data.frame(source22=source22,full_summary=full_file),
                     file.path(reports,'step22_source.tsv'),sep='\t')
  writeLines(c(
    'Step23 is a chromosome-8 exclusion sensitivity analysis of Step22 genome-wide LDSC rg.',
    'All 28 unordered pairs are re-estimated; no pair is selected on significance.',
    'Only chromosome 8 is removed. Step20 allele/QC decisions and all other chromosomes remain unchanged.',
    'The fixed CBIIT LDSC commit, free h2/gencov intercepts, allele checking, 200 jackknife blocks and two-step cutoff 30 are retained.',
    'A derived no-chr8 LD-score reference is created: chromosomes 1-7 and 9-22 point to the validated original files; chromosome 8 contains no LD-score rows and M_5_50=0.',
    'This prevents chromosome 8 from contributing either regression SNPs or the reference M denominator.',
    'BH-FDR is recalculated across the same fixed family of 28 no-chr8 tests.',
    'Full-genome versus no-chr8 changes are descriptive sensitivity comparisons only.',
    'Do NOT compute a naive independent-estimates Z test for the difference because the two rg estimates share nearly all SNPs.',
    'Loss of FDR significance alone is not treated as instability; direction, magnitude and uncertainty are reviewed together.',
    'No new chr8 regional exploration is performed.'),
    file.path(reports,'methods_and_interpretation.txt'))
  script <- file.path(base,'23_exclude_chr8_LDSC_rg_sensitivity.R')
  if(file.exists(script)) file.copy(script,file.path(reports,basename(script)),overwrite=FALSE)

  python_lines <- c(
"\"\"\"Step23: exclude chromosome 8 and re-run all 28 LDSC genetic correlations.\"\"\"",
"import sys, json, re, hashlib, traceback, time, importlib.util, gc, platform, gzip, os",
"from pathlib import Path",
"from itertools import combinations",
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
"class Log:",
"    def __init__(self,p):self.f=open(p,'w',encoding='utf-8')",
"    def log(self,s):self.f.write(str(s)+'\\n');self.f.flush();print(s,flush=True)",
"    def close(self):self.f.close()",
"",
"def write_table(rows,path):",
"    pd.DataFrame(rows).to_csv(path,sep='\\t',index=False,float_format='%.17g',na_rep='NA')",
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
"            if row[k+'_z'] is None or row[k+'_z']<4:flags.append(k+'_LOW_PRECISION_NOCHR8')",
"    valid=all(row[k] is not None for k in ['rg_internal','rg_se_internal','z_internal','p_internal'])",
"    valid=valid and row['rg_se_internal']>0 and 0<=row['p_internal']<=1",
"    if valid and not np.isclose(row['rg_internal']/row['rg_se_internal'],row['z_internal'],rtol=1e-8,atol=1e-10):",
"        valid=False;flags.append('Z_INCONSISTENCY')",
"    if getattr(fit,'_negative_hsq',False):valid=False;flags.append('NONPOSITIVE_PAIR_H2')",
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
"def build_nochr8_reference(cfg,out,report):",
"    src_prefix=cfg['reference_ld_chr']",
"    dst=out/'ld_reference_nochr8';dst.mkdir(exist_ok=True)",
"    rows=[]",
"    for ch in range(1,23):",
"        src_ld=Path(src_prefix+str(ch)+'.l2.ldscore.gz')",
"        src_m5=Path(src_prefix+str(ch)+'.l2.M_5_50')",
"        if not src_ld.is_file() or not src_m5.is_file():raise ValueError('Missing reference file for chr'+str(ch))",
"        dst_ld=dst/(str(ch)+'.l2.ldscore.gz');dst_m5=dst/(str(ch)+'.l2.M_5_50')",
"        src_m=Path(src_prefix+str(ch)+'.l2.M');dst_m=dst/(str(ch)+'.l2.M')",
"        if ch!=8:",
"            for src,dest in [(src_ld,dst_ld),(src_m5,dst_m5)]+([(src_m,dst_m)] if src_m.is_file() else []):",
"                if dest.exists() or dest.is_symlink():dest.unlink()",
"                dest.symlink_to(src.resolve())",
"            rows.append(dict(CHR=ch,mode='SYMLINK_ORIGINAL',n_ld_rows=None,M_5_50=float(np.loadtxt(src_m5))))",
"        else:",
"            hdr=pd.read_csv(src_ld,sep=r'\\s+',nrows=0)",
"            hdr.to_csv(dst_ld,sep='\\t',index=False,compression='gzip')",
"            dst_m5.write_text('0\\n')",
"            if src_m.is_file():dst_m.write_text('0\\n')",
"            rows.append(dict(CHR=ch,mode='EMPTY_EXCLUDED',n_ld_rows=0,M_5_50=0.0))",
"    write_table(rows,report/'nochr8_reference_manifest.tsv')",
"    check=pd.read_csv(dst/'8.l2.ldscore.gz',sep=r'\\s+')",
"    if len(check)!=0:raise ValueError('Derived chr8 LD score is not empty')",
"    if float(np.loadtxt(dst/'8.l2.M_5_50'))!=0:raise ValueError('Derived chr8 M_5_50 is not zero')",
"    return str(dst)+'/'",
"",
"def build_nochr8_sumstats(inputs,out,manifest,report,cfg):",
"    dest=out/'sumstats_nochr8';dest.mkdir(exist_ok=True)",
"    chr8_ref=pd.read_csv(cfg['reference_ld_chr']+'8.l2.ldscore.gz',sep=r'\\s+',usecols=['SNP'])",
"    if not chr8_ref.SNP.is_unique:raise ValueError('Reference chr8 SNP IDs are duplicated')",
"    chr8=set(chr8_ref.SNP.astype(str))",
"    sets={};rows=[]",
"    for rr in manifest.to_dict('records'):",
"        sid=rr['id'];trait,n=EXPECTED[sid]",
"        src=inputs/'sumstats'/(sid+'.sumstats.gz')",
"        if not src.is_file() or digest(src)!=rr['sumstats_sha256']:raise ValueError('Step20 sumstats hash mismatch: '+sid)",
"        x=pd.read_csv(src,sep=r'\\s+',dtype={'SNP':str,'A1':str,'A2':str})",
"        if list(x.columns)!=['SNP','A1','A2','Z','N']:raise ValueError('Unexpected step20 columns: '+sid)",
"        if len(x)!=int(rr['n_final']) or not x.SNP.is_unique:raise ValueError('Unexpected step20 row count/duplicates: '+sid)",
"        cfile=inputs/'reports'/(sid+'_chromosome_counts.tsv')",
"        cc=pd.read_csv(cfile,sep='\\t')",
"        hit=cc.loc[cc.CHR.astype(int)==8,'n']",
"        if len(hit)!=1:raise ValueError('Cannot recover step20 chr8 count: '+sid)",
"        expected_removed=int(hit.iloc[0])",
"        mask=x.SNP.isin(chr8)",
"        observed_removed=int(mask.sum())",
"        if observed_removed!=expected_removed:",
"            raise ValueError('chr8 removal count mismatch for %s: expected %d observed %d'%(sid,expected_removed,observed_removed))",
"        y=x.loc[~mask].copy()",
"        if len(y)<200000 or y.SNP.isin(chr8).any():raise ValueError('Invalid no-chr8 sumstats: '+sid)",
"        if not (y.N==n).all() or not np.isfinite(y.Z).all():raise ValueError('N/Z invalid after chr8 removal: '+sid)",
"        target=dest/(sid+'.sumstats.gz')",
"        y.to_csv(target,sep='\\t',index=False,compression='gzip',float_format='%.12g')",
"        # exact readback of retained values",
"        z=pd.read_csv(target,sep=r'\\s+',dtype={'SNP':str,'A1':str,'A2':str})",
"        if len(z)!=len(y) or not np.array_equal(z.SNP.to_numpy(),y.SNP.to_numpy()):raise ValueError('No-chr8 readback SNP mismatch: '+sid)",
"        if not np.array_equal(z.A1.to_numpy(),y.A1.to_numpy()) or not np.array_equal(z.A2.to_numpy(),y.A2.to_numpy()):raise ValueError('No-chr8 allele mismatch: '+sid)",
"        if not np.allclose(z.Z.to_numpy(),y.Z.to_numpy(),rtol=1e-10,atol=1e-12):raise ValueError('No-chr8 Z mismatch: '+sid)",
"        sets[sid]=set(y.SNP)",
"        rows.append(dict(id=sid,trait=trait,N=n,n_full=len(x),n_chr8_removed=observed_removed,n_nochr8=len(y),",
"                         source_sumstats=str(src),source_sha256=digest(src),nochr8_sumstats=str(target),nochr8_sha256=digest(target)))",
"    write_table(rows,report/'nochr8_input_manifest.tsv')",
"    overlaps=[]",
"    ids=list(EXPECTED.keys())",
"    for i,a in enumerate(ids):",
"        for b in ids[i+1:]:",
"            overlaps.append(dict(id1=a,id2=b,trait1=EXPECTED[a][0],trait2=EXPECTED[b][0],n_common=len(sets[a]&sets[b])))",
"    write_table(overlaps,report/'pairwise_SNP_overlap_nochr8.tsv')",
"    return dest,sets,{frozenset([r['trait1'],r['trait2']]):int(r['n_common']) for r in overlaps}",
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
"    for filename,col,kind in [('reference_file_checksums.tsv','sha256','sha256'),('source_code_checksums.tsv','md5','md5')]:",
"        for row in pd.read_csv(prior/filename,sep='\\t').to_dict('records'):",
"            p=Path(row['path']);ok=p.is_file() and digest(p,kind)==row[col]",
"            checks.append(dict(kind=kind,path=str(p),match=ok))",
"            if not ok:raise ValueError('Resource checksum mismatch: '+str(p))",
"    write_table(checks,report/'step19_resource_recheck.tsv')",
"    for rr in manifest.to_dict('records'):",
"        sid=rr['id'];trait,n=EXPECTED[sid]",
"        if rr['trait']!=trait or int(rr['N'])!=n or rr['status']!='PASS_INPUT_QC':raise ValueError('Step20 QC/trait mismatch: '+sid)",
"    manifest.to_csv(report/'step20_input_manifest.tsv',sep='\\t',index=False)",
"    # Build derived no-chr8 inputs and no-chr8 reference without touching prior steps.",
"    no_ref=build_nochr8_reference(cfg,out,report)",
"    no_sumstats,sets,expected_overlap=build_nochr8_sumstats(inputs,out,manifest,report,cfg)",
"    if len(expected_overlap)!=28:raise ValueError('Expected 28 no-chr8 overlap comparisons')",
"    sys.path.insert(0,str(Path(cfg['ldsc']).parent))",
"    spec=importlib.util.spec_from_file_location('fixed_ldsc23',cfg['ldsc'])",
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
"    h2.to_csv(report/'step21_fullgenome_h2_context.tsv',sep='\\t')",
"    results=[];commands=[];pairs=list(combinations(EXPECTED,2))",
"    write_table([dict(pair_id='pair_%02d'%(i+1),id1=a,id2=b,trait1=EXPECTED[a][0],trait2=EXPECTED[b][0])",
"                 for i,(a,b) in enumerate(pairs)],report/'planned_28_pairs.tsv')",
"    for i,(a,b) in enumerate(pairs,1):",
"        pair_id='pair_%02d'%i;ta,na=EXPECTED[a];tb,nb=EXPECTED[b]",
"        prefix=logs/(pair_id+'_');logpath=logs/(pair_id+'.log')",
"        argv=['--rg',str(no_sumstats/(a+'.sumstats.gz'))+','+str(no_sumstats/(b+'.sumstats.gz')),",
"              '--ref-ld-chr',no_ref,'--w-ld-chr',no_ref,'--out',str(prefix),",
"              '--n-blocks','200','--two-step','30','--print-delete-vals']",
"        args=ldsc.parser.parse_args(argv)",
"        if args.no_intercept or args.intercept_h2 is not None or args.intercept_gencov is not None or args.no_check_alleles or args.return_silly_things:",
"            raise ValueError('Unexpected intercept/allele/rg setting')",
"        commands.append(dict(pair_id=pair_id,argv=[cfg['python'],cfg['ldsc']]+argv,entrypoint='ldsc.sumstats.estimate_rg(args, log)'))",
"        (report/'commands.json').write_text(json.dumps(commands,indent=2,ensure_ascii=False))",
"        row=dict(pair_id=pair_id,id1=a,id2=b,trait1=ta,trait2=tb,N1=na,N2=nb,",
"                 step21_full_h2_z1=float(h2.loc[a,'h2_z']),step21_full_h2_z2=float(h2.loc[b,'h2_z']),",
"                 involves_liver_iron='liver_iron' in (ta,tb),involves_pancreas_fat='pancreas_fat' in (ta,tb),",
"                 run_status='FAILED',estimate_status='NOT_ESTIMATED',rg=None,rg_se=None,p=None,ci_lower=None,ci_upper=None,pair_QC_flag='NOT_EVALUATED')",
"        log=Log(logpath);start=time.time()",
"        try:",
"            log.log(ldsc.MASTHEAD);log.log('Step23 no-chr8 pair %d/28: %s vs %s'%(i,ta,tb));log.log(json.dumps(argv,ensure_ascii=False))",
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
"        results.append(row);write_table(results,report/'LDSC_rg_nochr8_raw_summary.tsv')",
"        print('FINISHED',i,'/28',ta,tb,row['run_status'],flush=True)",
"    messages=[]",
"    for p in sorted(logs.glob('*.log')):",
"        for line in p.read_text().splitlines():",
"            if any(k in line.lower() for k in ['warning','traceback','error','ratio <']):messages.append(dict(file=p.name,message=line))",
"    pd.DataFrame(messages,columns=['file','message']).to_csv(report/'LDSC_log_messages.tsv',sep='\\t',index=False)",
"    success=all(r['run_status']=='COMPLETED' for r in results)",
"    status.write_text('COMPLETED_NOCHR8_28_PAIR_FITS_REVIEW_ESTIMATES\\n' if success else 'PARTIAL_FAILURE_CHECK_LOGS\\n')",
"    if not success:raise RuntimeError('Some no-chr8 pairs failed; upload review package')",
"",
"if __name__=='__main__':",
"    try:main()",
"    except Exception:",
"        out=Path(sys.argv[3]);(out/'reports').mkdir(exist_ok=True)",
"        (out/'reports'/'python_error.txt').write_text(traceback.format_exc())",
"        raise"
  )

  helper <- file.path(out,'run_rg_nochr8.py')
  writeLines(python_lines,helper,useBytes=TRUE)
  log <- file.path(reports,'python_console.txt')
  cat('结果目录：',out,'\n')
  cat('第23步将先生成无chr8输入和参考，再顺序计算28组rg。详细日志：',log,'\n')
  flush.console()

  err <- NULL
  rc <- tryCatch(system2(cfg$python,
    args=vapply(c('-u',helper,prior,inputs,out,h2dir),shQuote,character(1)),
    stdout=log,stderr=log,
    env=c('PYTHONNOUSERSITE=1','PYTHONPATH=','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1')),
    error=function(e) {err <<- conditionMessage(e);127L})
  if(!is.null(err)) writeLines(err,file.path(reports,'R_error.txt'))

  post_error <- NULL
  tryCatch({
    raw <- file.path(reports,'LDSC_rg_nochr8_raw_summary.tsv')
    if(!file.exists(raw)) stop('未生成no-chr8 rg原始汇总表。')
    d <- as.data.frame(data.table::fread(raw,na.strings='NA'))
    for(k in c('p','rg','rg_se','ci_lower','ci_upper')) d[[k]] <- as.numeric(d[[k]])
    valid <- d$run_status=='COMPLETED' & d$estimate_status=='ESTIMABLE' & is.finite(d$p) & d$p>=0 & d$p<=1
    d$FDR_BH_28 <- NA_real_
    if(any(valid)) d$FDR_BH_28[valid] <- p.adjust(d$p[valid],method='BH',n=28L)
    d$FDR_significant <- ifelse(is.na(d$FDR_BH_28),NA,d$FDR_BH_28<0.05)
    complete <- nrow(d)==28L && all(valid)
    d$FDR_scope_status <- if(complete) 'ALL_28_ESTIMABLE_NOCHR8' else 'PROVISIONAL_FIXED_28_REVIEW_MISSING'
    data.table::fwrite(d,file.path(reports,'LDSC_rg_nochr8_summary_FDR.tsv'),sep='\t',na='NA')

    # 与第22步进行描述性比较；不进行错误的“独立估计差异Z检验”。
    full2 <- full
    keyfun <- function(a,b) vapply(seq_along(a),function(i) paste(sort(c(a[i],b[i])),collapse='__'),character(1))
    full2$pair_key <- keyfun(full2$trait1,full2$trait2)
    d$pair_key <- keyfun(d$trait1,d$trait2)
    if(anyDuplicated(full2$pair_key) || anyDuplicated(d$pair_key)) stop('配对键重复。')
    keep_full <- c('pair_key','rg','rg_se','ci_lower','ci_upper','p','FDR_BH_28','FDR_significant','n_regression')
    keep_no <- c('pair_key','trait1','trait2','rg','rg_se','ci_lower','ci_upper','p','FDR_BH_28','FDR_significant','n_regression','pair_QC_flag')
    ff <- full2[,keep_full,drop=FALSE]
    nn <- d[,keep_no,drop=FALSE]
    names(ff)[names(ff)!='pair_key'] <- paste0(names(ff)[names(ff)!='pair_key'],'_full')
    names(nn)[!names(nn)%in%c('pair_key','trait1','trait2')] <- paste0(names(nn)[!names(nn)%in%c('pair_key','trait1','trait2')],'_nochr8')
    cmp <- merge(nn,ff,by='pair_key',all.x=TRUE,sort=FALSE)
    if(nrow(cmp)!=28L || anyNA(cmp$rg_full)) stop('与第22步结果合并失败。')
    cmp$delta_rg_nochr8_minus_full <- cmp$rg_nochr8-cmp$rg_full
    cmp$direction_same <- sign(cmp$rg_nochr8)==sign(cmp$rg_full)
    cmp$abs_rg_retained_ratio <- ifelse(abs(cmp$rg_full)>1e-12,abs(cmp$rg_nochr8)/abs(cmp$rg_full),NA_real_)
    cmp$abs_rg_change_percent <- ifelse(abs(cmp$rg_full)>1e-12,
      100*(abs(cmp$rg_nochr8)-abs(cmp$rg_full))/abs(cmp$rg_full),NA_real_)
    cmp$FDR_status <- ifelse(cmp$FDR_significant_full & cmp$FDR_significant_nochr8,'significant_both',
      ifelse(cmp$FDR_significant_full & !cmp$FDR_significant_nochr8,'full_only',
        ifelse(!cmp$FDR_significant_full & cmp$FDR_significant_nochr8,'nochr8_only','nonsignificant_both')))
    cmp$interpretation_note <- 'Descriptive sensitivity comparison only; no independent-estimates difference test.'
    data.table::fwrite(cmp,file.path(reports,'full_vs_nochr8_rg_comparison.tsv'),sep='\t',na='NA')

    # 胰腺脂肪7组 + 胰腺内部三组重点比较。
    pf <- cmp[cmp$trait1=='pancreas_fat'|cmp$trait2=='pancreas_fat',,drop=FALSE]
    pf$other_trait <- ifelse(pf$trait1=='pancreas_fat',pf$trait2,pf$trait1)
    data.table::fwrite(pf,file.path(reports,'pancreas_fat_full_vs_nochr8.tsv'),sep='\t',na='NA')
    trio_pairs <- c('pancreas_fat__pancreas_iron','pancreas_fat__pancreas_volume','pancreas_iron__pancreas_volume')
    trio <- cmp[cmp$pair_key %in% trio_pairs,,drop=FALSE]
    data.table::fwrite(trio,file.path(reports,'pancreas_three_trait_full_vs_nochr8.tsv'),sep='\t',na='NA')

    ids <- c('liver_volume','pancreas_volume','visceral_adipose_tissue_volume',
      'abdominal_subcutaneous_adipose_tissue_volume','liver_fat','liver_iron','pancreas_fat','pancreas_iron')
    labels <- c('Liver volume','Pancreas volume','VAT volume','ASAT volume','Liver fat','Liver iron','Pancreas fat','Pancreas iron')
    mat <- matrix(NA_real_,8,8,dimnames=list(ids,ids));diag(mat)<-1
    qmat <- matrix(NA_real_,8,8,dimnames=list(ids,ids))
    for(i in seq_len(nrow(d))) {
      a <- match(d$trait1[i],ids);b <- match(d$trait2[i],ids)
      mat[a,b] <- mat[b,a] <- d$rg[i];qmat[a,b] <- qmat[b,a] <- d$FDR_BH_28[i]
    }
    data.table::fwrite(data.frame(trait=ids,mat,check.names=FALSE),file.path(reports,'rg_nochr8_matrix.tsv'),sep='\t',na='NA')
    data.table::fwrite(data.frame(trait=ids,qmat,check.names=FALSE),file.path(reports,'FDR_nochr8_matrix.tsv'),sep='\t',na='NA')

    tryCatch({
      local({
        pdf(file.path(reports,'rg_nochr8_heatmap.pdf'),width=9,height=9);on.exit(dev.off())
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
        title('LDSC genetic correlations after excluding chromosome 8')
        mtext('* BH-FDR <0.05 across 28 no-chr8 tests',side=3,line=.2,cex=.7)
      })
      local({
        z <- cmp[is.finite(cmp$rg_full)&is.finite(cmp$rg_nochr8),,drop=FALSE]
        pdf(file.path(reports,'full_vs_nochr8_rg_scatter.pdf'),width=7,height=7);on.exit(dev.off())
        lim <- range(c(-1,1,z$rg_full,z$rg_nochr8));
        plot(z$rg_full,z$rg_nochr8,xlim=lim,ylim=lim,pch=19,
             xlab='Full-genome rg (Step22)',ylab='rg after excluding chr8 (Step23)',
             main='Full genome vs chromosome-8-excluded rg')
        abline(0,1,lty=2,col='grey50');abline(h=0,v=0,lty=3,col='grey75')
        ix <- z$trait1=='pancreas_fat'|z$trait2=='pancreas_fat'|
              (z$pair_key %in% trio_pairs)
        text(z$rg_full[ix],z$rg_nochr8[ix],labels=z$pair_key[ix],pos=3,cex=.55)
      })
    },error=function(e) writeLines(conditionMessage(e),file.path(reports,'plot_error.txt')))

    writeLines(if(complete) 'ALL_28_NOCHR8_ESTIMABLE_REVIEW_QC_AND_FDR' else 'REVIEW_INCOMPLETE_OR_NONINTERPRETABLE_PAIRS',
               file.path(reports,'analysis_review_status.txt'))
    cat('\n===== 第23步重点：胰腺脂肪7组 full vs no-chr8 =====\n')
    print(pf[,c('other_trait','rg_full','rg_nochr8','delta_rg_nochr8_minus_full','abs_rg_change_percent',
                 'FDR_BH_28_full','FDR_BH_28_nochr8','direction_same'),drop=FALSE],row.names=FALSE)
  },error=function(e) {
    post_error <<- conditionMessage(e)
    writeLines(post_error,file.path(reports,'R_postprocess_error.txt'))
  })

  # 只打包小型结果/日志，不打包8个no-chr8 sumstats和派生LD参考。
  pack <- paste('import sys,zipfile,pathlib;',
    'o=pathlib.Path(sys.argv[1]);z=zipfile.ZipFile(o/"23_nochr8_rg_review.zip","w",zipfile.ZIP_DEFLATED);',
    '[z.write(p,str(p.relative_to(o))) for folder in ("reports","ldsc_logs") for p in sorted((o/folder).glob("*")) if p.is_file()];',
    'z.write(o/"run_rg_nochr8.py","run_rg_nochr8.py");z.close()')
  pack_result <- tryCatch(system2(cfg$python,args=vapply(c('-c',pack,out),shQuote,character(1)),stdout=TRUE,stderr=TRUE),
                          error=function(e) conditionMessage(e))

  cat('\n请上传核查包：',file.path(out,'23_nochr8_rg_review.zip'),'\n')
  if(!file.exists(file.path(out,'23_nochr8_rg_review.zip'))) cat('打包失败，请压缩本次结果目录中的reports和ldsc_logs上传。\n')
  if(rc!=0L || !is.null(post_error)) stop('第23步部分计算或整理未完成，请上传核查包。')
  check_status(file.path(reports,'completion_status.txt'),'COMPLETED_NOCHR8_28_PAIR_FITS_REVIEW_ESTIMATES')
  cat('\nCOMPLETED_NOCHR8_28_PAIR_FITS_REVIEW_ESTIMATES\n')
  cat('请先核查no-chr8输入计数、28组QC标记及full-vs-nochr8比较表，再解释敏感性结果。\n')
})()
