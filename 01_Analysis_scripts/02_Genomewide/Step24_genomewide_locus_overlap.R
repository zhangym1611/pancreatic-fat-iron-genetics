# 第24步：胰腺脂肪与胰腺铁全基因组独立显著位点识别及locus overlap分析
# 目的：
# 1) 从两份原始GWAS VCF中重新提取全基因组显著(P<5e-8)的双等位SNV；
# 2) 使用既有1000G EUR参考和PLINK进行独立位点clumping；
# 3) 以独立lead SNP为中心±500 kb构建操作性genomic loci并合并重叠窗口；
# 4) 比较胰腺脂肪与胰腺铁的exact shared variants、lead SNP及物理locus overlap；
# 5) chr8仅作为既有重点区域标记，不追加新的chr8机制探索。
#
# 重要：
# - 本步不使用第20步的HapMap3 LDSC子集来发现GWAS位点，避免漏掉非HapMap3的显著变异；
# - clumping仅用于定义独立lead SNP，不用于因果推断；
# - 物理窗口重叠不等于共定位，更不等于共享因果变异；
# - 后续local genetic correlation / cross-trait colocalization负责检验真正的局部共享信号。

(function() {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("缺少R包：data.table")

  base <- path.expand("~/胰腺")
  prior19 <- file.path(base, "19_LDSC_setup", "LDSC_setup_1566e430f1187")
  prior20 <- file.path(base, "20_LDSC_inputs", "GWAS_inputs_16c90c74e8c0e1")
  plink_bin <- path.expand("~/mr_tools/plink")
  ref_prefix <- file.path(base, "00_ld_reference", "EUR")

  check_status <- function(path, expected) {
    if (!file.exists(path)) stop("缺少完成标记：", path)
    x <- readLines(path, warn = FALSE)
    if (!length(x) || !identical(x[1], expected)) stop("完成标记不匹配：", path)
  }

  check_status(file.path(prior19, "completion_status.txt"),
               "COMPLETED_SOFTWARE_AND_REFERENCES_READY")
  check_status(file.path(prior20, "reports", "completion_status.txt"),
               "COMPLETED_INPUT_PREPARATION_REVIEW_QC_BEFORE_LDSC")

  cfg <- readRDS(file.path(prior19, "LDSC_paths.rds"))
  if (!identical(cfg$source_commit, "6c673952cee74bd5c57aef1555a03b1c015399a0"))
    stop("第19步LDSC源码版本与既往分析不一致。")
  if (!file.exists(cfg$python)) stop("第19步Python环境不存在：", cfg$python)
  if (!file.exists(plink_bin) || file.access(plink_bin, 1) != 0L)
    stop("PLINK不存在或没有执行权限：", plink_bin)
  if (!all(file.exists(paste0(ref_prefix, c(".bed", ".bim", ".fam")))))
    stop("1000G EUR PLINK参考文件不完整：", ref_prefix)

  ids <- c(
    pancreas_fat  = "ebi-a-GCST90016675",
    pancreas_iron = "ebi-a-GCST90016676"
  )
  N <- c(pancreas_fat = 25617L, pancreas_iron = 25617L)
  vcf <- file.path(base, paste0(ids, ".vcf.gz"))
  if (!all(file.exists(vcf))) stop("缺少原始VCF：", paste(vcf[!file.exists(vcf)], collapse = ", "))

  parent <- file.path(base, "24_genomewide_locus_overlap")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile("fat_iron_loci_", tmpdir = parent)
  if (!dir.create(out)) stop("无法建立输出目录：", out)

  reports <- file.path(out, "reports")
  extracted <- file.path(out, "GWS_variants")
  clumpdir <- file.path(out, "clumping")
  locidir <- file.path(out, "loci")
  logs <- file.path(out, "logs")
  for (d in c(reports, extracted, clumpdir, locidir, logs)) dir.create(d)

  writeLines("INCOMPLETE_CHECK_LOGS", file.path(reports, "completion_status.txt"))
  writeLines(capture.output(sessionInfo()), file.path(reports, "R_sessionInfo.txt"))

  script_source <- file.path(base, "24_genomewide_locus_overlap.R")
  if (file.exists(script_source))
    file.copy(script_source, file.path(reports, basename(script_source)), overwrite = FALSE)

  rules <- c(
    "Step24: genome-wide significant locus discovery and cross-trait physical locus overlap.",
    "Primary traits: pancreatic fat (GCST90016675) and pancreatic iron (GCST90016676).",
    "Genome-wide significance threshold: P < 5e-8, implemented as LP > -log10(5e-8).",
    "Raw GRCh37 GWAS VCFs are rescanned; discovery is NOT restricted to the HapMap3 LDSC subset.",
    "Eligible discovery records: autosomal biallelic SNVs, FILTER PASS or '.', finite ES/SE/LP/AF, SE>0, 0<AF<1.",
    "Variants without a usable rsID remain in the GWS audit but cannot be PLINK clumping index SNPs.",
    "Duplicate rsIDs among GWS records are conservatively excluded from clumping and explicitly reported.",
    "PLINK clumping reference: existing 1000 Genomes European reference at ~/胰腺/00_ld_reference/EUR.",
    "Clumping parameters: p1=5e-8, p2=1, r2=0.001, kb=10000.",
    "Reference matching for clumping requires rsID, chromosome, base-pair position and compatible nucleotide allele set.",
    "Palindromic SNVs are allowed for LD clumping because no effect-direction inference is made in clumping.",
    "Operational locus definition: each independent lead SNP +/-500 kb; overlapping windows within a trait are merged.",
    "Physical locus overlap is descriptive and does NOT establish colocalization or a shared causal variant.",
    "Exact shared GWS variants are reported separately from overlapping loci.",
    "Existing chr8 GRCh37 region 129047378-130047378 is flagged only for context; no new chr8 mechanistic analysis is performed.",
    "No association-based post hoc selection is used beyond the prespecified genome-wide threshold."
  )
  writeLines(rules, file.path(reports, "methods_and_interpretation.txt"))

  registry <- data.frame(
    trait = names(ids),
    id = unname(ids),
    N = unname(N[names(ids)]),
    vcf = unname(vcf[names(ids)]),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(registry, file.path(reports, "study_registry.tsv"), sep = "\t")

  # 固定输入和参考文件摘要，便于后续复核。
  checksum_files <- c(
    file.path(prior19, "LDSC_paths.rds"),
    file.path(prior20, "reports", "study_registry.tsv"),
    plink_bin,
    paste0(ref_prefix, c(".bed", ".bim", ".fam")),
    vcf
  )
  checksum_files <- checksum_files[file.exists(checksum_files)]
  checksum_tab <- data.frame(
    path = checksum_files,
    size_bytes = file.info(checksum_files)$size,
    md5 = unname(tools::md5sum(checksum_files)),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(checksum_tab, file.path(reports, "input_checksums.tsv"), sep = "\t")

  # Python负责流式读取两份约926万行VCF、参考BIM核查、PLINK clumping和locus overlap。
  py <- c(
'"""Step24: pancreatic fat/iron genome-wide GWS loci and locus overlap."""',
'import sys, os, re, gzip, json, math, hashlib, traceback, subprocess, shlex',
'from pathlib import Path',
'from collections import Counter',
'import pandas as pd',
'import numpy as np',
'',
'P_THRESHOLD = 5e-8',
'LP_THRESHOLD = -math.log10(P_THRESHOLD)',
'CLUMP_R2 = 0.001',
'CLUMP_KB = 10000',
'LOCUS_HALF_WINDOW = 500000',
'CHR8_START = 129047378',
'CHR8_END = 130047378',
'',
'TRAITS = {',
'  "pancreas_fat":  {"id":"ebi-a-GCST90016675","N":25617},',
'  "pancreas_iron": {"id":"ebi-a-GCST90016676","N":25617}',
'}',
'',
'def sha256(path):',
'    h=hashlib.sha256()',
'    with open(path,"rb") as f:',
'        for b in iter(lambda:f.read(2**20),b""): h.update(b)',
'    return h.hexdigest()',
'',
'def write_table(df,path):',
'    if isinstance(df,list): df=pd.DataFrame(df)',
'    df.to_csv(path,sep="\\t",index=False,na_rep="NA")',
'',
'def parse_sample(fmt,value):',
'    keys=fmt.split(":"); vals=value.split(":")',
'    if len(keys)!=len(vals): return None',
'    return dict(zip(keys,vals))',
'',
'def numeric(x):',
'    try:',
'        v=float(x)',
'        return v if math.isfinite(v) else None',
'    except Exception:',
'        return None',
'',
'def scan_gws(vcf,trait,outdir):',
'    rows=[]; counts=Counter(); format_patterns=Counter()',
'    duplicate_ids=Counter()',
'    with gzip.open(vcf,"rt") as f:',
'        for line in f:',
'            if line.startswith("##"): continue',
'            if line.startswith("#CHROM"):',
'                header=line.rstrip("\\n").lstrip("#").split("\\t")',
'                if len(header)<10: raise ValueError("Unexpected VCF header: "+vcf)',
'                continue',
'            counts["raw_records"] += 1',
'            p=line.rstrip("\\n").split("\\t")',
'            if len(p)<10: raise ValueError("Malformed VCF line")',
'            chrom=p[0].replace("chr",""); pos=p[1]; rsid=p[2]; ref=p[3].upper(); alt=p[4].upper()',
'            filt=p[6]; fmt=p[8]; val=p[9]',
'            if chrom not in [str(i) for i in range(1,23)]: continue',
'            counts["autosomal"] += 1',
'            if len(ref)!=1 or len(alt)!=1 or ref not in "ACGT" or alt not in "ACGT": continue',
'            counts["biallelic_SNV"] += 1',
'            if filt not in ("PASS","."): continue',
'            counts["PASS_or_unfiltered"] += 1',
'            sm=parse_sample(fmt,val); format_patterns[fmt]+=1',
'            if sm is None or not all(k in sm for k in ("ES","SE","LP","AF")): continue',
'            es=numeric(sm["ES"]); se=numeric(sm["SE"]); lp=numeric(sm["LP"]); af=numeric(sm["AF"])',
'            if es is None or se is None or lp is None or af is None or se<=0 or not (0<af<1): continue',
'            counts["valid_numeric"] += 1',
'            if lp <= LP_THRESHOLD: continue',
'            counts["P_lt_5e_8"] += 1',
'            posi=int(pos)',
'            pval=10.0**(-lp) if lp < 300 else 1e-300',
'            z=es/se',
'            rows.append(dict(trait=trait,CHR=int(chrom),BP=posi,SNP=rsid,REF=ref,ALT=alt,',
'                             BETA=es,SE=se,Z=z,LP=lp,P=pval,AF=af,MAF=min(af,1-af),',
'                             FILTER=filt,FORMAT=fmt))',
'            if rsid not in ("","."): duplicate_ids[rsid]+=1',
'    if not rows: raise ValueError("No genome-wide significant variants found for "+trait)',
'    d=pd.DataFrame(rows).sort_values(["CHR","BP","SNP"]).reset_index(drop=True)',
'    d["usable_rsid"]=~d.SNP.isin(["","."])',
'    dup={k for k,v in duplicate_ids.items() if v>1}',
'    d["duplicate_rsid"]=d.SNP.isin(dup)',
'    d["eligible_for_reference_match"]=d.usable_rsid & ~d.duplicate_rsid',
'    write_table(d,outdir/(trait+"_GWS_all.tsv"))',
'    write_table(d.loc[d.duplicate_rsid],outdir/(trait+"_GWS_duplicate_rsid.tsv"))',
'    write_table(pd.DataFrame([dict(stage=k,n=v) for k,v in counts.items()]),outdir/(trait+"_scan_counts.tsv"))',
'    write_table(pd.DataFrame([dict(FORMAT=k,n=v) for k,v in format_patterns.items()]),outdir/(trait+"_FORMAT_counts.tsv"))',
'    return d,counts',
'',
'def load_bim(prefix):',
'    bim=Path(str(prefix)+".bim")',
'    x=pd.read_csv(bim,sep=r"\\s+",header=None,names=["CHR","SNP","CM","BP","A1_REF","A2_REF"],dtype={"SNP":str})',
'    x=x.loc[x.CHR.between(1,22)].copy()',
'    if x.SNP.duplicated().any():',
'        dup=set(x.loc[x.SNP.duplicated(False),"SNP"])',
'        x=x.loc[~x.SNP.isin(dup)].copy()',
'    return x.set_index("SNP")',
'',
'def complement(a):',
'    return {"A":"T","T":"A","C":"G","G":"C"}.get(a,"N")',
'',
'def allele_compatible(ref,alt,a1,a2):',
'    s={ref,alt}; t={a1,a2}',
'    if s==t: return True,"same_strand"',
'    cs={complement(ref),complement(alt)}',
'    if cs==t: return True,"complement_strand"',
'    return False,"incompatible"',
'',
'def match_reference(d,bim,trait,outdir):',
'    rows=[]',
'    for r in d.loc[d.eligible_for_reference_match].itertuples(index=False):',
'        if r.SNP not in bim.index:',
'            rows.append(dict(SNP=r.SNP,CHR=r.CHR,BP=r.BP,status="NOT_IN_REFERENCE",allele_status="NA"))',
'            continue',
'        b=bim.loc[r.SNP]',
'        if int(b.CHR)!=int(r.CHR) or int(b.BP)!=int(r.BP):',
'            rows.append(dict(SNP=r.SNP,CHR=r.CHR,BP=r.BP,status="POSITION_MISMATCH",allele_status="NA",',
'                             ref_CHR=int(b.CHR),ref_BP=int(b.BP),ref_A1=b.A1_REF,ref_A2=b.A2_REF))',
'            continue',
'        ok,ast=allele_compatible(r.REF,r.ALT,str(b.A1_REF).upper(),str(b.A2_REF).upper())',
'        rows.append(dict(SNP=r.SNP,CHR=r.CHR,BP=r.BP,status="MATCHED" if ok else "ALLELE_MISMATCH",',
'                         allele_status=ast,ref_CHR=int(b.CHR),ref_BP=int(b.BP),ref_A1=b.A1_REF,ref_A2=b.A2_REF))',
'    a=pd.DataFrame(rows)',
'    write_table(a,outdir/(trait+"_reference_match.tsv"))',
'    matched=set(a.loc[a.status=="MATCHED","SNP"]) if len(a) else set()',
'    dm=d.loc[d.SNP.isin(matched)].copy()',
'    write_table(dm,outdir/(trait+"_GWS_reference_matched.tsv"))',
'    return dm,a',
'',
'def run_clump(d,trait,plink,ref_prefix,clumpdir,logsdir):',
'    assoc=clumpdir/(trait+"_clump_input.tsv")',
'    d[["SNP","P"]].to_csv(assoc,sep="\\t",index=False)',
'    prefix=clumpdir/trait',
'    cmd=[str(plink),"--bfile",str(ref_prefix),"--clump",str(assoc),',
'         "--clump-snp-field","SNP","--clump-field","P",',
'         "--clump-p1",str(P_THRESHOLD),"--clump-p2","1",',
'         "--clump-r2",str(CLUMP_R2),"--clump-kb",str(CLUMP_KB),',
'         "--allow-no-sex","--memory","2048","--threads","1","--out",str(prefix)]',
'    log=logsdir/(trait+"_plink_clump_console.txt")',
'    with open(log,"w") as fh:',
'        rc=subprocess.run(cmd,stdout=fh,stderr=subprocess.STDOUT).returncode',
'    if rc!=0: raise RuntimeError("PLINK clumping failed for "+trait+"; see "+str(log))',
'    cf=Path(str(prefix)+".clumped")',
'    if not cf.exists(): raise RuntimeError("PLINK .clumped file missing: "+str(cf))',
'    try:',
'        c=pd.read_csv(cf,sep=r"\\s+")',
'    except pd.errors.EmptyDataError:',
'        raise RuntimeError("Empty PLINK clump output for "+trait)',
'    if "SNP" not in c.columns: raise RuntimeError("Unexpected PLINK clump columns for "+trait)',
'    leads=c[["CHR","BP","SNP","P"]].copy()',
'    leads["CHR"]=leads.CHR.astype(int);leads["BP"]=leads.BP.astype(int)',
'    leads=leads.merge(d.drop_duplicates("SNP"),on=["SNP","CHR","BP"],how="left",suffixes=("_clump",""))',
'    leads=leads.sort_values(["CHR","BP","SNP"]).reset_index(drop=True)',
'    write_table(leads,clumpdir/(trait+"_independent_leads.tsv"))',
'    return leads,cmd',
'',
'def build_loci(leads,trait):',
'    tmp=[]',
'    for r in leads.itertuples(index=False):',
'        tmp.append(dict(trait=trait,CHR=int(r.CHR),start=max(1,int(r.BP)-LOCUS_HALF_WINDOW),',
'                        end=int(r.BP)+LOCUS_HALF_WINDOW,lead_SNP=r.SNP,lead_BP=int(r.BP),lead_P=float(r.P_clump if hasattr(r,"P_clump") else r.P)))',
'    x=pd.DataFrame(tmp).sort_values(["CHR","start","end"]).reset_index(drop=True)',
'    merged=[]',
'    for chrom,g in x.groupby("CHR",sort=True):',
'        cur=None',
'        for r in g.itertuples(index=False):',
'            if cur is None or r.start>cur["end"]:',
'                if cur is not None: merged.append(cur)',
'                cur=dict(trait=trait,CHR=int(chrom),start=int(r.start),end=int(r.end),',
'                         lead_SNPs=[r.lead_SNP],lead_BPs=[int(r.lead_BP)],lead_Ps=[float(r.lead_P)])',
'            else:',
'                cur["end"]=max(cur["end"],int(r.end));cur["lead_SNPs"].append(r.lead_SNP);cur["lead_BPs"].append(int(r.lead_BP));cur["lead_Ps"].append(float(r.lead_P))',
'        if cur is not None: merged.append(cur)',
'    out=[]',
'    for i,m in enumerate(merged,1):',
'        j=int(np.argmin(m["lead_Ps"]))',
'        out.append(dict(trait=trait,locus_id=trait+"_L%02d"%i,CHR=m["CHR"],start=m["start"],end=m["end"],',
'                        n_independent_leads=len(m["lead_SNPs"]),index_SNP=m["lead_SNPs"][j],index_BP=m["lead_BPs"][j],index_P=m["lead_Ps"][j],',
'                        lead_SNPs=";".join(m["lead_SNPs"]),',
'                        overlaps_known_chr8=(m["CHR"]==8 and m["start"]<=CHR8_END and m["end"]>=CHR8_START)))',
'    return pd.DataFrame(out)',
'',
'def locus_overlap(lf,li,gws_fat,gws_iron):',
'    rows=[]',
'    cols=["fat_locus","iron_locus","CHR","overlap_start","overlap_end","overlap_bp",',
'          "fat_index_SNP","iron_index_SNP","lead_distance_bp","exact_same_index_SNP",',
'          "n_exact_shared_GWS_variants","exact_shared_GWS_variants","overlaps_known_chr8"]',
'    sf=set(gws_fat.SNP)-{"","."};si=set(gws_iron.SNP)-{"","."};shared=sf&si',
'    for a in lf.itertuples(index=False):',
'        for b in li.itertuples(index=False):',
'            if a.CHR!=b.CHR: continue',
'            s=max(a.start,b.start);e=min(a.end,b.end)',
'            if s>e: continue',
'            shared_here=[]',
'            if shared:',
'                fsub=gws_fat.loc[(gws_fat.CHR==a.CHR)&(gws_fat.BP.between(s,e))&gws_fat.SNP.isin(shared),"SNP"]',
'                isub=gws_iron.loc[(gws_iron.CHR==a.CHR)&(gws_iron.BP.between(s,e))&gws_iron.SNP.isin(shared),"SNP"]',
'                shared_here=sorted(set(fsub)&set(isub))',
'            exact_lead=(a.index_SNP==b.index_SNP)',
'            rows.append(dict(fat_locus=a.locus_id,iron_locus=b.locus_id,CHR=a.CHR,',
'                             overlap_start=s,overlap_end=e,overlap_bp=e-s+1,',
'                             fat_index_SNP=a.index_SNP,iron_index_SNP=b.index_SNP,',
'                             lead_distance_bp=abs(a.index_BP-b.index_BP),',
'                             exact_same_index_SNP=exact_lead,',
'                             n_exact_shared_GWS_variants=len(shared_here),',
'                             exact_shared_GWS_variants=";".join(shared_here),',
'                             overlaps_known_chr8=(a.CHR==8 and s<=CHR8_END and e>=CHR8_START)))',
'    return pd.DataFrame(rows,columns=cols)',
'',
'def effect_shared(g1,g2):',
'    a=g1.loc[g1.usable_rsid].drop_duplicates("SNP")',
'    b=g2.loc[g2.usable_rsid].drop_duplicates("SNP")',
'    m=a.merge(b,on="SNP",suffixes=("_fat","_iron"))',
'    if not len(m): return m',
'    same=(m.REF_fat==m.REF_iron)&(m.ALT_fat==m.ALT_iron)',
'    swapped=(m.REF_fat==m.ALT_iron)&(m.ALT_fat==m.REF_iron)',
'    m["allele_relation"]=np.where(same,"SAME_REF_ALT",np.where(swapped,"SWAPPED_REF_ALT","OTHER"))',
'    m["iron_BETA_aligned_to_fat_ALT"]=np.where(same,m.BETA_iron,np.where(swapped,-m.BETA_iron,np.nan))',
'    m["effect_direction_concordant"]=np.where(np.isfinite(m.iron_BETA_aligned_to_fat_ALT),',
'        np.sign(m.BETA_fat)==np.sign(m.iron_BETA_aligned_to_fat_ALT),np.nan)',
'    return m',
'',
'def main():',
'    base=Path(sys.argv[1]);out=Path(sys.argv[2]);plink=Path(sys.argv[3]);ref_prefix=Path(sys.argv[4]);prior20=Path(sys.argv[5])',
'    report=out/"reports";gwsdir=out/"GWS_variants";clumpdir=out/"clumping";locidir=out/"loci";logsdir=out/"logs"',
'    for d in [report,gwsdir,clumpdir,locidir,logsdir]: d.mkdir(exist_ok=True)',
'    bim=load_bim(ref_prefix)',
'    allg={};matched={};leads={};loci={};commands=[];summary=[]',
'    for trait,meta in TRAITS.items():',
'        vcf=base/(meta["id"]+".vcf.gz")',
'        g,counts=scan_gws(vcf,trait,gwsdir);allg[trait]=g',
'        gm,a=match_reference(g,bim,trait,gwsdir);matched[trait]=gm',
'        if len(gm)==0: raise ValueError("No reference-matched GWS variants for "+trait)',
'        lead,cmd=run_clump(gm,trait,plink,ref_prefix,clumpdir,logsdir);leads[trait]=lead;commands.append(dict(trait=trait,command=cmd))',
'        loc=build_loci(lead,trait);loci[trait]=loc;write_table(loc,locidir/(trait+"_loci.tsv"))',
'        refcov=(a.status=="MATCHED").sum()/max(1,len(g.loc[g.eligible_for_reference_match]))',
'        summary.append(dict(trait=trait,N=meta["N"],n_GWS=len(g),n_GWS_with_usable_unique_rsid=int(g.eligible_for_reference_match.sum()),',
'                            n_GWS_reference_matched=len(gm),reference_match_fraction=float(refcov),',
'                            n_independent_leads=len(lead),n_merged_loci=len(loc),',
'                            n_loci_overlapping_known_chr8=int(loc.overlaps_known_chr8.sum())))',
'    Path(report/"plink_commands.json").write_text(json.dumps(commands,indent=2))',
'    write_table(pd.DataFrame(summary),report/"Step24_trait_summary.tsv")',
'    shared=effect_shared(allg["pancreas_fat"],allg["pancreas_iron"])',
'    write_table(shared,locidir/"exact_shared_GWS_variants.tsv")',
'    ov=locus_overlap(loci["pancreas_fat"],loci["pancreas_iron"],allg["pancreas_fat"],allg["pancreas_iron"])',
'    write_table(ov,locidir/"fat_iron_locus_overlap.tsv")',
'    # 独有locus列表',
'    fat_hit=set(ov.fat_locus) if len(ov) else set();iron_hit=set(ov.iron_locus) if len(ov) else set()',
'    write_table(loci["pancreas_fat"].loc[~loci["pancreas_fat"].locus_id.isin(fat_hit)],locidir/"pancreas_fat_nonoverlap_loci.tsv")',
'    write_table(loci["pancreas_iron"].loc[~loci["pancreas_iron"].locus_id.isin(iron_hit)],locidir/"pancreas_iron_nonoverlap_loci.tsv")',
'    # 核心汇总',
'    core=dict(',
'      n_fat_GWS=int(len(allg["pancreas_fat"])),',
'      n_iron_GWS=int(len(allg["pancreas_iron"])),',
'      n_exact_shared_GWS_variants=int(len(shared)),',
'      n_fat_independent_leads=int(len(leads["pancreas_fat"])),',
'      n_iron_independent_leads=int(len(leads["pancreas_iron"])),',
'      n_fat_loci=int(len(loci["pancreas_fat"])),',
'      n_iron_loci=int(len(loci["pancreas_iron"])),',
'      n_physical_overlapping_locus_pairs=int(len(ov)),',
'      n_overlap_pairs_with_exact_shared_GWS_variant=int((ov.n_exact_shared_GWS_variants>0).sum()) if len(ov) else 0,',
'      n_overlap_pairs_chr8=int(ov.overlaps_known_chr8.sum()) if len(ov) else 0)',
'    write_table(pd.DataFrame([core]),report/"Step24_core_summary.tsv")',
'    # 第20步文件完整性复核：本步不使用其HM3子集发现位点，但确认既往输入目录未被改动。',
'    src=prior20/"reports"/"study_registry.tsv"',
'    if not src.exists(): raise ValueError("Step20 registry missing")',
'    (report/"step20_registry_sha256.txt").write_text(sha256(src)+"  "+str(src)+"\\n")',
'    (report/"completion_status.txt").write_text("COMPLETED_STEP24_GWS_LOCI_AND_OVERLAP_REVIEW_RESULTS\\n")',
'    print("\\nSTEP24 COMPLETED",flush=True)',
'    print(json.dumps(core,indent=2),flush=True)',
'',
'if __name__=="__main__":',
'    try:',
'        main()',
'    except Exception:',
'        out=Path(sys.argv[2]);(out/"reports").mkdir(exist_ok=True)',
'        (out/"reports"/"python_error.txt").write_text(traceback.format_exc())',
'        raise'
  )

  helper <- file.path(out, "run_step24.py")
  writeLines(py, helper, useBytes = TRUE)

  console <- file.path(logs, "python_console.txt")
  cat("第24步结果目录：", out, "\n")
  cat("将顺序扫描两份完整VCF并进行PLINK clumping；进度日志：", console, "\n")
  cat("全基因组VCF较大，运行时间可能较长。\n")
  flush.console()

  err <- NULL
  rc <- tryCatch(
    system2(
      cfg$python,
      args = vapply(
        c("-u", helper, base, out, plink_bin, ref_prefix, prior20),
        shQuote, character(1)
      ),
      stdout = console, stderr = console,
      env = c(
        "PYTHONNOUSERSITE=1",
        "PYTHONPATH=",
        "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1",
        "MKL_NUM_THREADS=1"
      )
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      127L
    }
  )

  if (!is.null(err)) writeLines(err, file.path(reports, "R_error.txt"))

  status_file <- file.path(reports, "completion_status.txt")
  status <- if (file.exists(status_file)) readLines(status_file, warn = FALSE) else "MISSING"
  if (rc != 0L || !identical(status[1], "COMPLETED_STEP24_GWS_LOCI_AND_OVERLAP_REVIEW_RESULTS")) {
    if (file.exists(console)) cat(tail(readLines(console, warn = FALSE), 50L), sep = "\n")
    stop("第24步未成功完成。请检查：", out)
  }

  # R端二次检查及简明结果展示。
  summary_file <- file.path(reports, "Step24_trait_summary.tsv")
  core_file <- file.path(reports, "Step24_core_summary.tsv")
  overlap_file <- file.path(locidir, "fat_iron_locus_overlap.tsv")
  shared_file <- file.path(locidir, "exact_shared_GWS_variants.tsv")

  if (!all(file.exists(c(summary_file, core_file, overlap_file, shared_file))))
    stop("第24步关键输出缺失。")

  tab <- data.table::fread(summary_file)
  core <- data.table::fread(core_file)
  ov <- data.table::fread(overlap_file, na.strings = "NA")
  shared <- data.table::fread(shared_file, na.strings = "NA")

  # 参考覆盖率低时仅标记，不静默删除结果。
  tab[, reference_coverage_flag :=
        ifelse(reference_match_fraction < 0.80, "REVIEW_LT_80_PERCENT", "OK")]
  data.table::fwrite(tab, summary_file, sep = "\t", na = "NA")

  # chr8区域是否至少在物理overlap中被捕获；若无，要求人工核查而不是强制结论。
  chr8_n <- if ("overlaps_known_chr8" %in% names(ov))
    sum(as.logical(ov$overlaps_known_chr8), na.rm = TRUE) else 0L
  chr8_flag <- if (chr8_n >= 1L) "EXPECTED_CHR8_CONTEXT_PRESENT" else "REVIEW_CHR8_NOT_IN_PHYSICAL_OVERLAP"
  writeLines(chr8_flag, file.path(reports, "chr8_context_check.txt"))

  # 简单结果文本，不进行“新颖性”或“共定位”推断。
  lines <- c(
    "Step24 completed: genome-wide significant loci and descriptive locus overlap.",
    paste0("Pancreas fat GWS variants: ", core$n_fat_GWS),
    paste0("Pancreas iron GWS variants: ", core$n_iron_GWS),
    paste0("Exact shared GWS rsIDs: ", core$n_exact_shared_GWS_variants),
    paste0("Independent leads (fat / iron): ",
           core$n_fat_independent_leads, " / ", core$n_iron_independent_leads),
    paste0("Merged +/-500 kb loci (fat / iron): ",
           core$n_fat_loci, " / ", core$n_iron_loci),
    paste0("Physical overlapping locus pairs: ",
           core$n_physical_overlapping_locus_pairs),
    paste0("Overlapping locus pairs containing exact shared GWS variants: ",
           core$n_overlap_pairs_with_exact_shared_GWS_variant),
    paste0("Chr8 context check: ", chr8_flag),
    "",
    "Interpretation boundary:",
    "Physical locus overlap is a screening result only. It does not establish shared causality.",
    "Prioritize overlapping loci for the subsequent local genetic correlation and cross-trait colocalization steps."
  )
  writeLines(lines, file.path(reports, "Step24_readme_results.txt"))

  # 打包小型结果供复核；不包含原始VCF、PLINK参考或其他大型资源。
  zipfile <- file.path(out, "24_locus_overlap_review.zip")
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(out)
  include <- c(
    list.files("reports", full.names = TRUE),
    list.files("loci", full.names = TRUE),
    list.files("clumping", pattern = "independent_leads.tsv$|clump_input.tsv$", full.names = TRUE),
    list.files("GWS_variants", pattern = "scan_counts.tsv$|reference_match.tsv$|GWS_all.tsv$|GWS_duplicate_rsid.tsv$", full.names = TRUE),
    list.files("logs", full.names = TRUE)
  )
  include <- include[file.exists(include)]
  utils::zip(zipfile, files = include, flags = "-q")

  cat("\n===== 第24步完成 =====\n")
  print(tab)
  cat("\n核心汇总：\n")
  print(core)
  if (nrow(ov)) {
    cat("\n脂肪—铁物理重叠locus：\n")
    print(ov)
  } else {
    cat("\n未发现按预定义±500 kb窗口形成的物理locus overlap。\n")
  }
  cat("\n复核包：", zipfile, "\n")
  cat("请保留整个输出目录，并优先上传24_locus_overlap_review.zip用于下一步核查。\n")
})()
