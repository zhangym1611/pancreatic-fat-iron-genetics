#!/usr/bin/env python3
"""Static integrity/identifiers only. Does not execute or validate R code."""
from pathlib import Path
import csv,hashlib,sys
root=Path(__file__).resolve().parents[1]
manifest=root/'00_Documentation'/'SHA256_MANIFEST.csv'
assert manifest.exists(), 'Missing SHA256 manifest'
issues=[]; total=0
with manifest.open(encoding='utf8',newline='') as f:
    for row in csv.DictReader(f):
        total+=1;p=root/row['relative_path']
        if not p.is_file():issues.append('MISSING '+row['relative_path']);continue
        digest=hashlib.sha256(p.read_bytes()).hexdigest()
        if digest!=row['sha256']:issues.append('HASH MISMATCH '+row['relative_path'])
up=root/'01_Analysis_scripts'/'06_Glycaemic_followup'/'Step30a_glycaemic_ABF_colocalisation.R'
down=root/'01_Analysis_scripts'/'06_Glycaemic_followup'/'Step30b_glycaemic_SuSiE_colocalisation.R'
if not(up.is_file() and down.is_file()): issues.append('Step30a/30b absent')
else:
    a=up.read_text(encoding='utf8');b=down.read_text(encoding='utf8')
    for tok in ['COMPLETED_STEP30A_ABF_GLYCEMIC_LOCUS_FOLLOWUP','Step30a_SuSiE_followup_candidates.tsv','30a_glycaemic_locus_followup']:
        if tok not in a:issues.append('UPSTREAM MISSING: '+tok)
        if tok not in b:issues.append('DOWNSTREAM MISSING: '+tok)
    if 'any(grepl("^COMPLETED_STEP30A", x))' in b:issues.append('Unsafe prefix-based completed-run selection remains')
    if 'PANCREAS_STEP30A_RUN_DIR' not in b:issues.append('Missing upstream run-pinning option')
    if '31_glycemic_locus_followup' in a+b:issues.append('Unconverted old Step31 output path')
if issues:
    print('STATIC PACKAGE CHECK FAILED:')
    for issue in issues: print(' -',issue)
    sys.exit(1)
print(f'STATIC PACKAGE CHECK PASSED ({total} files match SHA256; Step30a/30b marker linkage present).')
print('This does NOT test R syntax, execution, scientific estimates, reference availability or reproducibility.')
