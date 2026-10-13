# ============================================================
# 第27步：chr19与chr8共享遗传信号的功能注释与生物学解释
# Functional annotation and biological interpretation
#
# 输入：
#   第26b步已完成的SuSiE-RSS + coloc.susie结果
#
# 目标：
#   1) 从chr19和chr8的credible sets中提取候选SNP
#   2) 识别两个表型credible set交集中的共享候选SNP
#   3) 使用Ensembl GRCh37 REST API进行VEP功能注释
#   4) 获取区域内基因并计算每个候选SNP的nearest gene
#   5) 获取Ensembl已有variant phenotype annotations（若API可返回）
#   6) 形成chr19 vs chr8的标准化比较表
#   7) 输出后续eQTL/文献/机制验证优先级表
#
# 注意：
#   - GWAS/coloc/credible-set SNP的nearest gene不是“因果基因”
#   - VEP功能后果不等于功能机制
#   - phenotype annotation仅为外部已有数据库注释，不能替代独立验证
#   - chr19来自Step25b探索性LAVA筛选，因此其功能注释属于机制跟进
#   - chr8为此前已建立的核心共享区域/阳性参照
#
# 基因组版本：
#   GRCh37 / hg19
#
# 公共API：
#   https://grch37.rest.ensembl.org
# ============================================================

(function() {

  options(stringsAsFactors = FALSE)
  base <- path.expand("~/胰腺")

  # ----------------------------------------------------------
  # 0. R packages
  # ----------------------------------------------------------
  lib <- path.expand("~/R_pancreas_step27")
  dir.create(lib, recursive=TRUE, showWarnings=FALSE)
  .libPaths(c(lib, .libPaths()))

  need <- function(pkg) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
      install.packages(pkg, lib=lib, repos="https://cloud.r-project.org")
    }
    if (!requireNamespace(pkg, quietly=TRUE)) {
      stop("无法加载R包：", pkg)
    }
  }

  for (p in c("data.table","httr","jsonlite")) need(p)

  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(httr))
  suppressPackageStartupMessages(library(jsonlite))

  # ----------------------------------------------------------
  # 1. 定位最新完成的Step26b
  # ----------------------------------------------------------
  p26b <- file.path(base, "26b_coloc_susie_fix")
  if (!dir.exists(p26b)) stop("未找到第26b目录：", p26b)

  d26b <- list.dirs(p26b, recursive=FALSE, full.names=TRUE)
  if (!length(d26b)) stop("第26b目录为空。")

  completed <- d26b[vapply(d26b, function(d) {
    f <- file.path(d, "reports", "completion_status.txt")
    file.exists(f) &&
      identical(
        readLines(f, warn=FALSE)[1],
        "COMPLETED_STEP26B_SUSIE_COLOCALIZATION"
      )
  }, logical(1))]

  if (!length(completed)) {
    stop("未找到完成标记为COMPLETED_STEP26B_SUSIE_COLOCALIZATION的第26b目录。")
  }

  step26b <- completed[which.max(file.info(completed)$mtime)]
  cat("使用Step26b目录：\n", step26b, "\n")

  # ----------------------------------------------------------
  # 2. 输出目录
  # ----------------------------------------------------------
  parent <- file.path(base, "27_functional_annotation")
  dir.create(parent, recursive=TRUE, showWarnings=FALSE)

  out <- tempfile("chr19_chr8_annotation_", tmpdir=parent)
  if (!dir.create(out)) stop("无法创建Step27输出目录。")

  reports <- file.path(out, "reports")
  rawdir  <- file.path(out, "raw_api")
  tables  <- file.path(out, "tables")
  logs    <- file.path(out, "logs")

  for (d in c(reports, rawdir, tables, logs))
    dir.create(d, recursive=TRUE, showWarnings=FALSE)

  writeLines(
    "INCOMPLETE_CHECK_LOGS",
    file.path(reports, "completion_status.txt")
  )

  # ----------------------------------------------------------
  # 3. 固定区域定义
  # ----------------------------------------------------------
  regions <- data.table(
    region_id=c(
      "chr19_LOC2340_primary",
      "chr8_positive_control"
    ),
    CHR=c(19L, 8L),
    START=c(33785836L, 129047378L),
    STOP=c(34633274L, 130047378L),
    role=c(
      "EXPLORATORY_SECONDARY_SHARED_LOCUS",
      "PRIMARY_SHARED_LOCUS_POSITIVE_CONTROL"
    )
  )

  fwrite(
    regions,
    file.path(tables, "Step27_regions.tsv"),
    sep="\t"
  )

  # ----------------------------------------------------------
  # 4. 读取credible sets
  # ----------------------------------------------------------
  read_cs <- function(region_id, trait) {

    f <- file.path(
      step26b,
      "susie",
      region_id,
      paste0(trait, "_credible_sets.tsv")
    )

    if (!file.exists(f)) {
      stop("缺少credible-set文件：", f)
    }

    x <- fread(f)

    req <- c("cs","SNP","BP")
    if (!all(req %in% names(x))) {
      stop(
        "credible-set文件字段异常：", f,
        "\n实际字段：", paste(names(x), collapse=",")
      )
    }

    x[, region_id := region_id]
    x[, trait := trait]

    x
  }

  cs_list <- list()

  for (rid in regions$region_id) {
    cs_list[[paste(rid,"fat",sep="__")]] <-
      read_cs(rid, "pancreas_fat")
    cs_list[[paste(rid,"iron",sep="__")]] <-
      read_cs(rid, "pancreas_iron")
  }

  cs_all <- rbindlist(cs_list, fill=TRUE)
  cs_all <- unique(cs_all)

  fwrite(
    cs_all,
    file.path(tables, "Step27_all_credible_set_variants_long.tsv"),
    sep="\t",
    na="NA"
  )

  # ----------------------------------------------------------
  # 5. 汇总每个SNP的credible-set证据
  # ----------------------------------------------------------
  snp_summary <- cs_all[
    ,
    .(
      BP=min(BP, na.rm=TRUE),
      traits=paste(sort(unique(trait)), collapse=";"),
      n_traits=uniqueN(trait),
      cs_fat=paste(sort(unique(cs[trait=="pancreas_fat"])), collapse=";"),
      cs_iron=paste(sort(unique(cs[trait=="pancreas_iron"])), collapse=";"),
      n_cs_fat=uniqueN(cs[trait=="pancreas_fat"]),
      n_cs_iron=uniqueN(cs[trait=="pancreas_iron"])
    ),
    by=.(region_id,SNP)
  ]

  snp_summary[, shared_between_traits := n_traits==2L]

  snp_summary <- merge(
    snp_summary,
    regions[,.(region_id,CHR,START,STOP,role)],
    by="region_id",
    all.x=TRUE
  )

  setorder(
    snp_summary,
    region_id,
    -shared_between_traits,
    BP
  )

  fwrite(
    snp_summary,
    file.path(tables, "Step27_candidate_variants.tsv"),
    sep="\t",
    na="NA"
  )

  shared_snps <- snp_summary[
    shared_between_traits==TRUE
  ]

  fwrite(
    shared_snps,
    file.path(tables, "Step27_shared_credible_set_variants.tsv"),
    sep="\t",
    na="NA"
  )

  # ----------------------------------------------------------
  # 6. Ensembl REST helper
  # ----------------------------------------------------------
  server <- "https://grch37.rest.ensembl.org"

  api_get_json <- function(path, query=list(), tries=4L) {

    for (i in seq_len(tries)) {

      r <- tryCatch(
        GET(
          paste0(server, path),
          query=query,
          add_headers(
            `Content-Type`="application/json",
            Accept="application/json"
          ),
          timeout(60)
        ),
        error=function(e) e
      )

      if (!inherits(r,"error") &&
          status_code(r)>=200 &&
          status_code(r)<300) {

        return(
          content(
            r,
            as="parsed",
            type="application/json",
            simplifyVector=FALSE
          )
        )
      }

      Sys.sleep(2*i)
    }

    NULL
  }

  api_post_json <- function(path, body, tries=4L) {

    for (i in seq_len(tries)) {

      r <- tryCatch(
        POST(
          paste0(server, path),
          body=body,
          encode="json",
          add_headers(
            `Content-Type`="application/json",
            Accept="application/json"
          ),
          timeout(120)
        ),
        error=function(e) e
      )

      if (!inherits(r,"error") &&
          status_code(r)>=200 &&
          status_code(r)<300) {

        return(
          content(
            r,
            as="parsed",
            type="application/json",
            simplifyVector=FALSE
          )
        )
      }

      Sys.sleep(2*i)
    }

    NULL
  }

  # ----------------------------------------------------------
  # 7. Ensembl VEP annotation
  # ----------------------------------------------------------
  rsids <- unique(
    snp_summary[
      grepl("^rs[0-9]+$", SNP),
      SNP
    ]
  )

  cat("待VEP注释rsID数量：", length(rsids), "\n")

  # VEP POST限制300 IDs；保守每批200
  chunks <- split(
    rsids,
    ceiling(seq_along(rsids)/200)
  )

  vep_raw <- list()

  for (i in seq_along(chunks)) {

    cat("VEP batch ", i, "/", length(chunks), "\n", sep="")

    res <- api_post_json(
      "/vep/human/id",
      list(ids=unname(chunks[[i]]))
    )

    if (!is.null(res)) {
      vep_raw <- c(vep_raw, res)
    }

    Sys.sleep(0.5)
  }

  saveRDS(
    vep_raw,
    file.path(rawdir, "Ensembl_GRCh37_VEP_raw.rds")
  )

  # ----------------------------------------------------------
  # 8. flatten transcript consequences
  # ----------------------------------------------------------
  tc_rows <- list()
  var_rows <- list()

  get1 <- function(x,nm,default=NA_character_) {
    if (!is.null(x[[nm]]) && length(x[[nm]])>0) {
      paste(unlist(x[[nm]]), collapse=";")
    } else {
      default
    }
  }

  for (v in vep_raw) {

    vid <- get1(v,"id")
    most_severe <- get1(v,"most_severe_consequence")
    assembly <- get1(v,"assembly_name")
    seq_region <- get1(v,"seq_region_name")
    start <- suppressWarnings(as.numeric(get1(v,"start",NA)))
    end <- suppressWarnings(as.numeric(get1(v,"end",NA)))
    allele_string <- get1(v,"allele_string")

    var_rows[[length(var_rows)+1L]] <- data.table(
      SNP=vid,
      assembly=assembly,
      CHR=seq_region,
      start=start,
      end=end,
      allele_string=allele_string,
      most_severe_consequence=most_severe
    )

    tcs <- v[["transcript_consequences"]]

    if (!is.null(tcs) && length(tcs)) {

      for (tc in tcs) {

        tc_rows[[length(tc_rows)+1L]] <- data.table(
          SNP=vid,
          gene_id=get1(tc,"gene_id"),
          gene_symbol=get1(tc,"gene_symbol"),
          transcript_id=get1(tc,"transcript_id"),
          biotype=get1(tc,"biotype"),
          consequence_terms=get1(tc,"consequence_terms"),
          canonical=get1(tc,"canonical"),
          impact=get1(tc,"impact"),
          hgvsc=get1(tc,"hgvsc"),
          hgvsp=get1(tc,"hgvsp"),
          distance=get1(tc,"distance"),
          strand=get1(tc,"strand")
        )
      }
    }
  }

  vep_variants <- if (length(var_rows))
    rbindlist(var_rows, fill=TRUE) else data.table()

  vep_tc <- if (length(tc_rows))
    rbindlist(tc_rows, fill=TRUE) else data.table()

  fwrite(
    vep_variants,
    file.path(tables, "Step27_Ensembl_VEP_variant_summary.tsv"),
    sep="\t",
    na="NA"
  )

  fwrite(
    vep_tc,
    file.path(tables, "Step27_Ensembl_VEP_transcript_consequences.tsv"),
    sep="\t",
    na="NA"
  )

  # ----------------------------------------------------------
  # 9. 获取区域内gene annotations
  # ----------------------------------------------------------
  gene_rows <- list()

  for (i in seq_len(nrow(regions))) {

    rr <- regions[i]

    region_string <- paste0(
      rr$CHR, ":",
      rr$START, "-",
      rr$STOP
    )

    res <- api_get_json(
      paste0(
        "/overlap/region/human/",
        region_string
      ),
      query=list(feature="gene")
    )

    if (is.null(res)) next

    for (g in res) {

      gene_rows[[length(gene_rows)+1L]] <- data.table(
        region_id=rr$region_id,
        gene_id=get1(g,"id"),
        gene_symbol=get1(g,"external_name"),
        gene_start=as.numeric(get1(g,"start",NA)),
        gene_end=as.numeric(get1(g,"end",NA)),
        strand=get1(g,"strand"),
        biotype=get1(g,"biotype"),
        description=get1(g,"description")
      )
    }

    Sys.sleep(0.5)
  }

  genes <- if (length(gene_rows))
    unique(rbindlist(gene_rows, fill=TRUE)) else data.table()

  fwrite(
    genes,
    file.path(tables, "Step27_region_genes_GRCh37.tsv"),
    sep="\t",
    na="NA"
  )

  # ----------------------------------------------------------
  # 10. nearest gene
  # ----------------------------------------------------------
  nearest_rows <- list()

  if (nrow(genes)) {

    for (i in seq_len(nrow(snp_summary))) {

      s <- snp_summary[i]
      gg <- genes[
        region_id==s$region_id &
        is.finite(gene_start) &
        is.finite(gene_end)
      ]

      if (!nrow(gg)) next

      gg[, distance_to_gene :=
        fifelse(
          s$BP < gene_start,
          gene_start-s$BP,
          fifelse(
            s$BP > gene_end,
            s$BP-gene_end,
            0
          )
        )
      ]

      md <- min(gg$distance_to_gene, na.rm=TRUE)

      nearest <- gg[
        distance_to_gene==md
      ]

      nearest_rows[[length(nearest_rows)+1L]] <-
        nearest[
          ,
          .(
            region_id=s$region_id,
            SNP=s$SNP,
            BP=s$BP,
            nearest_gene_id=gene_id,
            nearest_gene_symbol=gene_symbol,
            nearest_gene_biotype=biotype,
            distance_to_gene=distance_to_gene
          )
        ]
    }
  }

  nearest <- if (length(nearest_rows))
    rbindlist(nearest_rows, fill=TRUE) else data.table()

  fwrite(
    nearest,
    file.path(tables, "Step27_variant_nearest_genes.tsv"),
    sep="\t",
    na="NA"
  )

  # ----------------------------------------------------------
  # 11. Variant phenotype associations
  #     非关键模块，API失败不使Step27失败
  # ----------------------------------------------------------
  pheno_rows <- list()

  top_for_pheno <- unique(
    shared_snps[
      grepl("^rs[0-9]+$",SNP),
      SNP
    ]
  )

  if (length(top_for_pheno)) {

    pchunks <- split(
      top_for_pheno,
      ceiling(seq_along(top_for_pheno)/100)
    )

    for (i in seq_along(pchunks)) {

      res <- api_post_json(
        "/variation/human",
        list(ids=unname(pchunks[[i]]))
      )

      if (is.null(res)) next

      # POST variation returns named list keyed by rsID in many Ensembl versions
      if (!is.null(names(res))) {

        for (vid in names(res)) {

          obj <- res[[vid]]
          pa <- obj[["phenotypes"]]
          if (is.null(pa)) pa <- obj[["phenotype_associations"]]

          if (!is.null(pa) && length(pa)) {

            for (p in pa) {

              pheno_rows[[length(pheno_rows)+1L]] <- data.table(
                SNP=vid,
                phenotype=get1(p,"trait"),
                source=get1(p,"source"),
                study=get1(p,"study"),
                risk_allele=get1(p,"risk_allele"),
                associated_gene=get1(p,"associated_gene"),
                p_value=get1(p,"p_value")
              )
            }
          }
        }
      }

      Sys.sleep(0.5)
    }
  }

  phenos <- if (length(pheno_rows))
    rbindlist(pheno_rows, fill=TRUE) else data.table(
      SNP=character(),
      phenotype=character(),
      source=character(),
      study=character(),
      risk_allele=character(),
      associated_gene=character(),
      p_value=character()
    )

  fwrite(
    phenos,
    file.path(tables, "Step27_Ensembl_variant_phenotype_annotations.tsv"),
    sep="\t",
    na="NA"
  )

  # ----------------------------------------------------------
  # 12. 每个SNP的简化功能注释
  # ----------------------------------------------------------
  annot <- copy(snp_summary)

  if (nrow(vep_variants)) {
    annot <- merge(
      annot,
      vep_variants[
        ,
        .(
          SNP,
          most_severe_consequence
        )
      ],
      by="SNP",
      all.x=TRUE
    )
  } else {
    annot[, most_severe_consequence := NA_character_]
  }

  if (nrow(nearest)) {

    nearest_collapsed <- nearest[
      ,
      .(
        nearest_gene=paste(
          sort(unique(na.omit(nearest_gene_symbol))),
          collapse=";"
        ),
        nearest_gene_distance=min(
          distance_to_gene,
          na.rm=TRUE
        )
      ),
      by=.(region_id,SNP)
    ]

    annot <- merge(
      annot,
      nearest_collapsed,
      by=c("region_id","SNP"),
      all.x=TRUE
    )

  } else {

    annot[, nearest_gene := NA_character_]
    annot[, nearest_gene_distance := NA_real_]
  }

  # transcript-level genes touched by VEP
  if (nrow(vep_tc)) {

    tc_gene <- vep_tc[
      !is.na(gene_symbol) &
      gene_symbol!="",
      .(
        VEP_genes=paste(
          sort(unique(gene_symbol)),
          collapse=";"
        ),
        VEP_consequences=paste(
          sort(unique(consequence_terms)),
          collapse=";"
        )
      ),
      by=SNP
    ]

    annot <- merge(
      annot,
      tc_gene,
      by="SNP",
      all.x=TRUE
    )

  } else {

    annot[, VEP_genes := NA_character_]
    annot[, VEP_consequences := NA_character_]
  }

  # ----------------------------------------------------------
  # 13. 证据标签（不是“因果评分”）
  # ----------------------------------------------------------
  annot[, evidence_class :=
    fifelse(
      shared_between_traits &
        region_id=="chr19_LOC2340_primary",
      "CHR19_SHARED_CREDIBLE_SET_CANDIDATE",
      fifelse(
        shared_between_traits &
          region_id=="chr8_positive_control",
        "CHR8_SHARED_CREDIBLE_SET_CANDIDATE",
        "TRAIT_SPECIFIC_CREDIBLE_SET_MEMBER"
      )
    )
  ]

  annot[, followup_priority :=
    fifelse(
      shared_between_traits,
      "HIGH",
      "SECONDARY"
    )
  ]

  annot[, gwas_catalog_search :=
    paste0(
      "https://www.ebi.ac.uk/gwas/search?query=",
      SNP
    )
  ]

  annot[, pubmed_search :=
    paste0(
      "https://pubmed.ncbi.nlm.nih.gov/?term=",
      SNP
    )
  ]

  annot[, regulomedb_search :=
    paste0(
      "https://regulomedb.org/regulome-search/?regions=",
      SNP
    )
  ]

  setorder(
    annot,
    region_id,
    -shared_between_traits,
    BP
  )

  fwrite(
    annot,
    file.path(tables, "Step27_integrated_candidate_annotation.tsv"),
    sep="\t",
    na="NA"
  )

  # ----------------------------------------------------------
  # 14. 区域级比较
  # ----------------------------------------------------------
  status26b <- fread(
    file.path(
      step26b,
      "reports",
      "Step26b_ABF_plus_SuSiE_summary.tsv"
    )
  )

  region_compare <- annot[
    ,
    .(
      n_credible_set_variants=uniqueN(SNP),
      n_shared_between_traits=uniqueN(
        SNP[shared_between_traits==TRUE]
      ),
      n_trait_specific=uniqueN(
        SNP[shared_between_traits==FALSE]
      ),
      n_VEP_annotated=sum(
        !is.na(most_severe_consequence)
      ),
      n_region_genes=genes[
        region_id==.BY$region_id,
        uniqueN(gene_id)
      ]
    ),
    by=region_id
  ]

  # 尝试拼接ABF/SuSiE summary
  cols_keep <- intersect(
    c(
      "region_id",
      "ABF_PP_H3",
      "ABF_PP_H4",
      "n_cs_pancreas_fat",
      "n_cs_pancreas_iron",
      "n_coloc_susie_signal_pairs",
      "max_PP_H4_susie"
    ),
    names(status26b)
  )

  if (length(cols_keep)>1) {
    region_compare <- merge(
      region_compare,
      status26b[, ..cols_keep],
      by="region_id",
      all.x=TRUE
    )
  }

  region_compare <- merge(
    region_compare,
    regions[,.(region_id,role)],
    by="region_id",
    all.x=TRUE
  )

  fwrite(
    region_compare,
    file.path(tables, "Step27_chr19_vs_chr8_comparison.tsv"),
    sep="\t",
    na="NA"
  )

  # ----------------------------------------------------------
  # 15. 重点候选SNP表
  # ----------------------------------------------------------
  priority <- annot[
    shared_between_traits==TRUE
  ]

  # chr19明确已知的两个Step26b共享signal pair SNP作为核对标记
  priority[, step26b_known_shared_signal :=
    SNP %in% c("rs430752","rs33824")
  ]

  # chr8此前反复出现的核心候选作为核对标记
  priority[, prior_chr8_core_candidate :=
    SNP %in% c(
      "rs4733612",
      "rs1561927",
      "rs7824269",
      "rs7819154",
      "rs12676584"
    )
  ]

  setorder(
    priority,
    region_id,
    -step26b_known_shared_signal,
    -prior_chr8_core_candidate,
    BP
  )

  fwrite(
    priority,
    file.path(tables, "Step27_priority_shared_candidates.tsv"),
    sep="\t",
    na="NA"
  )

  # ----------------------------------------------------------
  # 16. 自动生成结果摘要
  # ----------------------------------------------------------
  n19_all <- annot[
    region_id=="chr19_LOC2340_primary",
    uniqueN(SNP)
  ]
  n19_shared <- annot[
    region_id=="chr19_LOC2340_primary" &
      shared_between_traits==TRUE,
    uniqueN(SNP)
  ]
  n8_all <- annot[
    region_id=="chr8_positive_control",
    uniqueN(SNP)
  ]
  n8_shared <- annot[
    region_id=="chr8_positive_control" &
      shared_between_traits==TRUE,
    uniqueN(SNP)
  ]

  lines <- c(
    "Step27 functional annotation summary",
    "",
    paste0("Step26b source: ", step26b),
    "Genome build: GRCh37/hg19",
    "",
    paste0(
      "chr19 LOC2340 credible-set variants: ",
      n19_all,
      "; shared between pancreatic fat and pancreatic iron credible sets: ",
      n19_shared,
      "."
    ),
    paste0(
      "chr8 positive-control credible-set variants: ",
      n8_all,
      "; shared between pancreatic fat and pancreatic iron credible sets: ",
      n8_shared,
      "."
    ),
    "",
    "Interpretation rules:",
    "1. Shared credible-set membership is evidence prioritization, not proof that the SNP is causal.",
    "2. Nearest-gene assignment is descriptive only and must not be interpreted as causal-gene assignment.",
    "3. Ensembl VEP annotations describe predicted molecular consequence and do not establish a tissue-specific mechanism.",
    "4. chr19 remains an exploratory secondary shared locus because it was selected in the Step25b sensitivity analysis.",
    "5. chr8 remains the primary shared locus/positive control because its evidence is stronger and was robust across earlier fine-mapping and colocalization analyses.",
    "",
    "Next use:",
    "Use Step27_priority_shared_candidates.tsv to choose variants/genes for tissue-specific eQTL follow-up, literature review, and regulatory annotation."
  )

  writeLines(
    lines,
    file.path(reports, "Step27_summary.txt")
  )

  # ----------------------------------------------------------
  # 17. QC
  # ----------------------------------------------------------
  qc <- data.table(
    metric=c(
      "n_regions",
      "n_all_unique_candidate_SNPs",
      "n_shared_candidate_SNPs",
      "n_VEP_variant_records",
      "n_VEP_transcript_records",
      "n_region_gene_records",
      "n_nearest_gene_records",
      "n_phenotype_annotation_records"
    ),
    value=c(
      nrow(regions),
      uniqueN(snp_summary$SNP),
      uniqueN(shared_snps$SNP),
      nrow(vep_variants),
      nrow(vep_tc),
      nrow(genes),
      nrow(nearest),
      nrow(phenos)
    )
  )

  fwrite(
    qc,
    file.path(reports, "Step27_QC_summary.tsv"),
    sep="\t"
  )

  # VEP可以部分失败，但至少要求共享SNP表和区域gene结果完整
  serious <- character()

  if (!nrow(snp_summary))
    serious <- c(serious, "NO_CREDIBLE_SET_VARIANTS")

  if (!nrow(shared_snps))
    serious <- c(serious, "NO_SHARED_CREDIBLE_SET_VARIANTS")

  if (!nrow(genes))
    serious <- c(serious, "NO_REGION_GENES_FROM_ENSEMBL")

  if (length(serious)) {
    writeLines(
      "COMPLETED_STEP27_WITH_REVIEW_FLAGS",
      file.path(reports, "completion_status.txt")
    )
    writeLines(
      serious,
      file.path(reports, "review_flags.txt")
    )
  } else {
    writeLines(
      "COMPLETED_STEP27_FUNCTIONAL_ANNOTATION",
      file.path(reports, "completion_status.txt")
    )
  }

  # ----------------------------------------------------------
  # 18. 打包精简review文件
  # ----------------------------------------------------------
  zipfile <- file.path(
    out,
    "27_functional_annotation_review.zip"
  )

  oldwd <- getwd()
  on.exit(setwd(oldwd), add=TRUE)
  setwd(out)

  files_to_zip <- c(
    list.files("reports", recursive=TRUE, full.names=TRUE),
    list.files("tables", recursive=TRUE, full.names=TRUE)
  )

  files_to_zip <- files_to_zip[
    file.exists(files_to_zip) &
      !dir.exists(files_to_zip)
  ]

  utils::zip(
    zipfile,
    files=files_to_zip,
    flags="-q"
  )

  # ----------------------------------------------------------
  # 19. 控制台输出
  # ----------------------------------------------------------
  cat("\n========================================\n")
  cat("第27步完成\n")
  cat("========================================\n\n")

  print(qc)

  cat("\nchr19 / chr8比较：\n")
  print(region_compare)

  cat("\n优先共享候选SNP（前30）：\n")
  print(head(priority,30))

  cat("\n结果目录：", out, "\n")
  cat("复核包：", zipfile, "\n")

})()
