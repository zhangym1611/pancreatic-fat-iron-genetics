# Run the entire script in the remote RStudio session.
# Reads original VCFs; writes only to a new audit directory.
(function() {
  pkgs <- c('data.table', 'Rsamtools', 'GenomicRanges', 'IRanges', 'coloc')
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) stop('Missing packages: ', paste(pkgs[!ok], collapse = ', '))
  base <- path.expand('~/胰腺')
  ids <- c(fat = 'ebi-a-GCST90016675', iron = 'ebi-a-GCST90016676')
  paths <- file.path(base, paste0(ids, '.vcf.gz'))
  if (!all(file.exists(c(paths, paste0(paths, '.tbi'))))) stop('VCF or TBI missing')
  parent <- file.path(base, '07_validation')
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  dest <- tempfile('chr8_audit_', tmpdir = parent)
  if (!dir.create(dest)) stop('Cannot create output directory')
  sink(file.path(dest, 'run_log.txt'), split = TRUE)
  on.exit(sink(), add = TRUE)
  cat('Output:', dest, '\n')
  cat('This is a robustness audit, not independent replication.\n')
  center <- 129547378L
  widths <- c(250000L, 500000L, 1000000L)
  read_region <- function(f, id) {
    con <- gzfile(f, 'rt')
    h <- tryCatch(readLines(con, n = 250L, warn = FALSE), finally = close(con))
    sample <- h[startsWith(h, '##SAMPLE=')]
    cat('\n', id, '\n', paste(sample, collapse = '\n'), '\n')
    if (length(sample) != 1L || !grepl(id, sample, fixed = TRUE) ||
        !grepl('TotalControls=25617', sample, fixed = TRUE)) stop('Unexpected study metadata')
    if (!any(grepl('assembly=HG19/GRCh37', h, fixed = TRUE))) stop('Check genome build')
    tb <- Rsamtools::TabixFile(f, index = paste0(f, '.tbi'))
    open(tb)
    on.exit(close(tb), add = TRUE)
    reg <- GenomicRanges::GRanges('8', IRanges::IRanges(center - max(widths), center + max(widths)))
    lines <- Rsamtools::scanTabix(tb, param = reg)[[1L]]
    if (!length(lines)) stop('No records returned')
    d <- data.table::fread(text = paste(lines, collapse = '\n'), sep = '\t',
                          header = FALSE, colClasses = 'character', na.strings = '.',
                          quote = '', showProgress = FALSE)
    if (ncol(d) != 10L) stop('Expected one-sample, ten-column VCF')
    data.table::setnames(d, c('chr','pos','rsid','REF','ALT','QUAL','FILTER','INFO','FORMAT','VALUE'))
    d$pos <- as.integer(d$pos)
    for (nm in c('ES','SE','LP','AF')) d[[nm]] <- NA_real_
    if (anyNA(d$FORMAT)) stop('Missing FORMAT')
    for (fmt in unique(d$FORMAT)) {
      ix <- which(d$FORMAT == fmt)
      keys <- strsplit(fmt, ':', fixed = TRUE)[[1L]]
      if (anyDuplicated(keys) || !all(c('ES','SE','LP','AF') %in% keys)) stop('Unsupported FORMAT: ', fmt)
      parts <- strsplit(d$VALUE[ix], ':', fixed = TRUE)
      if (any(lengths(parts) != length(keys))) stop('Unexpected sample field length')
      for (nm in c('ES','SE','LP','AF')) {
        v <- vapply(parts, `[`, character(1), match(nm, keys))
        v[v == '.'] <- NA_character_
        d[[nm]][ix] <- suppressWarnings(as.numeric(v))
      }
    }
    valid <- d$chr == '8' & is.finite(d$pos) &
      d$pos >= center - max(widths) & d$pos <= center + max(widths) &
      grepl('^[ACGT]+$', d$REF) & grepl('^[ACGT]+$', d$ALT) & d$REF != d$ALT &
      is.finite(d$ES) & is.finite(d$SE) & d$SE > 0 &
      is.finite(d$LP) & d$LP >= 0 & is.finite(d$AF) & d$AF > 0 & d$AF < 1 &
      (is.na(d$FILTER) | d$FILTER == 'PASS')
    valid[is.na(valid)] <- FALSE
    cat('Raw:', nrow(d), ' valid:', sum(valid), ' excluded:', sum(!valid), '\n')
    data.table::fwrite(d[!valid], file.path(dest, paste0(id, '_excluded.tsv')), sep = '\t')
    d <- d[valid]
    d$key <- paste(d$chr, d$pos, d$REF, d$ALT, sep = ':')
    if (anyDuplicated(d$key)) stop('Duplicate variant keys; inspect before continuing')
    d
  }
  fat <- read_region(paths[1L], ids[1L])
  iron <- read_region(paths[2L], ids[2L])
  k <- intersect(fat$key, iron$key)
  cat('\nShared variants:', length(k), '\n')
  cat('Fat-only:', sum(!fat$key %in% k), ' Iron-only:', sum(!iron$key %in% k), '\n')
  a <- as.data.frame(fat[match(k, fat$key)])
  b <- as.data.frame(iron[match(k, iron$key)])
  stopifnot(identical(a$key, b$key), identical(a$pos, b$pos))
  if (length(k) < 100L) stop('Too few shared variants; inspect input')
  saveRDS(list(fat = a, iron = b, N_assumed = 25617L, center = center,
               widths = widths), file.path(dest, 'matched_region_inputs.rds'))
  cat('Maximum paired AF difference:', max(abs(a$AF - b$AF)), '\n')
  # Dense regional coverage: no P-value, LD or positive-MAF thresholding.
  dataset <- function(d, fixed_sd) {
    x <- list(snp = d$key, position = d$pos, beta = d$ES, varbeta = d$SE^2,
              MAF = pmin(d$AF, 1-d$AF), N = 25617L, type = 'quant')
    if (fixed_sd) x$sdY <- 1
    x
  }
  results <- list(); fits <- list(); warning_log <- character()
  for (w in widths) for (fixed_sd in c(FALSE, TRUE)) for (p12 in c(1e-6, 1e-5, 5e-5)) {
    ix <- abs(a$pos - center) <= w
    label <- paste(w, if (fixed_sd) 'sdY1_sensitivity' else 'sdY_estimated', p12, sep = '_')
    cat('\nRunning:', label, ' n=', sum(ix), '\n')
    fit <- withCallingHandlers(
      coloc::coloc.abf(dataset(a[ix, ], fixed_sd), dataset(b[ix, ], fixed_sd),
                       p1 = 1e-4, p2 = 1e-4, p12 = p12),
      warning = function(e) {
        warning_log <<- c(warning_log, paste(label, conditionMessage(e), sep = ': '))
        invokeRestart('muffleWarning')
      })
    s <- fit$summary
    results[[label]] <- data.frame(half_window_bp = w,
      sdY_mode = if (fixed_sd) 'sdY1_sensitivity' else 'sdY_estimated',
      p1 = 1e-4, p2 = 1e-4, p12 = p12, nsnps = unname(s['nsnps']),
      PP_H0 = unname(s['PP.H0.abf']), PP_H1 = unname(s['PP.H1.abf']),
      PP_H2 = unname(s['PP.H2.abf']), PP_H3 = unname(s['PP.H3.abf']),
      PP_H4 = unname(s['PP.H4.abf']))
    fits[[label]] <- fit
    if (w == 500000L && !fixed_sd && p12 == 1e-5) {
      v <- fit$results[order(fit$results$SNP.PP.H4, decreasing = TRUE), ]
      v$rsid <- a$rsid[match(v$snp, a$key)]
      v$cumulative_PP_H4 <- cumsum(v$SNP.PP.H4)
      data.table::fwrite(v, file.path(dest, 'baseline_variant_posteriors.tsv'), sep = '\t')
      cat('Baseline H4:', unname(s['PP.H4.abf']), ' prior log value: about 0.99562\n')
      print(head(v[, c('rsid','snp','SNP.PP.H4','cumulative_PP_H4')], 10L))
    }
  }
  summary <- data.table::rbindlist(results)
  data.table::fwrite(summary, file.path(dest, 'coloc_validation_summary.tsv'), sep = '\t')
  saveRDS(fits, file.path(dest, 'coloc_validation_fits.rds'))
  writeLines(unique(warning_log), file.path(dest, 'warnings.txt'))
  # Diagnostic vector PDF; no LD colour scheme is claimed.
  plot_region <- function() {
    grDevices::pdf(file.path(dest, 'regional_association_diagnostic.pdf'), width = 9, height = 7)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mfrow = c(2,1), mar = c(4,4,2,1))
    for (z in c('fat','iron')) {
      d <- if (z == 'fat') a else b
      graphics::plot(d$pos/1e6, d$LP, pch = 16, cex = 0.45, col = '#32688A',
        xlab = 'Chromosome 8 position (Mb, GRCh37)', ylab = '-log10(P)', main = paste('Pancreas', z))
      graphics::abline(h = -log10(5e-8), lty = 2, col = 'grey50')
      graphics::abline(v = center/1e6, lty = 3, col = '#A54B3C')
    }
  }
  plot_region()
  cat('\n===== Validation summary =====\n')
  print(as.data.frame(summary), row.names = FALSE, digits = 6)
  cat('\nsdY=1 is a labelled sensitivity assumption, not a confirmed data property.\n')
  cat('Warnings are retained in warnings.txt. Sample overlap and multiple signals remain unresolved.\n')
  cat('Original data and earlier outputs were not modified.\nOutput:', dest, '\n')
  writeLines(capture.output(sessionInfo()), file.path(dest, 'sessionInfo.txt'))
})()
