# Run from the repository root: Rscript .tests/test_qc_threshold_defaults.R
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(dplyr))

# Evaluate the actual QC section without running the rest of the Seurat pipeline.
read_qc <- function(path) {
  code <- parse(path)
  start <- which(vapply(code, function(x) {
    is.call(x) && identical(x[[1]], as.name('<-')) &&
      identical(x[[2]], as.name('defaultThreshold'))
  }, logical(1)))
  end <- which(vapply(code, function(x) {
    is.call(x) && identical(x[[1]], as.name('write.csv')) &&
      identical(x[[2]], as.name('thresh'))
  }, logical(1)))
  stopifnot(length(start) == 1, length(end) == 1, start < end)
  code[start:end]
}
qc_code <- read_qc('workflow/scripts/seuratSampleQC.R')

run_qc <- function(metadata, csv=NULL, sample='sampleA', code=qc_code) {
  workdir <- tempfile('qc-threshold-test-')
  dir.create(workdir)
  previous <- setwd(workdir)
  on.exit({setwd(previous); unlink(workdir, recursive=TRUE)})
  counts <- matrix(1L, nrow=3, ncol=nrow(metadata),
                   dimnames=list(paste0('gene', 1:3), rownames(metadata)))
  seur <- CreateSeuratObject(Matrix::Matrix(counts, sparse=TRUE))
  for (name in names(metadata)) seur[[name]] <- metadata[[name]]
  filterfile <- NA_character_
  if (!is.null(csv)) {
    filterfile <- file.path(workdir, 'filter.csv')
    writeLines(csv, filterfile)
  }
  env <- new.env()
  env$seur <- seur
  env$opt <- list(filterfile=filterfile, sample=sample)
  eval(code, env)
  recorded <- read.csv('cell_filter_info.csv')
  stopifnot(
    isTRUE(all.equal(as.list(recorded), env$thresh, check.attributes=FALSE)),
    env$thresh$numCellsRemove == length(env$cellsToRemove),
    env$thresh$finalCellCount == nrow(metadata) - length(env$cellsToRemove)
  )
  list(thresh=env$thresh, removed=env$cellsToRemove, recorded=recorded)
}

metadata <- data.frame(
  nFeature_RNA=rep(1000, 22), nCount_RNA=rep(5000, 22),
  percent.mito=c(rep(3, 12), 4, 5, 6, 1, 2.5, 2.6, 3, 3, 3, 3),
  percent.rpl=rep(0, 22),
  row.names=c(paste0('ordinary', 1:12), 'rescued', 'mito_equal', 'mito_high',
              'mito_low', 'decimal_equal', 'decimal_high',
              'gene_low', 'gene_high', 'umi_low', 'umi_high')
)
metadata['gene_low', 'nFeature_RNA'] <- 500
metadata['gene_high', 'nFeature_RNA'] <- 1500
metadata['umi_low', 'nCount_RNA'] <- 2500
metadata['umi_high', 'nCount_RNA'] <- 7500
metadata['ordinary1', 'percent.rpl'] <- 10

automatic <- run_qc(metadata)
bounds <- c('nFeature_RNA_low', 'nFeature_RNA_high',
            'nCount_RNA_low', 'nCount_RNA_high', 'percent.mito_high')
stopifnot(
  identical(unlist(automatic$thresh[bounds], use.names=FALSE), c(1000, 1000, 5000, 5000, 3)),
  setequal(automatic$removed, c('gene_low', 'gene_high', 'umi_low', 'umi_high',
                               'rescued', 'mito_equal', 'mito_high'))
)

# Nonzero MADs: frozen expected cutoffs and rejected cells for this fixture.
varied <- data.frame(
  nFeature_RNA=c(rep(c(800, 900, 1000, 1100, 1200), 4), 10, 10000),
  nCount_RNA=c(rep(c(4000, 4500, 5000, 5500, 6000), 4), 50, 50000),
  percent.mito=c(rep(1:5, 4), 0, 90), row.names=paste0('cell', 1:22)
)
varied_result <- run_qc(varied)
stopifnot(
  identical(unlist(varied_result$thresh[bounds], use.names=FALSE), c(626, 1598, 3129, 7989, 13)),
  setequal(varied_result$removed, c('cell21', 'cell22'))
)
cat('PASS: automatic thresholds and filtering, including nonzero MADs\n')

partial_csv <- c('Sample,percent.mito_high', 'sampleA,5')
partial <- run_qc(metadata, partial_csv)
stopifnot(
  partial$thresh$percent.mito_high == 5,
  identical(partial$thresh[bounds[1:4]], automatic$thresh[bounds[1:4]]),
  setequal(partial$removed, c('gene_low', 'gene_high', 'umi_low', 'umi_high', 'mito_high'))
)
cat('PASS: mitochondrial-only override retains gene and UMI filtering and includes the boundary\n')

unlisted <- run_qc(metadata, partial_csv, sample='sampleB')
duplicate <- run_qc(metadata, c(partial_csv, 'sampleA,6'))
stopifnot(identical(unlisted, automatic), identical(duplicate, automatic))
cat('PASS: unlisted samples and existing duplicate-row fallback use defaults\n')

missing <- run_qc(metadata, c('Sample,percent.mito_high,nFeature_RNA_low,nCount_RNA_high',
                             'sampleA,5,,NA', 'sampleB,NA,900,6000'))
whitespace <- run_qc(metadata, c('Sample,percent.mito_high', 'sampleA," "', 'sampleB,5'))
stopifnot(identical(missing, partial), identical(whitespace, automatic))
cat('PASS: blank, NA, whitespace, and omitted values retain defaults\n')

zero <- run_qc(metadata, c('Sample,percent.mito_high', 'sampleA,0'))
stopifnot(zero$thresh$percent.mito_high == 0, length(zero$removed) == nrow(metadata))
strict <- run_qc(metadata, c('Sample,percent.mito_high', 'sampleA,2.5'))
stopifnot(
  strict$thresh$percent.mito_high == 2.5,
  !('decimal_equal' %in% strict$removed), 'decimal_high' %in% strict$removed,
  all(automatic$removed %in% strict$removed),
  !('rescued' %in% partial$removed), 'rescued' %in% automatic$removed
)
cat('PASS: zero, stricter, decimal, and more permissive overrides replace defaults\n')

multiple <- run_qc(metadata, c('Sample,percent.mito_HIGH,nCount_RNA_high', 'sampleA,5,10000'))
stopifnot(
  multiple$thresh$percent.mito_high == 5, multiple$thresh$nCount_RNA_high == 10000,
  !('umi_high' %in% multiple$removed), 'umi_low' %in% multiple$removed,
  identical(multiple$thresh[bounds[c(1, 2, 3)]], automatic$thresh[bounds[c(1, 2, 3)]])
)
lower <- run_qc(metadata, c('Sample,nFeature_RNA_low', 'sampleA,500'))
stopifnot(lower$thresh$nFeature_RNA_low == 500,
          !('gene_low' %in% lower$removed), 'gene_high' %in% lower$removed,
          identical(lower$thresh[bounds[-1]], automatic$thresh[bounds[-1]]))
extra <- run_qc(metadata, c('Sample,percent.rpl_high', 'sampleA,5'))
stopifnot(extra$thresh$percent.rpl_high == 5,
          setequal(extra$removed, union(automatic$removed, 'ordinary1')))
cat('PASS: multiple overrides, case-insensitive suffixes, and additional metadata cutoffs\n')

invalid <- tryCatch(run_qc(metadata, c('Sample,percent.mito_high', 'sampleA,invalid')),
                    error=function(e) conditionMessage(e))
stopifnot(is.character(invalid), grepl('sampleA', invalid),
          grepl('percent.mito_high', invalid), grepl('numeric', invalid))
cat('PASS: nonnumeric overrides fail with sample and threshold names\n')

# Evaluate report chunks against the CSV produced by the production QC code.
read_chunk <- function(path, pattern, section=NULL) {
  lines <- readLines(path)
  if (!is.null(section)) lines <- lines[seq.int(grep(section, lines), length(lines))]
  start <- grep(pattern, lines)
  if (!is.null(section)) start <- head(start, 1)
  stopifnot(length(start) == 1)
  remaining <- lines[(start + 1):length(lines)]
  parse(text=remaining[seq_len(which(remaining == '```')[1] - 1)])
}
report_test <- function() {
  workdir <- tempfile('qc-report-test-')
  dir.create(workdir)
  on.exit(unlink(workdir, recursive=TRUE))
  decimal <- run_qc(metadata, c('Sample,percent.mito_high', 'sampleA,2.555'))
  write.csv(decimal$recorded, file.path(workdir, 'cell_filter_info.csv'), row.names=FALSE)
  env <- new.env()
  env$params <- list(defaultfilter=FALSE)
  env$seuratdir <- workdir
  report_code <- read_chunk('workflow/scripts/seuratSampleQCReport.Rmd', '^```\\{r Cell Filter,')
  output <- capture.output(eval(report_code, env))
  stopifnot(any(grepl('Mitochondrial Percentage Thresholds: (NA, 2.555)', output, fixed=TRUE)),
            any(grepl('retain automatic defaults', output, fixed=TRUE)))

  # Run the production CSV merger for an overridden and an automatic sample.
  for (sample in c('sampleA', 'sampleB')) {
    dir.create(file.path(workdir, sample))
    values <- if (sample == 'sampleA') decimal$recorded else automatic$recorded
    write.csv(values, file.path(workdir, sample, 'cell_filter_info.csv'), row.names=FALSE)
  }
  merge_output <- system2(file.path(R.home('bin'), 'Rscript'),
    c(shQuote(normalizePath('workflow/scripts/cellFilterSummary.R')),
      '--datapath', shQuote(workdir), '--output', shQuote(file.path(workdir, 'summary.csv'))),
    stdout=TRUE, stderr=TRUE)
  stopifnot(is.null(attr(merge_output, 'status')))

  # Exercise the summary's actual numeric conversion, before HTML formatting.
  summary_code <- read_chunk('workflow/scripts/seuratSampleQCSummaryReport.Rmd',
                             '^```\\{r  echo=F ,results=', '^## Cell Filter Thresholds')
  env$params$cellfilter <- file.path(workdir, 'summary.csv')
  eval(summary_code[1:2], env)
  stopifnot(env$data$percent.mito_high[env$data$Sample == 'sampleA'] == 2.555,
            env$data$percent.mito_high[env$data$Sample == 'sampleB'] == 3,
            env$data$numCellsRemove[env$data$Sample == 'sampleA'] == length(decimal$removed))

  # The same thresholds are passed to the actual violin-plot function.
  plot_code <- parse('workflow/scripts/seuratSampleQC.R')
  plot_function <- which(vapply(plot_code, function(x) {
    is.call(x) && identical(x[[1]], as.name('<-')) &&
      identical(x[[2]], as.name('doVlnPlot'))
  }, logical(1)))
  env$VlnPlot <- function(...) ggplot2::ggplot()
  env$geom_boxplot <- ggplot2::geom_boxplot
  env$position_dodge <- ggplot2::position_dodge
  env$NoLegend <- function() ggplot2::theme(legend.position='none')
  env$geom_hline <- ggplot2::geom_hline
  eval(plot_code[plot_function], env)
  plot <- env$doVlnPlot('percent.mito', NULL, decimal$thresh)
  stopifnot(tail(plot$layers, 1)[[1]]$data$yintercept == 2.555)
}
report_test()
cat('PASS: final CSV values, decimal report cutoffs, summary values, and plot annotations agree\n')
