#!/usr/bin/env Rscript
# Export CellBender's filtered H5 as a Cell Ranger v3 matrix and keep diagnostics.

convert_cellbender <- function(input, output, summary_path, sample, raw_input, artifacts) {
    sce <- DropletUtils::read10xCounts(input, type = "HDF5", version = "3", row.names = "id")
    counts <- methods::as(SummarizedExperiment::assay(sce, "counts"), "dgCMatrix")
    features <- SummarizedExperiment::rowData(sce)
    barcodes <- as.character(SummarizedExperiment::colData(sce)$Barcode)

    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

    DropletUtils::write10xCounts(
        output, counts,
        type = "sparse", version = "3", overwrite = TRUE,
        gene.id = as.character(features$ID), gene.symbol = as.character(features$Symbol),
        gene.type = as.character(features$Type), barcodes = barcodes
    )

    diagnostics <- file.path(dirname(output), "diagnostics")
    dir.create(diagnostics, recursive = TRUE, showWarnings = FALSE)

    files <- file.path(artifacts, c(
        "corrected.log", "corrected.pdf",
        "corrected_metrics.csv", "corrected_cell_barcodes.csv", "corrected_report.html"
    ))

    files <- files[file.exists(files) & !dir.exists(files)]
    if (!all(file.copy(files, diagnostics, overwrite = TRUE))) {
        stop("Could not keep CellBender diagnostics")
    }

    record <- list(
        state = "matrix_export_complete", sample = sample,
        raw_input = raw_input, filtered_h5 = input, matrix_dir = output,
        genes = nrow(counts), cells = ncol(counts), stored_entries = length(counts@x),
        kept_diagnostics = as.list(basename(files)),
        module = "cellbender/0.3.0", cuda = TRUE, fpr = 0.01, epochs = 150,
        expected_cells = "automatic", total_droplets_included = "automatic",
        r_version = R.version.string, r_home = R.home(),
        DropletUtils = as.character(utils::packageVersion("DropletUtils")),
        Matrix = as.character(utils::packageVersion("Matrix"))
    )
    dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(record, summary_path, auto_unbox = TRUE, pretty = TRUE)
}

main <- function() {
    options <- list(
        optparse::make_option("--input", type = "character"),
        optparse::make_option("--output", type = "character"),
        optparse::make_option("--summary_path", type = "character"),
        optparse::make_option("--artifacts", type = "character",
            help = "Inference directory; keep small diagnostics alongside the matrix"
        ),
        optparse::make_option("--sample", type = "character"),
        optparse::make_option("--raw-input", type = "character", dest = "raw_input")
    )
    parser <- optparse::OptionParser(option_list = options)
    opt <- optparse::parse_args(parser)
    if (any(vapply(
        opt[c("input", "output", "summary_path", "sample", "raw_input", "artifacts")],
        is.null,
        logical(1)
    ))) {
        stop("Required options: --input, --output, --summary_path, --sample, --raw-input, --artifacts")
    }
    convert_cellbender(opt$input, opt$output, opt$summary_path, opt$sample, opt$raw_input, opt$artifacts)
}

if (sys.nframe() == 0L) main()
