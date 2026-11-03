# Optional per-capture correction. Source this file to test the matrix contract.
read_soupx_annotations <- function(path) {
    features <- read.delim(gzfile(file.path(path, "features.tsv.gz")),
        header = FALSE, colClasses = "character", check.names = FALSE,
        quote = "", comment.char = ""
    )
    barcodes <- readLines(gzfile(file.path(path, "barcodes.tsv.gz")))
    list(features = features, barcodes = barcodes)
}

write_soupx_result <- function(corrected, original, annotations, rho, sample, output) {
    input_umis <- sum(original)
    features <- annotations$features
    DropletUtils::write10xCounts(
        file.path(output, "filtered_feature_bc_matrix"), corrected,
        type = "sparse", version = "3", overwrite = TRUE,
        gene.id = features[[1]], gene.symbol = features[[2]], gene.type = features[[3]],
        barcodes = annotations$barcodes
    )
    corrected_umis <- sum(corrected)
    summary <- data.frame(
        sample = sample, soupx_version = as.character(utils::packageVersion("SoupX")),
        contamination_fraction = rho, input_umis = input_umis, corrected_umis = corrected_umis,
        removed_umis = input_umis - corrected_umis,
        removed_percent = 100 * (input_umis - corrected_umis) / input_umis,
        cells = ncol(corrected), features = nrow(corrected), method = "subtraction",
        round_to_integer = TRUE, seed = 1234
    )
    write.table(summary, file.path(output, "soupx_summary.tsv"),
        sep = "\t",
        quote = FALSE, row.names = FALSE
    )
    writeLines(capture.output(utils::sessionInfo()), file.path(output, "sessionInfo.txt"))
}

run_soupx <- function(outs, output, sample, commit) {
    outs <- normalizePath(outs, mustWork = TRUE)
    message(
        "Sample: ", sample, "\nInput: ", outs,
        "\nSource Git commit: ", commit, "\nSeed: 1234"
    )
    filtered_dir <- file.path(outs, "filtered_feature_bc_matrix")
    filtered <- read_soupx_annotations(filtered_dir)
    channel <- SoupX::load10X(outs, channelName = sample, readArgs = list(gene.column = 1))
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(file.path(output, "soupx_estimation.pdf"), width = 8, height = 6)
    on.exit(grDevices::dev.off(), add = TRUE)
    set.seed(1234)
    channel <- SoupX::autoEstCont(channel, doPlot = TRUE)
    rho <- unique(channel$metaData$rho)
    if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho >= 1) {
        stop("Automatic estimation did not produce a valid global contamination fraction")
    }
    corrected <- SoupX::adjustCounts(channel, method = "subtraction", roundToInt = TRUE)
    corrected <- corrected[filtered$features[[1]], filtered$barcodes, drop = FALSE]
    write_soupx_result(corrected, channel$toc, filtered, rho, sample, output)
}

main <- function() {
    options <- list(
        optparse::make_option("--outs", type = "character"),
        optparse::make_option("--output", type = "character"),
        optparse::make_option("--sample", type = "character"),
        optparse::make_option("--commit", type = "character", default = "unknown")
    )
    opt <- optparse::parse_args(optparse::OptionParser(option_list = options))
    run_soupx(opt$outs, opt$output, opt$sample, opt$commit)
}

if (sys.nframe() == 0L) main()
