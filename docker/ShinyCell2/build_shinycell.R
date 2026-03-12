#!/usr/bin/env Rscript
prompt_for_install <- function(pkg) {
  cat(paste0(pkg, " is not installed. Would you like to install it? (y/n) "))
  response <- tolower(readLines("stdin", n = 1))
  if (response == "y") {
    if (pkg == "ShinyCell2") {
      remotes::install_github("OpenOmics/ShinyCell2", quiet = TRUE)
    } else {
      install.packages(pkg, quiet = TRUE)
    }
  }
}

if (!suppressPackageStartupMessages(require("remotes", quietly = TRUE))) prompt_for_install("remotes")
if (!suppressPackageStartupMessages(require("ShinyCell2", quietly = TRUE))) prompt_for_install("ShinyCell2")
if (!suppressPackageStartupMessages(require("Seurat", quietly = TRUE))) prompt_for_install("Seurat")
if (!suppressPackageStartupMessages(require("argparse", quietly = TRUE))) prompt_for_install("argparse")

# give good tracebacks on errors
options(error = function() traceback(2))
err <- function(...) {
  cat(sprintf(...), sep = "\n", file = stderr())
}
fatal <- function(...) {
  err(...)
  quit(status = 1)
}

# Command line arguments
parser <- ArgumentParser()
parser$add_argument(
  "-j", "--obj",
  dest = "object",
  help = "RDS file created with saveRDS containing a seurat object",
  required = TRUE
)
parser$add_argument(
  "-o",
  "--outdir",
  dest = "outdir",
  help = "Output directory for shiny app",
  required = TRUE
)
parser$add_argument(
  "--proj",
  dest = "project",
  help = "Project name, becomes title of the app, example: NCBR-34",
  required = TRUE
)
parser$add_argument(
  "--markers",
  dest = "markers",
  help = "Path to a plain text file with marker genes (one per line).",
  default = NULL
)
parser$add_argument(
  "--cluster_labels",
  dest = "cluster_labels",
  help = "Name of metadata column in seurat_obj@meta.data used as cluster labels",
  default = NULL
)
parser$add_argument(
  "--rmmeta",
  dest = "meta.to.rm",
  help = "Comma delimited list of names in seurat_object@meta.data that will be dropped prior to deploying to ShinyCell2",
  default = NULL,
  metavar = "META.VAR.TO.RM, META.VAR.TO.RM, META.VAR.TO.RM, ..."
)
parser$add_argument(
  "--silent",
  dest = "silent",
  action = "store_true",
  help = "Toggle to suppress verbose output.",
  default = FALSE
)
parser$add_argument("--defred",
  dest = "defaultreduction",
  help = paste0("The default reduction to use with ShinyCell2, this value must exist in Seurat::DefaultDimReduc(obj).\n
    The two major principal components for this reductions must be labeled the same label with a 1 and a 2 trailing it.\n
    i.e. defaultreduction = UMAP, UMAP1 and UMAP2 are the two components that must exist"), default = NULL
)
parser$add_argument(
  "-l", "--maxlevels",
  dest = "max.levels",
  type = "integer",
  help = "The maximum allowable amount levels/factors of categorical values (int).",
  default = NULL
)
parser$add_argument(
  "--filesonly",
  dest = "files.only",
  action = "store_true",
  help = "Toggle for only running makeShinyFiles() without makeShinyCodes().",
  default = FALSE
)
parser$add_argument(
  "--codesonly",
  dest = "codes.only",
  action = "store_true",
  help = "Toggle for only running makeShinyCodes() without makeShinyFiles().",
  default = FALSE
)
parser$add_argument(
  "-a",
  "--assay",
  dest = "assaytouse",
  metavar = "<ASSAY1>, <ASSAY2>, <ASSAY3>, ...",
  help = "The assay to utilize for ShinyCell2 web application. \n
          Comma delimit multiple assays in a single string e.g.: RNA,spatial,ATAC,etc.",
  default = NULL
)
args <- parser$parse_args()

if (args$silent) {
  options(error = function() traceback(1))
  options(warn = -1)
}

# setup opt parse variables for downstream usage into shinycell2
rds_file <- args$object
project_name <- args$project
shiny_app_dir <- normalizePath(args$outdir, mustWork = FALSE)

seurat_obj <- readRDS(rds_file)

# Validate mutually inclusive args: cluster_labels and markers
if (!is.null(args$cluster_labels) && is.null(args$markers)) {
  fatal("--cluster_labels requires --markers to be set as well")
}
if (!is.null(args$markers) && is.null(args$cluster_labels)) {
  fatal("--markers requires --cluster_labels to be set as well")
}

# If markers file provided, check it exists and is readable
markers_list <- NULL
if (!is.null(args$markers)) {
  if (!file.exists(args$markers)) {
    fatal(paste0("Markers file not found: ", args$markers))
  }
  if (file.access(args$markers, mode = 4) != 0) {
    fatal(paste0("Markers file is not readable: ", args$markers))
  }
  # read markers as plain text, one per line
  markers_list <- tryCatch(
    {
      scan(args$markers, what = character(), sep = "\n", quiet = TRUE)
    },
    error = function(e) {
      fatal(paste0("Failed to read markers file: ", e$message))
    }
  )
}
# Validate cluster_labels exists in seurat object metadata
if (!is.null(args$cluster_labels)) {
  if (!(args$cluster_labels %in% colnames(seurat_obj@meta.data))) {
    fatal(paste0("cluster_labels '", args$cluster_labels, "' not found in seurat_obj@meta.data"))
  }
}
if (is.null(args$meta.to.rm) || is.na(args$meta.to.rm) || args$meta.to.rm == "" || args$meta.to.rm == "NA") {
  rmmeta <- NULL
} else {
  if ("," %in% args$meta.to.rm) {
    rmmeta <- unlist(strsplit(args$meta.to.rm, ",", fixed = TRUE))
  } else {
    rmmeta <- c(trimws(gsub("[\r\n]", "", args$meta.to.rm)))
  }
}
if (is.null(args$assaytouse) || is.na(args$assaytouse) || args$assaytouse == "" || args$assaytouse == "NA") {
  assaytouse <- NULL
} else {
  if ("," %in% args$assaytouse) {
    assaytouse <- unlist(strsplit(args$assaytouse, ",", fixed = TRUE))
  } else {
    assaytouse <- c(trimws(gsub("[\r\n]", "", args$assaytouse)))
  }
}
if (is.null(args$defaultreduction) || is.na(args$defaultreduction) || args$defaultreduction == "" || args$defaultreduction == "NA") {
  defaultreduction <- NULL
} else {
  if (args$defaultreduction %in% names(seurat_obj@reductions)) {
    this_key <- seurat_obj@reductions[[args$defaultreduction]]@key
    defaultreduction <- c(paste0(this_key, "1"), paste0(this_key, "2"))
  } else {
    fatal(paste0("`", args$defaultreduction, "` reduction not found in seurat object!"))
  }
}

valid_max_levels <- !is.null(args$max.levels) &&
  length(args$max.levels) == 1 &&
  !is.na(args$max.levels)

if (valid_max_levels) {
  max.levels <- args$max.levels
} else {
  max.levels <- NULL
}
if (!dir.exists(shiny_app_dir)) {
  dir.create(shiny_app_dir, recursive = TRUE, showWarnings = FALSE)
}

# Sanity check: Does the RDS file
# actually contain a seurat object?
if (class(seurat_obj) == "SeuratObject") {
  # It doesn't look like it...
  err("Fatal Error: Failed to provide an RDS file with a Seuart Object.")
  fatal(" └── Please create a new RDS file with a seurat object!")
}

# Remove unsupported assay
# ShinyCell2 supports:
#   - CITEseq
#   - spatial
#   - scATAC (Signac)
# ShinyCell2 does NOT support:
#   - celltype "assays" added by azimuth (prediction.score.celltype.l*)
#   - HTO
#   - sketch assays for large datasets
for (assay in Assays(seurat_obj)) {
  if (startsWith(assay, "prediction.score.celltype.l")) {
    seurat_obj[[assay]] <- NULL
  }
}

unsupported_assays <- c("HTO", "sketch")
for (assay in unsupported_assays) {
  if (assay %in% names(seurat_obj@assays)) {
    cat(paste0(assay, "unsupported assay removed!", sep = " "))
    seurat_obj[[assay]] <- NULL
  }
}

# Create ShinyCell config file
# to make the application
config_params <- list()
if (!is.null(max.levels)) {
  config_params$maxLevels <- max.levels
}

# Determine which metadata to include based on remove_metas list
remove_metas <- c()

if (!is.null(rmmeta)) {
  remove_metas <- c(remove_metas, rmmeta)
}

# Add any metadata columns related to unsupported assays to remove_metas
all_meta_cols <- colnames(seurat_obj@meta.data)
for (meta_col in all_meta_cols) {
  for (assay in unsupported_assays) {
    if (grepl(assay, meta_col, fixed = TRUE)) {
      remove_metas <- c(remove_metas, meta_col)
    }
  }
}

# Remove duplicates from remove_metas
remove_metas <- unique(remove_metas)

# Create meta.to.include: all metadata columns NOT in remove_metas
meta_to_include <- setdiff(all_meta_cols, remove_metas)

# Add meta.to.include to config parameters if there are columns to include
if (length(meta_to_include) > 0) {
  # config_params$meta.to.del <- remove_metas
  config_params$meta.to.include <- meta_to_include
}

shinycell_config <- do.call(
  createConfig,
  c(seurat_obj, config_params)
)

# if (!is.null(length(meta_to_include) > 0)) {
#   shinycell_config$meta.to.del <- remove_metas
#   shinycell_config <- delMeta(shinycell_config, remove_metas)
# }

# Build the Shiny Application,
# in the default location for
# Shiny/Posit server: i.e.
# /srv/shiny-server/${app_name}
files_params <- list(
  seurat_obj,
  shinycell_config,
  shiny.dir = shiny_app_dir,
  shiny.prefix = "sc1"
)

if (!is.null(defaultreduction)) {
  files_params$dimred.to.use <- args$defaultreduction
  files_params$default.dimred <- defaultreduction
}

if (!is.null(args$markers) && !is.null(args$cluster_labels)) {
  files_params$precomputed.deg <- args$markers
  files_params$clusters <- args$cluster_labels
}

if (!args$codes.only) {
  start_time <- Sys.time()
  cat(paste0("INFO - ", start_time, " - start building ShinyCell2 files\n"))
  do.call(
    makeShinyFiles,
    files_params
  )
  end_time <- Sys.time()
  cat(paste0("INFO - ", end_time, " - stop building ShinyCell2 files\n"))
  cat(paste0("INFO - Total time taken to build ShinyCell2 files: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes\n"))
}
cat("\n")
if (!args$files.only) {
  start_time <- Sys.time()
  cat(paste0("INFO - ", start_time, " - start building ShinyCell2 codes\n"))
  makeShinyCodes(
    shiny.title = project_name,
    shiny.dir = shiny_app_dir,
    shiny.prefix = "sc1"
  )
  end_time <- Sys.time()
  cat(paste0("INFO - ", end_time, " - stop building ShinyCell2 codes\n"))
  cat(paste0("INFO - Total time taken to build ShinyCell2 codes: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes\n"))
}
