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
if (is.null(args$meta.to.rm) || is.na(args$meta.to.rm) || args$meta.to.rm == "" || args$meta.to.rm == "NA" || tolower(args$meta.to.rm) == "none") {
  rmmeta <- NULL
} else {
  if (grepl(",", args$meta.to.rm, fixed = TRUE)) {
    rmmeta <- unlist(strsplit(args$meta.to.rm, ",", fixed = TRUE))
  } else {
    rmmeta <- c(trimws(gsub("[\r\n]", "", args$meta.to.rm)))
  }
}
if (is.null(args$assaytouse) || is.na(args$assaytouse) || args$assaytouse == "" || args$assaytouse == "NA" || tolower(args$assaytouse) == "none") {
  assaytouse <- NULL
} else {
  if (grepl(",", args$assaytouse, fixed = TRUE)) {
    assaytouse <- unlist(strsplit(args$assaytouse, ",", fixed = TRUE))
  } else {
    assaytouse <- c(trimws(gsub("[\r\n]", "", args$assaytouse)))
  }
}


if (is.null(args$defaultreduction) || is.na(args$defaultreduction) || args$defaultreduction == "" || args$defaultreduction == "NA" || tolower(args$defaultreduction) == "none") {
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
    cat(paste(assay, " unsupported assay removed!"))
    seurat_obj[[assay]] <- NULL
  }
}

# If --assay was provided, validate the requested assays exist, drop all others,
# and promote the first listed assay to the default.  This must happen BEFORE
# createConfig so the config only reflects the retained assays.
if (!is.null(assaytouse)) {
  missing_assays <- setdiff(assaytouse, names(seurat_obj@assays))
  if (length(missing_assays) > 0) {
    fatal(paste0(
      "Pre-flight ERROR: requested assay(s) not found in seurat object: ",
      paste(missing_assays, collapse = ", "), "\n",
      " └── Available assays: ",
      paste(names(seurat_obj@assays), collapse = ", ")
    ))
  }
  assays_to_drop <- setdiff(names(seurat_obj@assays), assaytouse)
  for (a in assays_to_drop) {
    seurat_obj[[a]] <- NULL
    cat(paste0("INFO - Assay '", a, "' dropped (not listed in --assay).\n"))
  }
  SeuratObject::DefaultAssay(seurat_obj) <- assaytouse[1]
  cat(paste0("INFO - Default assay set to '", assaytouse[1], "'.\n"))
}

# Pre-flight check: every assay in the object must cover the same set of cells.
# A mismatch would corrupt the ShinyCell2 data files, which assume cell ordering
# is aligned across assays.
active_assays <- names(seurat_obj@assays)
assay_ncells  <- vapply(active_assays,
                        function(a) ncol(seurat_obj@assays[[a]]),
                        numeric(1))
if (length(unique(assay_ncells)) > 1) {
  cell_report <- paste(
    paste0("    '", names(assay_ncells), "': ", assay_ncells, " cells"),
    collapse = "\n"
  )
  if (!is.null(assaytouse)) {
    fatal(paste0(
      "Pre-flight ERROR: assays have inconsistent cell counts:\n",
      cell_report, "\n",
      " └── All assays must contain the same set of cells.\n",
      " └── Fix: re-integrate or subset assays to a common cell set, or use\n",
      "     --assay to restrict to a single consistent assay (e.g. --assay RNA)."
    ))
  } else {
    cat(paste0(
      "WARN - Pre-flight: assays have inconsistent cell counts:\n",
      cell_report, "\n",
      " └── Continuing with all assays since --assay was not specified.\n"
    ))
  }
} else {
  cat(paste0(
    "INFO - Pre-flight: ", length(active_assays), " assay(s) [",
    paste(active_assays, collapse = ", "), "] each have ",
    unique(assay_ncells), " cells. OK.\n"
  ))
}

# Pre-flight check: reconcile the 'sample' and 'orig.ident' metadata columns.
# These are conventionally redundant per-cell sample identifiers, and downstream
# ShinyCell2 building / the Shiny app expect a 'sample' column to be present.
# Reconcile the two so both are available:
#   - both present      -> leave as-is
#   - only 'orig.ident' -> copy it to 'sample'
#   - only 'sample'     -> copy it to 'orig.ident'
#   - neither present   -> fatal error
# This must run before metadata columns are collected for createConfig so that
# any copied column is included in the app.
if (class(seurat_obj)[1] == "Seurat") {
  has_sample <- "sample" %in% colnames(seurat_obj@meta.data)
  has_orig   <- "orig.ident" %in% colnames(seurat_obj@meta.data)
  if (has_sample && has_orig) {
    cat("INFO - Pre-flight: both 'sample' and 'orig.ident' metadata columns present. OK.\n")
  } else if (has_orig && !has_sample) {
    seurat_obj@meta.data[["sample"]] <- seurat_obj@meta.data[["orig.ident"]]
    cat("INFO - Pre-flight: 'sample' column missing; copied 'orig.ident' to 'sample'.\n")
  } else if (has_sample && !has_orig) {
    seurat_obj@meta.data[["orig.ident"]] <- seurat_obj@meta.data[["sample"]]
    cat("INFO - Pre-flight: 'orig.ident' column missing; copied 'sample' to 'orig.ident'.\n")
  } else {
    fatal(paste0(
      "Pre-flight ERROR: Seurat object is missing both the 'sample' and ",
      "'orig.ident' metadata columns.\n",
      " └── At least one of 'sample' or 'orig.ident' must be present in\n",
      "     seurat_obj@meta.data to identify each cell's sample of origin.\n",
      " └── Available metadata columns: ",
      paste(colnames(seurat_obj@meta.data), collapse = ", ")
    ))
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
# Only remove assay-specific meta columns (nCount_*, nFeature_*), not clustering columns
all_meta_cols <- colnames(seurat_obj@meta.data)
for (meta_col in all_meta_cols) {
  for (assay in unsupported_assays) {
    pattern <- paste0("^(nCount_|nFeature_)", assay)
    if (grepl(pattern, meta_col)) {
      remove_metas <- c(remove_metas, meta_col)
    }
  }
}

# Remove duplicates from remove_metas
remove_metas <- unique(remove_metas)

# Never remove the cluster_labels column if specified
if (!is.null(args$cluster_labels)) {
  remove_metas <- setdiff(remove_metas, args$cluster_labels)
}

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

# Pre-flight check: detect cell-count mismatch between @meta.data and reductions.
# The default reduction is treated as required: any missing cells there are a
# fatal error because the app cannot be rendered without a complete default view.
# Secondary reductions are non-fatal (cells will be dropped with a warning).
if (class(seurat_obj)[1] == "Seurat") {
  # Mirror the same reduction-selection logic used inside makeShinyFilesGEX.
  # check_reduc[1] is the default reduction; the rest are secondary.
  check_reduc <- if (!is.null(files_params$dimred.to.use)) {
    files_params$dimred.to.use
  } else {
    dr_all <- setdiff(names(seurat_obj@reductions), c("pca", "ref.pca", "lsi"))
    c(SeuratObject::DefaultDimReduc(seurat_obj),
      setdiff(dr_all, SeuratObject::DefaultDimReduc(seurat_obj)))
  }
  default_reduc <- check_reduc[1]
  meta_cells    <- rownames(seurat_obj@meta.data)

  # --- Validate the default reduction first ---
  if (!default_reduc %in% names(seurat_obj@reductions)) {
    fatal(paste0(
      "Pre-flight ERROR: default reduction '", default_reduc,
      "' not found in seurat object.\n",
      " └── Available reductions: ",
      paste(names(seurat_obj@reductions), collapse = ", ")
    ))
  }
  default_emb   <- seurat_obj@reductions[[default_reduc]]@cell.embeddings
  default_ncols <- ncol(default_emb)
  default_nrows <- nrow(default_emb)

  if (default_ncols < 2) {
    fatal(paste0(
      "Pre-flight ERROR: default reduction '", default_reduc,
      "' has only ", default_ncols, " dimension(s); at least 2 are required."
    ))
  }

  missing_in_default <- setdiff(meta_cells, rownames(default_emb))
  extra_in_default   <- setdiff(rownames(default_emb), meta_cells)

  if (length(missing_in_default) > 0) {
    fatal(paste0(
      "Pre-flight ERROR: ", length(missing_in_default), " of ", length(meta_cells),
      " cell(s) in @meta.data have no embeddings in the default reduction '",
      default_reduc, "'.\n",
      " └── This indicates the Seurat object is inconsistent (e.g. cells were\n",
      "     added to @meta.data after the reduction was computed, or the object\n",
      "     was created via reference-based mapping that excludes reference cells).\n",
      " └── Fix: re-run the reduction on the full object, or subset @meta.data\n",
      "     to only the cells present in the reduction before calling this script."
    ))
  }

  cat(paste0(
    "INFO - Pre-flight: default reduction '", default_reduc,
    "' OK (", default_nrows, " cells x ", default_ncols, " dims).\n"
  ))

  if (length(extra_in_default) > 0) {
    cat(paste0(
      "WARN - Pre-flight: ", length(extra_in_default),
      " cell(s) in '", default_reduc,
      "' embeddings are absent from @meta.data and will be ignored.\n"
    ))
  }

  # --- Check secondary reductions (non-fatal) ---
  for (dr in check_reduc[-1]) {
    if (!dr %in% names(seurat_obj@reductions)) {
      cat(paste0("WARN - Pre-flight: secondary reduction '", dr,
                 "' not found in seurat object and will be skipped.\n"))
      next
    }
    emb_cells     <- rownames(seurat_obj@reductions[[dr]]@cell.embeddings)
    missing_cells <- setdiff(meta_cells, emb_cells)
    extra_cells   <- setdiff(emb_cells, meta_cells)
    if (length(missing_cells) > 0) {
      cat(paste0("WARN - Pre-flight: ", length(missing_cells),
                 " cell(s) in @meta.data have no embeddings in secondary reduction '",
                 dr, "' and will be dropped from that view.\n"))
    }
    if (length(extra_cells) > 0) {
      cat(paste0("WARN - Pre-flight: ", length(extra_cells),
                 " cell(s) in '", dr,
                 "' embeddings are absent from @meta.data and will be ignored.\n"))
    }
    if (length(missing_cells) == 0 && length(extra_cells) == 0) {
      cat(paste0("INFO - Pre-flight: secondary reduction '", dr, "' OK.\n"))
    }
  }
}

# Pre-flight check: spatial objects require a "sample" metadata column so that
# slide display names can be resolved in the Shiny app.  Mirror the same
# detection logic used by makeShinyFiles / makeShinyFilesSpatial.
if (class(seurat_obj)[1] == "Seurat" &&
    .hasSlot(seurat_obj, "images") &&
    length(seurat_obj@images) > 0) {
  cat("INFO - Pre-flight: spatial data detected (",
      length(seurat_obj@images), " slide(s)).\n")
  if (!"sample" %in% colnames(seurat_obj@meta.data)) {
    fatal(paste0(
      "Pre-flight ERROR: spatial Seurat object is missing the required 'sample' ",
      "metadata column.\n",
      " └── seurat_obj@meta.data must contain a 'sample' column that maps each\n",
      "     cell/spot to its sample identifier (e.g. K100066, K24685, ...).\n",
      " └── Available metadata columns: ",
      paste(colnames(seurat_obj@meta.data), collapse = ", ")
    ))
  }
  sample_vals <- seurat_obj@meta.data[["sample"]]
  n_na <- sum(is.na(sample_vals) | sample_vals == "")
  if (n_na == length(sample_vals)) {
    fatal(paste0(
      "Pre-flight ERROR: the 'sample' metadata column exists but is entirely ",
      "empty (all NA or blank).\n",
      " └── Every cell/spot must have a non-empty sample identifier."
    ))
  }
  if (n_na > 0) {
    cat(paste0(
      "WARN - Pre-flight: ", n_na, " of ", length(sample_vals),
      " cell(s) have NA or blank 'sample' values.\n"
    ))
  } else {
    cat(paste0(
      "INFO - Pre-flight: 'sample' metadata column OK (",
      length(unique(sample_vals)), " unique sample(s)).\n"
    ))
  }
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
