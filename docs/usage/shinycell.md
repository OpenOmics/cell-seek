# <code>cell-seek <b>shinycell</b></code>

## 1. About
The `cell-seek` executable is composed of several inter-related sub commands. Please see `cell-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>cell-seek <b>shinycell</b></code> sub command in more detail. With minimal configuration, the **`shinycell`** sub command enables you to build a ShinyCell2 application from a Seurat RDS object.

The ShinyCell2 application creates an interactive web application for exploring single-cell data. The command validates inputs, creates configuration files, and writes the files and/or Shiny code required to deploy an interactive single-cell web app.

You can find more information about ShinyCell2 at its [GitHub repository](https://github.com/the-ouyang-lab/ShinyCell2) and the OpenOmics fork of ShinyCell2 [here](https://github.com/OpenOmics/ShinyCell2).
Setting up the ShinyCell2 pipeline is fast and easy! In its most basic form, <code>cell-seek <b>shinycell</b></code> only has *3 required inputs*:

- Seurat RDS object file. (with a set of valid cluster identities labels) 
- A project title.
- An output directory to write the shiny application files.

A user can optionally provide a TSV file with marker genes to add a DEG page to the application.

!!! important "Supported Assays"
    `ShinyCell2` does not support all `Seurat` object assays. **Mileage may vary with the Seurat Object you use.**
    **_HTO assays_** and **_Azimuth cell.annotation_** `(prediction.score.celltype.l*)` assays will be removed from your object when using the ShinyCell2 pipeline.
    
    These assay types are unsupported currently by the original ShinyCell2.

    For more information about support for assays and extended format support please consult the original ShinyCell2 repository: https://github.com/the-ouyang-lab/ShinyCell2

## 2. Synopsis
```text
$ cell-seek shinycell [-h] \
      -j OBJ -o OUTDIR --proj PROJECT \
      [--markers MARKERS] [--cluster_labels CLUSTER_LABELS] \
      [--rmmeta META1,META2] [--defred DEFAULTRED] \
      [-l MAXLEVELS] [-a ASSAYS] \
      [--mode {slurm,local}] [--threads THREADS] \
      [--tmp-dir TMP_DIR] [--job-name JOB_NAME] [--silent]
```

The synopsis for this command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a path to a Seurat RDS object via `-j` or `--obj` argument, an output directory to store results via `-o` or `--outdir` argument, and a project name via `--proj` argument.

You can always use the `-h` option for information on a specific command.

### 2.1 Required Arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `-j OBJ, --obj OBJ`  
> **Path to Seurat RDS object file.**  
> *type: file*  
>
> Path to an RDS file produced with saveRDS() that contains a Seurat object. The script will attempt to read and validate the Seurat object.
>
> ***Example:*** `-j seurat_object.rds`

---  
  `-o OUTDIR, --outdir OUTDIR`
> **Path to an output directory.**   
> *type: path*
>   
> Output directory where the Shiny app files and generated code will be written. Directory will be created if it does not exist.
>
> ***Example:*** `-o /srv/shiny-server/my_app`

---  
  `--proj PROJECT`
> **Project name for the application.**   
> *type: string*
>   
> Project name used as the app title and prefix for generated assets.
>
> ***Example:*** `--proj NCBR-34`

### 2.2 Analysis Options

Each of the following arguments are optional, and do not need to be provided.

  `--markers MARKERS`
> **Path to a marker genes file.**   
> *type: file*
>   
> Path to a plain text file with marker genes, one gene per line. When provided, requires `--cluster_labels` to be set as well. Used to generate marker-based DEG pages or panels. If not provided, no DEG pages/panels will be created.
>
> ***Example:*** `--markers markers.txt`

---  
  `--cluster_labels CLUSTER_LABELS`
> **Metadata column name for cluster identities.**   
> *type: string*
>   
> Metadata column name in seurat_obj@meta.data that contains cluster identities corresponding to the markers file. Required when `--markers` is provided. If not provided along with markers, no DEG pages/panels will be created.
>
> ***Example:*** `--cluster_labels seurat_clusters`

---  
  `--rmmeta META1,META2`
> **Metadata columns to remove.**   
> *type: string*
>   
> Comma-delimited list of metadata column names to drop from the Seurat object before deploying. Useful to remove large or sensitive columns. If not provided, no metadata columns are removed.
>
> ***Example:*** `--rmmeta barcode,unwanted_col`

---  
  `--defred DEFAULTRED`
> **Default dimensionality reduction.**   
> *type: string*
>   
> Default dimensionality reduction to use for the app (must exist in Seurat object reductions). The code expects the reduction to expose a key so UMAP1/UMAP2 (or <KEY>1/<KEY>2) coordinate names can be derived. If not provided, uses all active reductions.
>
> ***Example:*** `--defred UMAP`

---  
  `-l MAXLEVELS, --maxlevels MAXLEVELS`
> **Maximum levels for categorical metadata.**   
> *type: int*  
> *default: 50*
>   
> Maximum allowed number of factor levels for categorical metadata fields. Fields exceeding this will be treated as continuous / filtered out.
>
> ***Example:*** `-l 50`

---  
  `-a ASSAYS, --assay ASSAYS`
> **Assays to include in the app.**   
> *type: string*
>   
> Comma-delimited list of assays to include in the Shiny app (e.g. RNA,spatial,ATAC). If omitted, default assays present in the Seurat object are used when compatible.
>
> ***Example:*** `-a RNA,spatial`

### 2.3 Orchestration Options

Each of the following arguments are optional, and do not need to be provided.

  `--mode {slurm,local}`  
> **Execution Method.**  
> *type: string*  
> *default: slurm*
>
> Execution Method. Defines the mode or method of execution. Valid mode options include: slurm or local.
>
> ***slurm***    
> The slurm execution method will submit jobs to the [SLURM workload manager](https://slurm.schedmd.com/). It is recommended running cell-seek in this mode as execution will be significantly faster in a distributed environment. This is the default mode of execution.
>
> ***local***  
> Local executions will run serially on compute instance. This is useful for testing, debugging, or when a users does not have access to a high performance computing environment. If this option is not provided, it will default to a local execution mode.
>
> ***Example:*** `--mode slurm`

---  
  `--job-name JOB_NAME`  
> **Set the name of the pipeline's master job.**  
> *type: string*
> *default: pl:cell-seek:shinycell*
>
> When submitting the pipeline to a job scheduler, like SLURM, this option allows you to set the name of the pipeline's master job. By default, the name of the pipeline's master job is set to "pl:cell-seek:shinycell".
>
> ***Example:*** `--job-name pl_shinycell_id-42`

---  
  `--silent`            
> **Silence standard output.**  
> *type: boolean flag*
>
> Reduces the amount of information directed to standard output when submitting master job to the job scheduler. Only the job id of the master job is returned.
>
> ***Example:*** `--silent`

---  
  `--threads THREADS`   
> **Max number of threads for each process.**  
> *type: int*  
> *default: 2*
>
> Max number of threads for each process. This option is more applicable when running the pipeline with `--mode local`. It is recommended setting this value to the maximum number of CPUs available on the host machine.
>
> ***Example:*** `--threads 12`

---  
  `--tmp-dir TMP_DIR`   
> **Path for temporary files.**  
> *type: path*  
> *default: `/lscratch/$SLURM_JOBID`*
>
> Path on the file system for writing temporary output files. By default, the temporary directory is set to '/lscratch/$SLURM_JOBID' for backwards compatibility with the NIH's Biowulf cluster; however, if you are running the pipeline on another cluster, this option will need to be specified. Ideally, this path should point to a dedicated location on the filesystem for writing tmp files. On many systems, this location is set to somewhere in /scratch. If you need to inject a variable into this string that should NOT be expanded, please quote this options value in single quotes.
>
> ***Example:*** `--tmp-dir /scratch/$USER/`

### 2.4 Miscellaneous Options  
Each of the following arguments are optional, and do not need to be provided.

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
>
> Shows command's synopsis, help message, and an example command
>
> ***Example:*** `--help`

## 3. Example

### 3.1 Basic ShinyCell2 Application

Create a basic ShinyCell2 app from a Seurat RDS file:

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2.) Run the shinycell subcommand
./cell-seek shinycell -j seurat_obj.rds \
                       -o ./ShinyCell2_app \
                       --proj NCBR-34 \
                       --mode slurm
```

### 3.2 Including Marker Genes

Create a ShinyCell2 app with marker genes and cluster labels:

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2.) Run the shinycell subcommand with markers
./cell-seek shinycell -j seurat_obj.rds \
                       -o ./ShinyCell2_markers \
                       --proj NCBR-34 \
                       --markers markers.txt \
                       --cluster_labels seurat_clusters \
                       --mode slurm
```

### 3.3 Advanced Configuration

Create a ShinyCell2 app with specific assays and metadata filtering:

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2.) Run the shinycell subcommand with advanced options
./cell-seek shinycell -j seurat_obj.rds \
                       -o ./ShinyCell2_advanced \
                       --proj NCBR-34 \
                       -a RNA,spatial \
                       --rmmeta large_matrix,sensitive_col \
                       --defred UMAP \
                       -l 100 \
                       --mode slurm
```

### 3.4 Local Execution Mode

Run ShinyCell2 in local mode for testing or debugging:

```bash
# Step 1.) Ensure you have sufficient resources
# on your current node

# Step 2.) Run the shinycell subcommand in local mode
./cell-seek shinycell -j seurat_obj.rds \
                       -o ./ShinyCell2_local \
                       --proj NCBR-34 \
                       --mode local \
                       --threads 8
```

## 4. Marker File Format

The `--markers` option expects a TSV (tab-separated values) file in the same format as produced by Seurat's `FindAllMarkers()` function. This file contains differential expression results for marker genes across clusters.

### 4.1 Required Columns

The marker file must contain the following columns:

- **cluster**: The cluster identity (must match values in the metadata column specified by `--cluster_labels`)
- **gene**: Gene name or identifier
- **p_val**: Unadjusted p-value
- **avg_log2FC**: Average log2 fold-change
- **pct.1**: Percentage of cells expressing the gene in the cluster
- **pct.2**: Percentage of cells expressing the gene in all other clusters
- **p_val_adj**: Adjusted p-value (e.g., Bonferroni correction)

### 4.2 Example Format

```tsv
p_val	avg_log2FC	pct.1	pct.2	p_val_adj	cluster	gene
0	2.5	0.95	0.25	0	0	CD3D
0	2.3	0.92	0.18	0	0	CD3E
1.2e-250	1.8	0.85	0.35	2.4e-246	1	CD14
```

### 4.3 FindAllMarkers vs FindMarkers

Seurat provides two main functions for identifying marker genes:

#### FindAllMarkers (Recommended)

`FindAllMarkers()` automatically identifies markers for **all clusters** and includes the `cluster` column in its output, making it directly compatible with ShinyCell2:

```r
# Run FindAllMarkers - produces ShinyCell2-compatible output
markers <- FindAllMarkers(seurat_obj, 
                         only.pos = TRUE,
                         min.pct = 0.25,
                         logfc.threshold = 0.25)

# Save directly - cluster column is already included
write.table(markers, 
           file = "markers.tsv",
           sep = "\t",
           quote = FALSE,
           row.names = FALSE)
```

#### FindMarkers (Requires Post-Processing)

`FindMarkers()` compares gene expression between **two specific groups** but does **not** include cluster labels in its output. To use `FindMarkers()` results with ShinyCell2, you must manually add the `cluster` column.

The following script demonstrates how to run `FindMarkers()` on each cluster and create `FindAllMarkers`-like output:

```r
#!/usr/bin/env Rscript
# generate_markers_from_findmarkers.R
# Creates FindAllMarkers-compatible output using FindMarkers for each cluster

library(Seurat)

# Load your Seurat object
seurat_obj <- readRDS("seurat_object.rds")

# Specify the metadata column containing cluster identities
cluster_column <- "seurat_clusters"  # Change this to match your metadata column

# Get unique cluster identities
clusters <- unique(seurat_obj@meta.data[[cluster_column]])

# Initialize an empty list to store results
all_markers <- list()

# Loop through each cluster and run FindMarkers
for (cluster_id in clusters) {
  cat(paste0("Processing cluster: ", cluster_id, "\n"))
  
  # Run FindMarkers comparing this cluster vs all others
  tryCatch({
    markers <- FindMarkers(
      seurat_obj,
      ident.1 = cluster_id,
      ident.2 = NULL,  # Compare against all other cells
      group.by = cluster_column,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      test.use = "wilcox"  # Wilcoxon rank sum test (default)
    )
    
    # Add gene names as a column (from rownames)
    markers$gene <- rownames(markers)
    
    # Add cluster identity column
    markers$cluster <- cluster_id
    
    # Store in list
    all_markers[[as.character(cluster_id)]] <- markers
    
  }, error = function(e) {
    cat(paste0("Error processing cluster ", cluster_id, ": ", e$message, "\n"))
  })
}

# Combine all results into a single data frame
combined_markers <- do.call(rbind, all_markers)

# Reset row names
rownames(combined_markers) <- NULL

# Reorder columns to match FindAllMarkers output
combined_markers <- combined_markers[, c("p_val", "avg_log2FC", "pct.1", "pct.2", 
                                         "p_val_adj", "cluster", "gene")]

# Optional: Filter by significance and fold-change
combined_markers <- combined_markers[combined_markers$p_val_adj < 0.05 & 
                                     combined_markers$avg_log2FC > 0.25, ]

# Optional: Sort by cluster and adjusted p-value
combined_markers <- combined_markers[order(combined_markers$cluster, 
                                           combined_markers$p_val_adj), ]

# Save to TSV file
write.table(combined_markers,
           file = "markers_from_findmarkers.tsv",
           sep = "\t",
           quote = FALSE,
           row.names = FALSE)

cat(paste0("\nGenerated ", nrow(combined_markers), " marker entries across ", 
          length(unique(combined_markers$cluster)), " clusters\n"))
cat("Output saved to: markers_from_findmarkers.tsv\n")
```

#### Key Differences Summary

| Feature | FindAllMarkers | FindMarkers |
|---------|---------------|-------------|
| **Cluster column** | ✅ Included automatically | ❌ Must be added manually |
| **Scope** | All clusters at once | One comparison at a time |
| **Use case** | Quick marker discovery | Specific pairwise comparisons |
| **ShinyCell2 compatibility** | ✅ Direct | ⚠️ Requires post-processing |

### 4.4 Usage Example

Once you have generated your marker file using either method:

```bash
./cell-seek shinycell -j seurat_obj.rds \
                       -o ./ShinyCell2_app \
                       --proj MyProject \
                       --markers markers.tsv \
                       --cluster_labels seurat_clusters \
                       --mode slurm
```

## 5. Notes

- When using `--markers` you **MUST** also specify `--cluster_labels`.
- The `--defred` option must correspond to an existing reduction name in the Seurat object; the script will resolve the reduction key to form e.g. UMAP1/UMAP2.
- Use `--rmmeta` to remove columns that are large, private, or not appropriate for public deployment.

