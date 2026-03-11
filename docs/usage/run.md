# <code>cell-seek <b>run</b></code>

## 1. About
The `cell-seek` executable is composed of several inter-related sub commands. Please see `cell-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>cell-seek <b>run</b></code> sub command in more detail. The following page lists the options applicable to each of the different pipelines available to select via the <code><b>--pipeline</b></code> flag. With minimal configuration, the **`run`** sub command enables you to start running cell-seek pipeline.

Setting up the cell-seek pipeline is fast and easy! In its most basic form, <code>cell-seek <b>run</b></code> only has *five required inputs*.

## 2. Synopsis
```text
$ cell-seek run [--help] \
      [--dry-run] [--job-name JOB_NAME] [--mode {{slurm,local}}] \
      [--sif-cache SIF_CACHE] [--singularity-cache SINGULARITY_CACHE] \
      [--silent] [--threads THREADS] [--tmp-dir TMP_DIR] \
      [--aggregate {{mapped, none}}] [--exclude-introns] \
      [--chain {{auto, TR, IG}}] \
      [--libraries LIBRARIES] [--features FEATURES] \
      [--cmo-sample CMOSAMPLE] [--cmo-reference CMOREFERENCE] \
      [--hto-sample HTOSAMPLE] \
      [--ocm-sample OCMSAMPLE] \
      [--probe-set PROBESET] [--probe-sample PROBESAMPLE] \
      [--filter FILTER] [--metadata METADATA] [--create-bam] \
      [--rename RENAME] \
      --input INPUT [INPUT ...] \
      --output OUTPUT \
      --pipeline {gex, ...} \
      --genome {hg38, ...} \
      --cellranger {8.0.0, ...}
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of FastQ (globbing is supported) to analyze via `--input` argument, an output directory to store results via `--output` argument, the version of the pipeline to run via `--pipeline` argument, the reference genome to use via `--genome` argument, and the version of cellranger to use via `--cellranger` argument.

You can always use the `-h` option for information on a specific command.

### 2.1 Pipelines 

The following is an overview of the available pipelines. A detailed breakdown of the required and optional arguments for each of the pipelines can be found within their corresponding pages.

 * [<code><b>GEX</b></code>: Run the cell-seek pipeline for GEX (scRNA) only data.](gex.md)
 * [<code><b>VDJ</b></code>: Run the cell-seek pipeline for VDJ only data.](vdj.md)
 * [<code><b>CITE</b></code>: Run the cell-seek pipeline for CITE-seq data. This can handle GEX + CITE or CITE-seq only data.](cite.md)
 * [<code><b>MULTI</b></code>: Run the cell-seek pipeline using Multi. This can handle a range of data combinations. Use this whenever cell multiplexing is performed using CellRanger (HTO, OCM, Flex), or GEX with VDJ data.](multi.md)
 * [<code><b>ATAC</b></code>: Run the cell-seek pipeline for ATAC only data.](atac.md)
 * [<code><b>Multiome</b></code>: Run the cell-seek pipeline for multiomic (ATAC + GEX) data.](multiome.md)


