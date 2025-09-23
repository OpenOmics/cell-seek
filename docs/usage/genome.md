# <code>cell-seek <b>genome</b></code>

## 1. About
The `cell-seek` executable is composed of several inter-related sub commands. Please see `cell-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>cell-seek <b>genome</b></code> sub command in more detail. With minimal configuration, the **`genome`** sub command enables you to create a custom reference genome JSON file that can be used in the **`run`** sub command.

If a reference genome that does not come with the pipeline, then a custom json file needs to be provided to run.

This command does not help with creating the 10x compatible reference itself, that would need to be done separately. 10x documentation about the process can be found for [GEX](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-3p-references), [VDJ](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-5p-references), [ATAC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/analysis/inputs/creating-a-reference-package-mkref), and [Multiome](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/inputs/mkref)


Creating a custom reference genome file is fast and easy! In its most basic form, <code>cell-seek <b>genome</b></code> only has *one required input* with the optional arguments supplying the reference paths.

## 2. Synopsis
```text
$ ./cell-seek genome [-h] \
           [--gex GEX_PATH] [--vdj VDJ_PATH] \
           [--atac ATAC_PATH] [--multiome MULTIOME_PATH] \
           --name NAME
```

The synopsis for this command shows its parameters and their usage. Optional parameters are shown in square brackets.

A user **must** provide the name of the reference via `--name` argument. After running the genome sub command, you can use the generated json file as the reference when running the pipeline.

Use you can always use the `-h` option for information on a specific command.

### 2.1 Required Arguments

  `--name NAME`
> **Name of the reference.**
> *type: string*
>
> The name of the reference, which will also be used as the name of the generated reference file.
> ***Example:*** `--name NAME`

### 2.2 Reference Options

Each of the following arguments are optional and do not need to be provided. At least one needs to be used for a reference path to be included in the generated file.

`--gex GEX_PATH`
> **Path to the gene expression reference.**
> *type: path*
>
> This reference is the same one used for the gene expression capture for the GEX, CITE, and MULTI pipelines. This flag is the only one needed to supply the value for all three pipelines.
>
> ***Example:*** `--gex /data/GEX_PATH`

---  
`--vdj VDJ_PATH`
> **Path to the vdj reference.**
> *type: path*
>
> This reference is the one used for the VDJ capture for the VDJ and MULTI pipelines. This flag is the only one needed to supply the value for the two pipelines.
>
> ***Example:*** `--vdj /data/VDJ_PATH`

---  
`--atac ATAC_PATH`
> **Path to the atac reference.**
> *type: path*
>
> This reference is used in the ATAC pipeline.
>
> ***Example:*** `--atac /data/ATAC_PATH`

---  
`--multiome MULTIOME_PATH`
> **Path to the multiome reference.**
> *type: path*
>
> This reference is used in the multiome pipeline.
>
> ***Example:*** `--multiome /data/MULTIOME_PATH`

### 2.3 Miscellaneous Options  
Each of the following arguments are optional, and do not need to be provided.

  `-h, --help`
> **Display Help.**
> *type: boolean*
>
> Shows command's synopsis, help message, and an example command
>
> ***Example:*** `--help`

## 3. Example
```bash
# Step 0.) Grab an interactive node (do not run on head node)
srun -N 1 -n 1 --time=12:00:00 -p interactive --mem=8gb  --cpus-per-task=4 --pty bash
module purge
module load singularity snakemake

# Step 1.) Create a custom NAME.json file with provided reference(s)
cell-seek genome --name NAME --gex /data/GEX_PATH
```
