# <code>cell-seek <b>run</b> --pipeline cite</code>

## 1. About
The `cell-seek` executable is composed of several inter-related sub commands. Please see `cell-seek -h` for all available options.

This part of the documentation describes options and concepts for [<code>cell-seek <b>run</b></code>](run.md) sub command for the **CITE** pipeline which an be selected via the <code><b>--pipeline</b></code> flag in more detail. With minimal configuration, the **`run`** sub command enables you to start running cell-seek pipeline.

Setting up the cell-seek pipeline is fast and easy! In its most basic form, <code>cell-seek <b>run</b></code> only has *five required inputs*.

### 1.1 Use Case

The CITE pipeline should be used when there is antibody capture (CITE-seq) data that will be analyzed. This CITE-seq data can be processed alongside gene expression data or on its own. If modalities other than gene expression were also captured, then [the MULTI pipeline](multi.md) should be used.

This pipeline can also be used if cell multiplexing was performed, but the demultiplexing will not be performed by Cell Ranger. In that scenario then the CITE pipeline should be used and the cell multiplexing tags treated as a normal antibody capture.

A basic guideline of which pipeline should be used for different modalities can be found in the [synopsis section of the <code><b>run</b></code> documentation.](run.md#21-pipelines)

## 2. Synopsis
```text
$ cell-seek run [--help] \
      [--dry-run] [--job-name JOB_NAME] [--mode {{slurm,local}}] \
      [--sif-cache SIF_CACHE] [--singularity-cache SINGULARITY_CACHE] \
      [--silent] [--threads THREADS] [--tmp-dir TMP_DIR] \
      [--exclude-introns] \
      [--libraries LIBRARIES] [--features FEATURES] \
      [--filter FILTER] [--metadata METADATA] [--create-bam] \
      [--forcecells FORCECELLS] \
      --input INPUT [INPUT ...] \
      --output OUTPUT \
      --pipeline cite \
      --genome {hg38, ...} \
      --cellranger {8.0.0, ...}
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of FastQ (globbing is supported) to analyze via `--input` argument, an output directory to store results via `--output` argument, the version of the pipeline to run via `--pipeline` argument, the reference genome to use via `--genome` argument, and the version of cellranger to use via `--cellranger` argument.

You can always use the `-h` option for information on a specific command.


### 2.1 Required Arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input FastQ file(s) or Cell Ranger folder(s).**  
> *type: file(s) or folder(s)*  
>
> FastQ Input: One or more FastQ files can be provided. The pipeline does NOT support single-end data. From the command-line, each input file should separated by a space. Multiple input FastQ files per sample can be provided. Globbing is supported! This makes selecting FastQ files easy. Input FastQ files should always be gzipp-ed.
>
> ***Example:*** `--input .tests/*_R?_fastq.gz`
>
>
> Cell Ranger Input: Cell Ranger output folders can be provided. It is expected that the outs folder is contained within the Cell Ranger output folders, and keep the normal output folder structure. Globbing is supported!
>
> ***Example:*** `--input .tests/*/`

---  
  `--output OUTPUT`
> **Path to an output directory.**   
> *type: path*
>   
> This location is where the pipeline will create all of its output files, also known as the pipeline's working directory. If the provided output directory does not exist, it will be created automatically.
>
> ***Example:*** `--output /data/$USER/cell-seek_out`

---  
  `--pipeline cite`
> **The pipeline to run.**   
> *type: string*
>   
> This option selects the version of the pipeline to run. The documentation provided is based on choosing the option for CITE-seek.
>
> ***Example:*** `--pipeline cite`

---  
  `--genome {hg38, mm10, hg2024, mm2024, custom.json}`
> **Reference genome.**   
> *type: string*
>   
> This option defines the reference genome of the samples. cell-seek does comes bundled with prebuilt reference files for human and mouse samples, The options hg38 or mm10 would select the 2020 release of the reference. The options hg2024 or mm2024 would select the 2024 release of the reference. More information about the officially released references can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-reference-release-notes).
>
> A custom reference genome can also be provided via a json file. Additional information for creating this json file can be found in [<code>cell-seek <b>genome</b></code>](genome.md).
>
> For prebuilt references please select one of the following options: hg38, mm10, hg2024, mm2024
>
> ***Example:*** `--genome hg38`

---  
  `--cellranger {7.1.0, 7.2.0, 8.0.0, 9.0.0, 10.0.0}`
> **The version of Cell Ranger to run.**   
> *type: string*
>   
> This option specifies which version of Cell Ranger to use when running GEX, VDJ, CITE, or MULTI pipelines. Please select one of the following options: 7.1.0, 7.2.0, 8.0.0, 9.0.0, 10.0.0
>
> ***Example:*** `--cellranger 7.1.0`


### 2.2 Conditionally Required Arguments

The following arguments are only required when FastQ files are used as input. They are not required when Cell Ranger output file is used as input.

`--libraries LIBRARIES`
> **Libraries file.**   
> *type: file*
>   
> A CSV file containing information about each library. It contains each sample's name, flowcell, demultiplexed name, and library type. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-feature-bc-analysis).

> *Here is an example libraries.csv file:*
> ```
> Name,Flowcell,Sample,Type
> IL15_LNs,H7CNNBGXG,IL15_LNs,Gene Expression
> IL15_LNs,H7CT7BGXG,IL15_LNs_ADT,Antibody Capture
> EGR1,H7CNNBGXG,EGR1_GEX,Gene Expression
> EGR1,H7CT7BGXG,EGR1_ADT,Antibody Capture
> ```

> *Where:*

> - *Name:* name of the sample passed to Cell Ranger.  
> - *Flowcell:* A unique identifier in the path that the FASTQ files are located.
> - *Sample:* Name that was used when demultiplexing, this should match the FASTQ files.  
> - *Type:* library type for each sample. List of supported options:  
>        * Gene Expression
>        * CRISPR Guide Capture
>        * Antibody Capture
>        * Custom
>
> ***Example:*** `--libraries libraries.csv`


---  
`--features FEATURES`
> **Features file.**   
> *type: file*
>   
> A feature reference CSV file containing information for processing a feature barcode data. This file should contain a unique ID for the feature, a human readable name, sequence, feature type, read, and pattern. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-feature-bc-analysis).

> *Here is an example features.csv file:*
> ```
> id,name,sequence,feature_type,read,pattern
> CITE_CD64,CD64,AGTGGG,Antibody Capture,R2,5PNN(BC)
> CITE_CD8,CD8,TCACCGT,Antibody Capture,R2,5PNNN(BC)
> ```

> *Where:*  

> - *id:* Unique ID for this feature. Must not contain whitespace, quote or comma characters. Each ID must be unique and must not collide with a gene identifier from the transcriptome.
> - *name:* Human-readable name for this feature. Must not contain whitespace.
> - *sequence:* Nucleotide barcode sequence associated with this feature, e.g. the antibody barcode or sgRNA protospacer sequence.
> - *read:* Specifies which RNA sequencing read contains the Feature Barcode sequence. Must be R1 or R2, but in most cases R2 is the correct read.
> - *pattern:* Specifies how to extract the sequence of the feature barcode from the read.

> - *Type:* Type of the feature. List of supported options:  
>        * Antibody Capture
>        * CRISPR Guide Capture
>        * Antigen Capture
>        * Custom
>
> ***Example:*** `--features features.csv`

### 2.3 Analysis Options

`--exclude-introns`
> **Exclude introns from the count alignment.**   
> *type: boolean flag*
>   
> Turns off the option of including introns when performing alignment. This flag is only applicable when dealing with gene expression related data.
>
> ***Example:*** `--exclude-introns`


---  
  `--filter FILTER`
> **Filter threshold file.**   
> *type: file*
>   
> Filter threshold file. A CSV file containing the different thresholds to be applied for individual samples within the project during the QC analysis. The file should contain a header row with Sample as the column name for the sample IDs, and the name of each metric that will be filtered along with if it is the high or low threshold for that metric. Each row is then the entries for each sample that the manual thresholds will be applied. If no file is provided then the default thresholds will be used. If a cell is left blank for a sample then that sample would not be filtered based on that criteria. 
>
> *Here is an example filter.csv file:*
> ```
> Sample,nFeature_RNA_low,nFeature_RNA_high,percent.mito_high
> sample1,500,6000,15
> sample2,500,6000,5
> sample4,500,6000,5
> ```
>
> *Where:*  
>
> - *Sample:* Unique sample ID that should match the sample name used for Cell Ranger count.
> - *nFeature_RNA_low,nFeature_RNA_high,percent.mito_high:* Example entries that can be used for manual thresholding. The column names need to be formatted as metadataname_high/low. Entries that ends with high will be treated as the upper threshold. Entries that ends with low will be treated as the lower threshold. Valid metadata names include nCount_RNA, nFeature_RNA, and percent.mito.
>
> ***Example:*** `--filter filter.csv`

---  
  `--metadata METADATA`
> **Sample metadata file.**   
> *type: file*
>   
> Sample metadata file. A CSV file containing sample level metadata information that will be included as new metadata columns during QC analysis. The file should contain a header row with Sample as the column name for the sample IDs, and the name of each metadata column that will be added to their associated samples. Each row is then the entries for each sample with the values that will be included as metadata. If no file is provided then no metadata will be added to the samples. If a cell is left blank for a sample then the metadata column for that sample would be an empty string. 
>
> *Here is an example metadata.csv file:*
> ```
> Sample,Type,batch
> sample1,tumor,1
> sample2,normal,1
> sample4,tumor,2
> ```
>
> *Where:*  
>
> - *Sample:* Unique sample ID that should match the sample name used for Cell Ranger count.
> - *Type,batch:* Example metadata entries that can be applied to the created Seurat object. The column names will be the resulting metadata entries, so it is recommended to use ones that do not overlap with column names that would be created automatically.
>
> ***Example:*** `--metadata metadata.csv`


---  
  `--create-bam`
> **Create bam files.**   
> *type: boolean flag*
>   
> By default the no-bam flag is used when running Cell Ranger. Use this flag to ensure that a bam file is created for each sample during analysis. This flag is only applicable when dealing with gene expression related data.
>
> ***Example:*** `--create-bam`

---
`--forcecells FORCECELLS`
> **Force cells file.**  
> *type: file*
>
> Force cells file. A CSV file containing the name of the sample (the Cell Ranger outputted name) and the number of cells to force the sample to. It will generally be used if the first analysis run appears to do a poor job at estimating the number of cells, and a re-run is needed to adjust the number of cells in the sample.
>
> *Here is an example forcecells.csv file:*
> ```
> Sample,Cells
> Sample1,3000
> Sample2,5000
> ```
>
> *Where:*
>
> - *Sample:* The sample name used as the Cell Ranger output
> - *Cells:* The number of cells the sample should be forced to
>
> In this example, Sample1 and Sample2 will be run while being forced to have 3000 and 5000 cells respectively. Any other samples that are processed will be run without using the force cells flag and will use the default cell calling algorithm.
>
> ***Example:*** `--forcecells forcecells.csv`


### 2.4 Orchestration Options

Each of the following arguments are optional, and do not need to be provided.

  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean flag*
>
> Displays what steps in the pipeline remain or will be run. Does not execute anything!
>
> ***Example:*** `--dry-run`

---  
  `--silent`            
> **Silence standard output.**  
> *type: boolean flag*
>
> Reduces the amount of information directed to standard output when submitting master job to the job scheduler. Only the job id of the master job is returned.
>
> ***Example:*** `--silent`

---  
  `--mode {slurm,local}`  
> **Execution Method.**  
> *type: string*  
> *default: slurm*
>
> Execution Method. Defines the mode or method of execution. Vaild mode options include: slurm or local.
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
> *default: pl:cell-seek*
>
> When submitting the pipeline to a job scheduler, like SLURM, this option always you to set the name of the pipeline's master job. By default, the name of the pipeline's master job is set to "pl:cell-seek".
>
> ***Example:*** `--job-name pl_id-42`

---  
  `--singularity-cache SINGULARITY_CACHE`  
> **Overrides the $SINGULARITY_CACHEDIR environment variable.**  
> *type: path*  
> *default: `--output OUTPUT/.singularity`*
>
> Singularity will cache image layers pulled from remote registries. This ultimately speeds up the process of pull an image from DockerHub if an image layer already exists in the singularity cache directory. By default, the cache is set to the value provided to the `--output` argument. Please note that this cache cannot be shared across users. Singularity strictly enforces you own the cache directory and will return a non-zero exit code if you do not own the cache directory! See the `--sif-cache` option to create a shareable resource.
>
> ***Example:*** `--singularity-cache /data/$USER/.singularity`

---  
  `--sif-cache SIF_CACHE`
> **Path where a local cache of SIFs are stored.**  
> *type: path*  
>
> Uses a local cache of SIFs on the filesystem. This SIF cache can be shared across users if permissions are set correctly. If a SIF does not exist in the SIF cache, the image will be pulled from Dockerhub and a warning message will be displayed. The `cell-seek cache` subcommand can be used to create a local SIF cache. Please see `cell-seek cache` for more information. This command is extremely useful for avoiding DockerHub pull rate limits. It also remove any potential errors that could occur due to network issues or DockerHub being temporarily unavailable. We recommend running cell-seek with this option when ever possible.
>
> ***Example:*** `--singularity-cache /data/$USER/SIFs`

---  
  `--threads THREADS`   
> **Max number of threads for each process.**  
> *type: int*  
> *default: 2*
>
> Max number of threads for each process. This option is more applicable when running the pipeline with `--mode local`.  It is recommended setting this vaule to the maximum number of CPUs available on the host machine.
>
> ***Example:*** `--threads 12`


---  
  `--tmp-dir TMP_DIR`   
> **Max number of threads for each process.**  
> *type: path*  
> *default: `/lscratch/$SLURM_JOBID`*
>
> Path on the file system for writing temporary output files. By default, the temporary directory is set to '/lscratch/$SLURM_JOBID' for backwards compatibility with the NIH's Biowulf cluster; however, if you are running the pipeline on another cluster, this option will need to be specified. Ideally, this path should point to a dedicated location on the filesystem for writing tmp files. On many systems, this location is set to somewhere in /scratch. If you need to inject a variable into this string that should NOT be expanded, please quote this options value in single quotes.
>
> ***Example:*** `--tmp-dir /scratch/$USER/`

### 2.5 Miscellaneous Options  
Each of the following arguments are optional, and do not need to be provided.

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
>
> Shows command's synopsis, help message, and an example command
>
> ***Example:*** `--help`

## 3. Example

### 3.1 Basic Pipeline Usage

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline cite \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline cite \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --mode slurm
```

### 3.2 Run Downstream on Existing Cell Ranger Output

It is possible to use cell-seek to perform the initial downstream analysis on existing Cell Ranger output. The files are expected to be in the Cell Ranger outputted format with the outs folder present.  The sample level folders should be provided as the input for cell-seek.

For example, if sample1 was run in Cell Ranger then sample1/outs/ contains the final pipeline output files, and sample1 should be provided as input to cell-seek.

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*/ \
                  --output /data/$USER/output \
                  --pipeline cite \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*/ \
                  --output /data/$USER/output \
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --mode slurm
```

### 3.3 Run with Custom Filters

In order to use this option, a CSV file should be created. In the CSV there should be a row for each sample that will have custom filters applied. If these custom filters are being applied to a project that was already processed by the cell-seek pipeline, then the created `Project_Cell_Filters.csv` can be used as a template.

The following is an example of a created filter CSV file.

```
Sample,nFeature_RNA_low,nFeature_RNA_high,percent.mito_high
sample1,500,6000,15
sample2,500,6000,5
sample4,500,6000,
```

Based on this file, custom filters will only be applied to samples sample1, sample2, and sample4. These samples will have the following filters applied to only keep cells meeting the following criteria:

- sample1:  nFeature_RNA > 500 and nFeature_RNA < 6000 and percent.mito < 15
- sample2:  nFeature_RNA > 500 and nFeature_RNA < 6000 and percent.mito < 5
- sample4:  nFeature_RNA > 500 and nFeature_RNA < 6000


```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline cite \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --filter filter.csv
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline cite \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --filter filter.csv
                  --mode slurm
```

### 3.3 Adding Sample Level Metadata

Sample level information can be provided in a CSV file, which would be added to the generated Seurat objects. In the CSV file there should be a row for each sample, and a column for each aspect of information that would be added. If a cell is left blank for a sample then the metadata column for that sample would be an empty string.

The following is an example of a created metadata CSV file.

```
Sample,Type,batch
sample1,tumor,1
sample2,normal,1
sample4,tumor,2
```

Based on this file, columns named Type, batch, and tissue will be added to the generated Seurat objects for sample1, sample2, and sample4. Any other samples processed in the pipeline would not have metadata columns added. 

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline cite \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --metadata metadata.csv
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline cite \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --metadata metadata.csv
                  --mode slurm
```

### 3.4 Creating BAM Files

Cell Ranger BAM file generation is turned off by default to save space. However, if the BAM file is required for downstream analysis, then it will need to be turned on.

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline cite \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --create-bam
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline cite \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --create-bam
                  --mode slurm
```

### 3.5 Running while Forcing Cell Call


It is possible to force the number of cells that are called by Cell Ranger. In this situation, Cell Ranger will call the top X cell barcodes with the highest UMI count as cells, where X the number cells that the sample is forced to. This is generally used if the first analysis run appears to do a poor job at estimating the number of cells, and a re-run while adjusting the number of cells in the sample is helpful.

A CSV file needs to be created with the first column containing the name of the sample (the Cell Ranger outputted name) and the second column containing the number of cells to force the sample to. Only the samples included in the CSV file will be run while forcing the cell call. Any other samples that are processed will use the default cell calling algorithm.

The following is an example of a force cells CSV file.

```
Sample,Cells
Sample1,3000
Sample2,5000
```

Based on this file, Sample1 and Sample 2 will be run while being forced to have 3000 and 5000 cells respectively.

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline cite \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --forcecells forcecells.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline cite \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --forcecells forcecells.csv \
                  --mode slurm
```