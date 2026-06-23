# <code>cell-seek <b>run</b> --pipeline gex</code>

## 1. About
The `cell-seek` executable is composed of several inter-related sub commands. Please see `cell-seek -h` for all available options.

This part of the documentation describes options and concepts for [<code>cell-seek <b>run</b></code>](run.md) sub command for the **GEX** pipeline which an be selected via the <code><b>--pipeline</b></code> flag in more detail. With minimal configuration, the **`run`** sub command enables you to start running cell-seek pipeline.

Setting up the cell-seek pipeline is fast and easy! In its most basic form, <code>cell-seek <b>run</b></code> only has *five required inputs*.

### 1.1 Use Case

The GEX pipeline should be used when only single cell RNA (scRNA) gene expression data will be analyzed. If other modalities were captured and will be analyzed alongside scRNA, then one of the other pipelines should be used.

A basic guideline of which pipeline should be used for different modalities can be found in the [synopsis section of the <code><b>run</b></code> documentation.](run.md#21-pipelines)

## 2. Synopsis
```text
$ cell-seek run [--help] \
      [--dry-run] [--job-name JOB_NAME] [--mode {{slurm,local}}] \
      [--sif-cache SIF_CACHE] [--singularity-cache SINGULARITY_CACHE] \
      [--silent] [--threads THREADS] [--tmp-dir TMP_DIR] \
      [--aggregate {{mapped, none}}] [--exclude-introns] \
      [--filter FILTER] [--metadata METADATA] [--create-bam] \
      [--rename RENAME] [--forcecells FORCECELLS] \
      --input INPUT [INPUT ...] \
      --output OUTPUT \
      --pipeline gex \
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
  `--pipeline gex`
> **The pipeline to run.**   
> *type: string*
>   
> This option selects the pipeline to run. The documentation provided is based on choosing the option for gene expression (GEX).
>
> ***Example:*** `--pipeline gex`

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
> ***Example:*** `--genome hg2024`

---  
  `--cellranger {7.1.0, 7.2.0, 8.0.0, 9.0.0, 10.0.0}`
> **The version of Cell Ranger to run.**   
> *type: string*
>   
> This option specifies which version of Cell Ranger to use when running GEX, VDJ, CITE, or MULTI pipelines. Please select one of the following options: 7.1.0, 7.2.0, 8.0.0, 9.0.0, 10.0.0
>
> ***Example:*** `--cellranger 7.1.0`

### 2.2 Analysis Options

Each of the following arguments are optional, and do not need to be provided.

`--aggregate {mapped, none}`
> **Cell Ranger aggregate normalization.**   
> *type: string*
>  
> This option defines the normalization mode that should be used. Mapped is what Cell Ranger would run by default, which subsamples reads from higher depth samples until each library type has an equal number of reads per cell that are confidently mapped.  None means to not normalize at all. If this flag is not used then aggregate will not be run. Aggregate analysis is generally not needed, but it can be used to generate a Loupe Browser file for interactive exploration of the data. To run Cell Ranger aggregate, please select one of the following options: mapped, none.
>
> ***Example:*** `--aggregate mapped`


---  
  `--exclude-introns`
> **Exclude introns from the count alignment.**   
> *type: boolean flag*
>   
> Turns off the option of including introns when performing alignment. This flag is only applicable when dealing with gene expression related data.
>
> ***Example:*** `--exclude-introns`


---  
  `--create-bam`
> **Create bam files.**   
> *type: boolean flag*
>   
> By default the no-bam flag is used when running Cell Ranger. Use this flag to ensure that a bam file is created for each sample during analysis. This flag is only applicable when dealing with gene expression related data.
>
> ***Example:*** `--create-bam`

---  
  `--filter FILTER`
> **Filter threshold file.**   
> *type: file*
>   
> Filter threshold file. A CSV file containing the different thresholds to be applied for individual samples within the project during the QC analysis. The file should contain a header row with Sample as the column name for the sample IDs, and the name of each metric that will be filtered along with if it is the high or low threshold for that metric. Each row is then the entries for each sample that the manual thresholds will be applied. If no file is provided then the default thresholds will be used. If a cell is left blank for a sample then that sample would not be filtered based on that criteria. This flag is currently only applicable when dealing with GEX projects.
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
  `--rename RENAME`
> **Rename sample file.**   
> *type: file*
>
> Rename sample file. A CSV file containing the name of the FASTQ file and the new name of the sample. Only the samples listed in the CSV files will be run.
>
> *Here is an example rename.csv file:*
> ```
> FASTQ,Name
> original_name1,new_name1
> original_name2,new_name2
> original_name3,new_name3
> original_name3-2,new_name3
> original_name4,original_name4
> ```
>
> *Where:*
>
> - *FASTQ:* The name that is used in the FASTQ file
> - *Name:* Unique sample ID that is the sample name used for Cell Ranger count.
>
> In this example, new_name3 has FASTQ files with two different names. With this input, both sets of FASTQ files will be used when processing the sample as new_name3. original_name4 will not be renamed. Any FASTQ file that does not have the name original_name1, original_name2, original_name3, or original_name4 will not be run.
>
> ***Example:*** `--rename rename.csv`

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

### 2.3 Orchestration Options

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

### 3.1 Basic Pipeline Usage

#### 3.1.1 Single Input Folder

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --mode slurm
```

#### 3.1.2 Multiple Input Folders

FASTQ files from multiple flowcells can be provided as input. The paths to different FASTQ files can be provided in a space separated format.

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*_R?_fastq.gz .tests2/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz .tests2/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
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
                  --pipeline gex \
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
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --filter filter.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --filter filter.csv \
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
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --metadata metadata.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --metadata metadata.csv \
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
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --create-bam
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --create-bam
                  --mode slurm
```

### 3.5 Renaming Samples

It is possible to rename samples to different names from the ones that are used in the FASTQ files. This function can be used to change the sample names to something more informative, or it could be used if FASTQ files names changed for a specific sample.  

In order to use this option, a CSV file should be created. In the CSV there should be a row for each FASTQ file that will be processed. The first column contains the name of the FASTQ file while the second column contains the output sample name. Only the FASTQ files and samples that are listed within the CSV file will be processed.

Cell Ranger expects the FASTQ file to have a format of

`[Sample Name]_S[Sample Number]_L00[Lane Number]_[Read Type]_001.fastq.gz`

or

`[Sample Name]_S[Sample Number]_[Read Type]_001.fastq.gz`

The FASTQ name that is used in the rename CSV file should match the section in the `[Sample Name]` listed above. An example file could be `Sample_GEX_S1_L001_R1_001.fastq.gz`. In this file, the sample name would be Sample_GEX. More information about the FASTQ naming format that Cell Ranger expects can be found at the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-specifying-fastqs).

The following is an example of a rename CSV file.

```
FASTQ,Name
sample1_run1,sample1
sample2_run1,sample2
sample3_run1,sample3
sample3_run2,sample3
sample4,sample4
```

Based on this file, the FASTQ files with the name sample1_run1, sample2_run1, sample3_run1, sample3_run2, and sample4 will be processed, and any FASTQ file that does not match these names will not be run. The outputted sample names would be sample1, sample2, sample3, and sample4. The FASTQ files with the name sample3_run1 and sample3_run2 will be processed as sample3. 


```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --rename rename.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --rename rename.csv \
                  --mode slurm
```

### 3.6 Running while Forcing Cell Call

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
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --forcecells forcecells.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --forcecells forcecells.csv \
                  --mode slurm
```