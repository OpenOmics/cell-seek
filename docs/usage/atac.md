# <code>cell-seek <b>run</b> --pipeline atac</code>

## 1. About
The `cell-seek` executable is composed of several inter-related sub commands. Please see `cell-seek -h` for all available options.

This part of the documentation describes options and concepts for [<code>cell-seek <b>run</b></code>](run.md) sub command for the **ATAC** pipeline which an be selected via the <code><b>--pipeline</b></code> flag in more detail. With minimal configuration, the **`run`** sub command enables you to start running cell-seek pipeline.

Setting up the cell-seek pipeline is fast and easy! In its most basic form, <code>cell-seek <b>run</b></code> only has *five required inputs*.

### 1.1 Use Case

The ATAC pipeline should be used when only ATAC data will be analyzed. If gene expression data was captured and will be analyzed alongside ATAC, then [the Multiome pipeline](multiome.md) should be used.

A basic guideline of which pipeline should be used for different modalities can be found in the [synopsis section of the <code><b>run</b></code> documentation.](run.md#21-pipelines)


## 2. Synopsis
```text
$ cell-seek run [--help] \
      [--dry-run] [--job-name JOB_NAME] [--mode {{slurm,local}}] \
      [--sif-cache SIF_CACHE] [--singularity-cache SINGULARITY_CACHE] \
      [--silent] [--threads THREADS] [--tmp-dir TMP_DIR] \
      [--rename RENAME] [--forcecells FORCECELLS] \
      --input INPUT [INPUT ...] \
      --output OUTPUT \
      --pipeline atac \
      --genome {hg38, ...} \
      --cellranger {2.1.0, 2.2.0}
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of FastQ (globbing is supported) to analyze via `--input` argument, an output directory to store results via `--output` argument, the version of the pipeline to run via `--pipeline` argument, the reference genome to use via `--genome` argument, and the version of cellranger-atac to use via `--cellranger` argument.

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
  `--pipeline atac`
> **The pipeline to run.**   
> *type: string*
>   
> This option selects the version of the pipeline to run. The documentation provided is based on choosing the option for ATAC.
>
> ***Example:*** `--pipeline atac`

---  
  `--genome {hg38, mm10, hg2024, mm2024, custom.json}`
> **Reference genome.**   
> *type: string*
>   
> This option defines the reference genome of the samples. cell-seek does comes bundled with prebuilt reference files for human and mouse samples, e.g. hg38 or mm10.
>
> A custom reference genome can also be provided via a json file. Additional information for creating this json file can be found in [<code>cell-seek <b>genome</b></code>](genome.md).
>
> For prebuilt references please select one of the following options: hg38, mm10
>
> ***Example:*** `--genome hg38`

---  
  `--cellranger {2.1.0, 2.2.0}`
> **The version of Cell Ranger ATAC to run.**   
> *type: string*
>   
> This option specifies which version of Cell Ranger ATAC to use when running the ATAC pipeline. Please select one of the following options: 2.1.0, 2.2.0
>
> ***Example:*** `--cellranger 2.2.0`


### 2.2 Analysis Options
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
                  --pipeline atac \
                  --genome hg38 \
                  --cellranger 2.2.0 \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline atac \
                  --genome hg38 \
                  --cellranger 2.2.0 \
                  --mode slurm
```

#### 3.1.2 Multiple Input Folders

FASTQ files from multiple folders can be provided as input. The paths to different FASTQ files can be provided in a space separated format.

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*_R?_fastq.gz .tests2/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline atac \
                  --genome hg38 \
                  --cellranger 2.2.0 \
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
                  --cellranger 2.2.0 \
                  --mode slurm
```

### 3.2 Run Downstream on Existing Cell Ranger Output

It is possible to use cell-seek to perform the initial downstream analysis on existing Cell Ranger ATAC output. The files are expected to be in the Cell Ranger ATAC outputted format with the outs folder present.  The sample level folders should be provided as the input for cell-seek.

For example, if sample1 was run in Cell Ranger ATAC then sample1/outs/ contains the final pipeline output files, and sample1 should be provided as input to cell-seek.

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*/ \
                  --output /data/$USER/output \
                  --pipeline atac \
                  --genome hg38 \
                  --cellranger 2.2.0 \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*/ \
                  --output /data/$USER/output \
                  --pipeline atac \
                  --genome hg38 \
                  --cellranger 2.2.0 \
                  --mode slurm
```


### 3.3 Renaming Samples

It is possible to rename samples to different names from the ones that are used in the FASTQ files. This function can be used to change the sample names to something more informative, or it could be used if FASTQ files names changed for a specific sample.  

In order to use this option, a CSV file should be created. In the CSV there should be a row for each FASTQ file that will be processed. The first column contains the name of the FASTQ file while the second column contains the output sample name. Only the FASTQ files and samples that are listed within the CSV file will be processed.

Cell Ranger expects the FASTQ file to have a format of

`[Sample Name]_S[Sample Number]_L00[Lane Number]_[Read Type]_001.fastq.gz`

or

`[Sample Name]_S[Sample Number]_[Read Type]_001.fastq.gz`

The FASTQ name that is used in the rename CSV file should match the section in the `[Sample Name]` listed above. An example file could be `Sample_ATAC_S1_L001_R1_001.fastq.gz`. In this file, the sample name would be Sample_ATAC. More information about the FASTQ naming format that Cell Ranger expects can be found at the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/analysis/inputs/specifying-input-fastq-files).

The following is a potential example of a rename CSV file.

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
                  --pipeline atac \
                  --genome hg38 \
                  --cellranger 2.2.0 \
                  --rename rename.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline atac \
                  --genome hg38 \
                  --cellranger 2.2.0 \
                  --rename rename.csv \
                  --mode slurm
```

### 3.4 Running while Forcing Cell Call

It is possible to force the number of cells that are called by Cell Ranger ATAC. In this situation, Cell Ranger ATAC will call the top X cell barcodes with the highest number of fragments overlapping peaks as cells, where X the number cells that the sample is forced to. This is generally used if the first analysis run appears to do a poor job at estimating the number of cells, and a re-run while adjusting the number of cells in the sample is helpful.

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
                  --pipeline atac \
                  --genome hg38 \
                  --cellranger 2.2.0 \
                  --forcecells forcecells.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline atac \
                  --genome hg38 \
                  --cellranger 2.2.0 \
                  --forcecells forcecells.csv \
                  --mode slurm
```