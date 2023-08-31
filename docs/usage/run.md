# <code>cell-seek <b>run</b></code>

## 1. About 
The `cell-seek` executable is composed of several inter-related sub commands. Please see `cell-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>cell-seek <b>run</b></code> sub command in more detail. With minimal configuration, the **`run`** sub command enables you to start running cell-seek pipeline. 

Setting up the cell-seek pipeline is fast and easy! In its most basic form, <code>cell-seek <b>run</b></code> only has *two required inputs*.

## 2. Synopsis
```text
$ cell-seek run [--help] \
      [--mode {slurm,local}] [--job-name JOB_NAME] [--batch-id BATCH_ID] \
      [--tmp-dir TMP_DIR] [--silent] [--sif-cache SIF_CACHE] \ 
      [--singularity-cache SINGULARITY_CACHE] \
      [--dry-run] [--threads THREADS] \
      [--libraries LIBRARIES] [--features FEATURES] \
      [--cmo-reference CMOREFERENCE] [--cmo-sample CMOSAMPLE] \
      [--exclude-introns] \
      --input INPUT [INPUT ...] \
      --output OUTPUT \
      --version {gex, ...} \
      --genome {hg38, ...}
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of FastQ (globbing is supported) to analyze via `--input` argument, an output directory to store results via `--output` argument, the version of the pipeline to run via `--version` argument, and the reference genome to use via `--genome` argument.

Use you can always use the `-h` option for information on a specific command. 

### 2.1 Required arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input FastQ file(s).**  
> *type: file(s)*  
> 
> One or more FastQ files can be provided. The pipeline does NOT support single-end data. From the command-line, each input file should seperated by a space. Globbing is supported! This makes selecting FastQ files easy. Input FastQ files should always be gzipp-ed.
> 
> ***Example:*** `--input .tests/*.R?.fastq.gz`

---  
  `--output OUTPUT`
> **Path to an output directory.**   
> *type: path*
>   
> This location is where the pipeline will create all of its output files, also known as the pipeline's working directory. If the provided output directory does not exist, it will be created automatically.
> 
> ***Example:*** `--output /data/$USER/cell-seek_out`

---  
  `--version {gex, vdj, cite, multi, atac, multiome}`
> **The version of the pipeline to run.**   
> *type: string*
>   
> This option selects the version of the pipeline to run. Can be chosen from the options gene expression (GEX), VDJ, feature barcode (CITE), Cell Ranger multi (Multi), ATAC, or multiome.
> 
> ***Example:*** `--version gex`

---  
  `--genome {hg38, mm10}`
> **Reference genome.**   
> *type: string*
>   
> This option defines the reference genome of the samples. cell-seek does comes bundled with prebuilt reference files for human and mouse samples, e.g. hg38 or mm10. Please select one of the following options: hg38, mm10
> 
> ***Example:*** `--genome hg38`

### 2.2 Analysis options

Each of the following arguments are optional, and do not need to be provided. 

  `--libraries LIBRARIES`
> **Libraries file.**   
> *type: file*
>   
> A CSV file containing information about each library. This file is used in feature barcode (cite), multi, and multiome analysis. It contains each sample's name, flowcell, demultiplexed name, and library type. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).

> *Here is an example libraries.csv file:*
> ```
> Name,Flowcell,Sample,Type
> IL15_LNs,H7CNNBGXG,IL15_LNs,Gene Expression
> IL15_LNs,H7CT7BGXG,IL15_LNs_BC,Antibody Capture
> ``` 

> *Where:* 

> - *Name:* name of the sample passed to CellRanger.  
> - *Flowcell:* The flowcell ID that contains the FASTQ files for this set of data.  
> - *Sample:* Name that was used when demultiplexing, this should match the FASTQ files.  
> - *Type:* library type for each sample. List of supported options:  
>        * Gene Expression
>        * CRISPR Guide Capture
>        * Antibody Capture
>        * Custom

> ***Example:*** `--libraries libraries.csv`


---  
  `--features FEATURES`
> **Features file.**   
> *type: file*
>   
> A feature reference CSV file containing information for processing a feature barcode data. This file is used in feature barcode, and may be used in multi analysis. This file should contain a unique ID for the feature, a human readable name, sequence, feature type, read, and pattern. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).

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
>        * Gene Expression
>        * CRISPR Guide Capture
>        * Antibody Capture
>        * Custom

> ***Example:*** `--features features.csv`


---  
  `--cmo-reference CMOREFERENCE`
> **CMO reference file.**   
> *type: file*
>   
> A CMO reference CSV file containing information for processing hashtag data, which is used in multi analysis if custom hashtags will be processed. This file should contain a unique ID for the hashtag, a human readable name, sequence, feature type, read, and pattern. More information about the cmo reference file and its requirements can be found on the [10x Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi).

> *Here is an example cmo_reference.csv file:*
> ```
> id,name,sequence,feature_type,read,pattern
> CMO301,CMO301,ATGAGGAATTCCTGC,Multiplexing Capture,R2,5P(BC)
> CMO302,CMO302,CATGCCAATAGAGCG,Multiplexing Capture,R2,5P(BC)
> ``` 

> *Where:*  

> - *id:* Unique ID for this feature. Must not contain whitespace, quote or comma characters. Each ID must be unique and must not collide with a gene identifier from the transcriptome. 
> - *name:* Human-readable name for this feature. Must not contain whitespace. 
> - *sequence:* Nucleotide barcode sequence associated with this hashtag
> - *feature_type: Type of the feature. This should always be multiplexing capture. 
> - *read:* Specifies which RNA sequencing read contains the Feature Barcode sequence. Must be R1 or R2, but in most cases R2 is the correct read. 
> - *pattern:* Specifies how to extract the sequence of the feature barcode from the read.
> ***Example:*** `--cmo-reference cmo_reference.csv`

---  
  `--cmo-sample CMOSAMPLE`
> **CMO sample file.**   
> *type: file*
>   
> A CMO sample CSV file containing sample to hashtag information used for multi analysis. CMO IDs should match the ones used in the CMO reference file. If no CMO reference is provided then CMO ID should match the ones used by 10x. The same CMO sample will be used on all multi libraries. This file should contain a unique sample ID for the sample, and the CMO IDs associated with that sample. If more than one CMO ID is associated with a sample then a | should be used to separate the tags. More information and examples about the samples section of the multi config file and its requirements can be found on the [10x Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi).

> *Here is an example cmo_sample.csv file:*
> ```
> sample_id,cmo_ids
> sample1,CMO301
> sample2,CMO302|CMO303
> ``` 

> *Where:*  

> - *sample_id:* Unique sample ID for this hashtagged sample. Must not contain whitespace, quote or comma characters. Each sample ID must be unique. 
> - *cmo_ids:* Unique CMO ID(s) that the sample is hashtagged with. Must match either entries in the cmo_reference.csv file or 10x CMO IDs.
> ***Example:*** `--cmo-sample cmo_sample.csv`

---  
  `--exclude-introns`
> **Exclude introns from the count alignment.**   
> *type: boolean flag*
>   
> Turns off the option of including introns when performing alignment. This flag is only applicable when dealing with gene expression related data.
> ***Example:*** `--exclude-introns`
                                      
### 2.3 Orchestration options

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

### 2.4 Miscellaneous options  
Each of the following arguments are optional, and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

## 3. Example
```bash 
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --version gex \
                  --genome hg38 \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to 
# the cluster. It is recommended running 
# the pipeline in this mode.
./cell-seek run --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --version gex \
                  --genome hg38 \
                  --mode slurm
```
