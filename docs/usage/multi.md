# <code>cell-seek <b>run</b> --pipeline multi</code>

## 1. About
The `cell-seek` executable is composed of several inter-related sub commands. Please see `cell-seek -h` for all available options.

This part of the documentation describes options and concepts for [<code>cell-seek <b>run</b></code>](run.md) sub command for the **MULTI** pipeline which an be selected via the <code><b>--pipeline</b></code> flag in more detail. With minimal configuration, the **`run`** sub command enables you to start running cell-seek pipeline.

Setting up the cell-seek pipeline is fast and easy! In its most basic form, <code>cell-seek <b>run</b></code> only has *five required inputs*.

### 1.1 Use Case

The MULTI pipeline uses the Cell Ranger multi function. There are multiple situation where this should be used. These include the following situations

- Any combination that contains both gene expression and VDJ data 
- Any time sample multiplexing was performed and Cell Ranger will be used for demultiplexing
    - This can include on-chip multiplexing (OCM), hashing with antibody capture (HTO), and CellPlex (CMO)
- Flex (fixed RNA) capture 

If the data does not include these modalities, then another pipeline should be considered.

A basic guideline of which pipeline should be used for different modalities can be found in the [synopsis section of the <code><b>run</b></code> documentation.](run.md#21-pipelines)


## 2. Synopsis
```text
$ cell-seek run [--help] \
      [--dry-run] [--job-name JOB_NAME] [--mode {{slurm,local}}] \
      [--sif-cache SIF_CACHE] [--singularity-cache SINGULARITY_CACHE] \
      [--silent] [--threads THREADS] [--tmp-dir TMP_DIR] \
      [--exclude-introns] \
      [--libraries LIBRARIES] [--features FEATURES] \
      [--cmo-sample CMOSAMPLE] [--cmo-reference CMOREFERENCE] \
      [--hto-sample HTOSAMPLE] \
      [--ocm-sample OCMSAMPLE] \
      [--probe-set PROBESET] [--probe-sample PROBESAMPLE] \
      [--filter FILTER] [--metadata METADATA] [--create-bam] \
      [--forcecells FORCECELLS] \
      --input INPUT [INPUT ...] \
      --output OUTPUT \
      --pipeline multi \
      --genome {hg38, ...} \
      --cellranger {8.0.0, ...}
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of FastQ (globbing is supported) to analyze via `--input` argument, an output directory to store results via `--output` argument, the `multi` pipeline via `--pipeline` argument, the reference genome to use via `--genome` argument, and the version of cellranger to use via `--cellranger` argument.

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
  `--pipeline multi`
> **The pipeline to run.**   
> *type: string*
>   
> This option selects the version of the pipeline to run. The documentation provided is based on choosing the option for multi analysis.
>
> ***Example:*** `--pipeline multi`

---  
  `--genome {hg38, mm10, hg2024, mm2024, custom.json}`
> **Reference genome.**   
> *type: string*
>   
> This option defines the reference genome of the samples. cell-seek does comes bundled with prebuilt reference files for human and mouse samples, The options hg38 or mm10 would select the 2020 release of the reference. The options hg2024 or mm2024 would select the 2024 release of the reference. More information about the officially released references can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-reference-release-notes). Since there is no 2024 released VDJ reference, if hg2024 or mm2024 is selected in a run that includes VDJ data, the VDJ reference CR 7.1 release will be used for human, and CR 7.0 release will be used for mouse.
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
> IL15_LNs,H7CT7BGXG,IL15_LNs_BC,Antibody Capture
> ```

> *Where:*

> - *Name:* name of the sample passed to Cell Ranger.  
> - *Flowcell:* A unique identifier in the path that the FASTQ files are located.
> - *Sample:* Name that was used when demultiplexing, this should match the FASTQ files.  
> - *Type:* library type for each sample. List of supported options:  
>        * Gene Expression
>        * CRISPR Guide Capture
>        * Antibody Capture
>        * Multiplexing Capture
>        * VDJ
>        * VDJ-B
>        * VDJ-T
>        * Custom
>
> ***Example:*** `--libraries libraries.csv`

### 2.3 Demultiplexing Options

This section contains the flags that are applicable when sample multiplexing was performed. Different flags should be used depending on the multiplexing method performed.

#### 2.3.1 CellPlex (CMO) Options

`--cmo-sample CMOSAMPLE`
> **CMO sample file.**   
> *type: file*
>   
> A CMO sample CSV file containing sample to CMO information used for multi analysis. CMO IDs should match the ones used in the CMO reference file. If no CMO reference is provided then CMO ID should match the ones used by 10x. The same CMO sample will be used on all multi libraries. This file should contain a unique sample ID for the sample, and the CMO ID(s) associated with that sample. If more than one CMO ID is associated with a sample then a | should be used to separate the tags. More information and examples about the samples section of the multi config file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

> *Here is an example cmo_sample.csv file:*
> ```
> sample_id,cmo_ids
> sample1,CMO301
> sample2,CMO302|CMO303
> ```

> *Where:*  

> - *sample_id:* Unique sample ID for this tagged sample. Must not contain whitespace, quote or comma characters. Each sample ID must be unique.
> - *cmo_ids:* Unique CMO ID(s) that the sample is tagged with. Must match either entries in the cmo_reference.csv file or 10x CMO IDs.
>
> ***Example:*** `--cmo-sample cmo_sample.csv`

---
`--cmo-reference CMOREFERENCE`
> **CMO reference file.**   
> *type: file*
>   
> A CMO reference CSV file containing information for processing CMO data, which is used in multi analysis if custom Cell Multiplexing oligos will be processed. This file should contain a unique ID for the CMO, a human readable name, sequence, feature type, read, and pattern. More information about the cmo reference file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

> *Here is an example cmo_reference.csv file:*
> ```
> id,name,sequence,feature_type,read,pattern
> CMO301,CMO301,ATGAGGAATTCCTGC,Multiplexing Capture,R2,5P(BC)
> CMO302,CMO302,CATGCCAATAGAGCG,Multiplexing Capture,R2,5P(BC)
> ```

> *Where:*  

> - *id:* Unique ID for this feature. Must not contain whitespace, quote or comma characters. Each ID must be unique and must not collide with a gene identifier from the transcriptome.
> - *name:* Human-readable name for this feature. Must not contain whitespace.
> - *sequence:* Nucleotide barcode sequence associated with this CMO
> - *feature_type:* Type of the feature. This should always be multiplexing capture.
> - *read:* Specifies which RNA sequencing read contains the Feature Barcode sequence. Must be R1 or R2, but in most cases R2 is the correct read.
> - *pattern:* Specifies how to extract the sequence of the feature barcode from the read.
>
> ***Example:*** `--cmo-reference cmo_reference.csv`

#### 2.3.2 Hashing with Antibody Capture (HTO) Arguments

`--hto-sample HTOSAMPLE`
> **HTO sample file.**   
> *type: file*
>   
> A HTO sample CSV file containing sample to hashtag information used for multi analysis. This should always be used in conjunction with the feature file provided using the <code>feature</code> flag. More information about that flag can be found under [Analysis Options](multi.md#25-analysis-options).

> HTO IDs should match the ones used in the feature reference file. The same HTO sample will be used on all multi libraries. This file should contain a unique sample ID for the sample, and the HTO ID(s) associated with that sample. If more than one HTO ID is associated with a sample then a | should be used to separate the tags. More information and examples about the samples section of the multi config file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

> *Here is an example hto_sample.csv file:*
> ```
> sample_id,hashtag_ids
> sample1,HTO1
> sample2,HTO2|HTO3
> ```

> *Where:*  

> - *sample_id:* Unique sample ID for this hashtagged sample. Must not contain whitespace, quote or comma characters. Each sample ID must be unique.
> - *hashtag_ids:* Unique HTO ID(s) that the sample is hashtagged with. Must match entries in the features.csv file.
>
> ***Example:*** `--hto-sample hto_sample.csv`

#### 2.3.3 On-Chip Multiplexing (OCM) Arguments

`--ocm-sample OCMSAMPLE`
> **OCM sample file.**   
> *type: file*
>   
> An OCM sample CSV file containing sample to OCM barcode used for multi analysis. The same OCM sample will be used on all multi libraries. This file should contain a unique sample ID for the sample, and the OCM ID(s) associated with that sample. If more than one OCM ID is associated with a sample then a | should be used to separate the tags. More information and examples about the samples section of the multi config file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

> *Here is an example ocm_sample.csv file:*
> ```
> sample_id,ocm_barcode_ids
> sample1,OB1
> sample2,OB2
> sample3,OB3|OB4
> ```

> *Where:*  

> - *sample_id:* Unique sample ID for this tagged sample. Must not contain whitespace, quote or comma characters. Each sample ID must be unique.
> - *ocm_barcode_ids:* Unique OCM ID(s) that the sample is tagged with.
>
> ***Example:*** `--ocm-sample ocm_sample.csv`

### 2.4 Flex (Fixed RNA) Arguments

`--probe-set PROBESET`
> **Probe set file.**   
> *type: file*
>   
> The probe set reference CSV file that matches the probes that were used in the capture. These are provided by 10x and can be found at [their downloads page](https://www.10xgenomics.com/support/software/cell-ranger/downloads). 
> 
> For those running the pipeline on the NIH's Biowulf, these references can also be found in the cell-seek references section in OpenOmics. The probe set reference files currently available in OpenOmics are:
> 
> - Human
>     - Chromium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv
>     - Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv
>     - Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv
> - Mouse
>     - Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv
>     - Chromium_Mouse_Transcriptome_Probe_Set_v1.1.1_GRCm39-2024-A.csv
> 
>
> ***Example:*** `--probe-set Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv`

---
`--probe-sample PROBESAMPLE`
> **Probe sample file.**   
> *type: file*
>   
> A probe sample CSV file containing sample to probe barcode information used for multi analysis. This file should be used if multiplexing was performed. This file should contain a unique sample ID for the sample, and the probe barcode(s) associated with that sample. If more than one probe barcode is associated with a sample then a | should be used to separate the tags. More information and examples about the samples section of the multi config file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

> *Here is an example probe_sample.csv file:*
> ```
> sample_id,probe_barcode_ids
> sample1,BC001
> sample2,BC002|BC003
> ```

> *Where:*  

> - *sample_id:* Unique sample ID for this tagged sample. Must not contain whitespace, quote or comma characters. Each sample ID must be unique.
> - *probe_barcode_ids:* Unique probe barcode(s) that the sample is captured with.
>
> ***Example:*** `--probe-sample probe_sample.csv`

---  

### 2.5 Analysis Options

Each of the following arguments are optional, and do not need to be provided.

`--features FEATURES`
> **Features file.**   
> *type: file*
>   
> A feature reference CSV file containing information for processing a feature barcode data. This file is used in feature barcode, and may be used in multi analysis. This file should contain a unique ID for the feature, a human readable name, sequence, feature type, read, and pattern. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-feature-bc-analysis).

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
>        * CRISPR Guide Capture
>        * Antibody Capture
>        * Custom
>
> ***Example:*** `--features features.csv`

---  
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
> This flag is currently not properly tested for situations where Cell Ranger demultiplexing is performed. Please double-check the output if applying this flag to those runs.
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
> This flag is currently not properly tested for situations where Cell Ranger demultiplexing is performed. Please double-check the output if applying this flag to those runs.
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
> This file can created in two different formats. The first one will contain the name of the sample and the number of cells to be forced to.
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
> The second format is only compatible when hashtag multiplexing is used and the number of cells needs to be forced for a specific hashtagged sample.
>
> *Here is an example forcecells.csv file:*
> ```
> Name,Sample,Cells
> Library1,Sample1,3000
> Library1,Sample2,5000
> ```
>
> *Where:*
>
> - *Library:* The name of the library that is provided as to Cell Ranger when running multi analysis. This should match the name that is given in the libraries.csv file.
> - *Sample:* The sample ID used for the associated hashtag. This will have to match the value used in the CMO/HTO/OCM sample file or the CMO reference file that is provided as input. If only a CMO reference file is provided, the pipeline default assigns each sample with the IDs used in the reference file.
> - *Cells:* The number of cells the sample should be forced to
>
> In this example, the hashtags associated with Sample1 and Sample2 in Library 1 will be run while being forced to 3000 and 5000 cells respectively. Any other libraries or samples that are processed will be run without using the force cells flag.
>
> ***Example:*** `--forcecells forcecells.csv`


### 2.6 Orchestration Options

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

### 2.7 Miscellaneous Options  
Each of the following arguments are optional, and do not need to be provided.

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
>
> Shows command's synopsis, help message, and an example command
>
> ***Example:*** `--help`

## 3. Example

### 3.1 GEX and VDJ

**Preparing Libraries File**

A CSV file needs to be created that specifies the modalities for each set of FASTQ files associated with each set of libraries to be processed by Cell Ranger. The file contains each Cell Ranger output sample name, FASTQ location information, FASTQ sample name, and library type. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-feature-bc-analysis).

The following an example libraries CSV file.

```
Name,Flowcell,Sample,Type
IL15_LNs,H7CNNBGXG,IL15_LNs,Gene Expression
IL15_LNs,H7CT7BGXG,IL15_LNs_VDJ,VDJ
WT,H7CNNBGXG,Control,Gene Expression
WT,H7CT7BGXG,Control,VDJ

```

The Flowcell column is used to identify the path the FASTQ files are located. This is in case there is a situation where the same FASTQ name is used for different modalities across two different runs. In the example above, the name Control was used for the FASTQ files in two different sequencing runs. The files located in the path H7CNNBGXG are associated with a gene expression capture, while the files located in the path H7CT7BGXG are associated with VDJ.

**Running the Pipeline**

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --mode slurm
```


### 3.2 Run with CMO

**Preparing Libraries File**

A CSV file needs to be created that specifies the modalities for each set of FASTQ files associated with each set of libraries to be processed by Cell Ranger. The file contains each Cell Ranger output sample name, FASTQ location information, FASTQ sample name, and library type. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-feature-bc-analysis).

The following an example libraries CSV file.

```
Name,Flowcell,Sample,Type
IL15_LNs,H7CNNBGXG,IL15_LNs,Gene Expression
IL15_LNs,H7CT7BGXG,IL15_LNs_CMO,Multiplexing Capture
WT,H7CNNBGXG,Control,Gene Expression
WT,H7CT7BGXG,Control,Multiplexing Capture

```

The Flowcell column is used to identify the path the FASTQ files are located. This is in case there is a situation where the same FASTQ name is used for different modalities across two different runs. In the example above, the name Control was used for the FASTQ files in two different sequencing runs. The files located in the path H7CNNBGXG are associated with a gene expression capture, while the files located in the path H7CT7BGXG are associated with CMO.

**Preparing CMO Sample File**

A CSV file needs to be created that specifies the sample to CMO information used for Cell Ranger demultiplexing. The CMO IDs should match the ones used by 10x. The same CMO sample will be used on all multi libraries. This file should contain a unique sample ID for the sample, and the CMO ID(s) associated with that sample. If more than one CMO ID is associated with a sample then a | should be used to separate the tags. More information and examples about the samples section of the multi config file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

The following is an example cmo_sample CSV file.

```
sample_id,cmo_ids
sample1,CMO301
sample2,CMO302|CMO303
```

**Running the Pipeline**

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --cmo-sample cmo_sample.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --cmo-sample cmo_sample.csv \
                  --mode slurm
```

### 3.3 Run with HTO

There can be a slight difference in the cells that are called when performing demultiplexing on HTO data when using the CMO flag or the HTO flag. The differences should be minor enough that it does not matter which method is used. If the data being processed is a continuation of an earlier project, then the method selected should match what was previously used.

#### 3.3.1 Using the CMO flag

**Preparing Libraries File**

A CSV file needs to be created that specifies the modalities for each set of FASTQ files associated with each set of libraries to be processed by Cell Ranger. The file contains each Cell Ranger output sample name, FASTQ location information, FASTQ sample name, and library type. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-feature-bc-analysis).

The following an example libraries CSV file.

```
Name,Flowcell,Sample,Type
IL15_LNs,H7CNNBGXG,IL15_LNs,Gene Expression
IL15_LNs,H7CT7BGXG,IL15_LNs_HTO,Multiplexing Capture
WT,H7CNNBGXG,Control,Gene Expression
WT,H7CT7BGXG,Control,Multiplexing Capture

```

The Flowcell column is used to identify the path the FASTQ files are located. This is in case there is a situation where the same FASTQ name is used for different modalities across two different runs. In the example above, the name Control was used for the FASTQ files in two different sequencing runs. The files located in the path H7CNNBGXG are associated with a gene expression capture, while the files located in the path H7CT7BGXG are associated with HTO.

**Preparing CMO Reference File**

A CSV file needs to be created that contains information about each of the hashtags that was used on the data. This file should contain a unique ID for the hashtag, a human readable name, sequence, feature type, read, and pattern. More information about the cmo reference file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

The following is an example cmo_reference.csv file:

```
id,name,sequence,feature_type,read,pattern
HTO1,HTO1,GTCAACTCTTTAGCG,Multiplexing Capture,R2,5P(BC)
HTO2,HTO2,TGATGGCCTATTGGG,Multiplexing Capture,R2,5P(BC)
HTO3,HTO3,CTTGCCGCATGTCAT,Multiplexing Capture,R2,5P(BC)
```

If no CMO sample file is provided, then the pipeline would use the hashtag ID as the sample ID when demultiplexing.



**Running the Pipeline**

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --cmo-reference cmo_reference.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --cmo-reference cmo_reference.csv \
                  --mode slurm
```
The following is an example of running `multi` while providing hashtag information.

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --cmo-reference cmo_reference.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --cmo-reference cmo_reference.csv \
                  --mode slurm
```

#### 3.3.2 Using the HTO flag

**Preparing Libraries File**

A CSV file needs to be created that specifies the modalities for each set of FASTQ files associated with each set of libraries to be processed by Cell Ranger. The file contains each Cell Ranger output sample name, FASTQ location information, FASTQ sample name, and library type. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-feature-bc-analysis).

The following an example libraries CSV file.

```
Name,Flowcell,Sample,Type
IL15_LNs,H7CNNBGXG,IL15_LNs,Gene Expression
IL15_LNs,H7CT7BGXG,IL15_LNs_HTO,Multiplexing Capture
WT,H7CNNBGXG,Control,Gene Expression
WT,H7CT7BGXG,Control,Multiplexing Capture

```

The Flowcell column is used to identify the path the FASTQ files are located. This is in case there is a situation where the same FASTQ name is used for different modalities across two different runs. In the example above, the name Control was used for the FASTQ files in two different sequencing runs. The files located in the path H7CNNBGXG are associated with a gene expression capture, while the files located in the path H7CT7BGXG are associated with HTO.

**Preparing Features File**

A CSV file needs to be created that contains information about each of the hashtags that was used on the data. This file should contain a unique ID for the hashtag, a human readable name, sequence, feature type, read, and pattern. More information about the cmo reference file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

If a features files already exists containing surface antibody tags that were used on the samples, then the hashtag information should be included in addition to this surface antibody information.

The following is an example features.csv file:

```
id,name,sequence,feature_type,read,pattern
CITE_CD64,CD64,AGTGGG,Antibody Capture,R2,5PNN(BC)
CITE_CD8,CD8,TCACCGT,Antibody Capture,R2,5PNNN(BC)
HTO1,HTO1,GTCAACTCTTTAGCG,Antibody Capture,R2,5P(BC)
HTO2,HTO2,TGATGGCCTATTGGG,Antibody Capture,R2,5P(BC)
HTO3,HTO3,CTTGCCGCATGTCAT,Antibody Capture,R2,5P(BC)
```


**Preparing HTO Sample File**

A CSV file needs to be created that specifies the sample to HTO information used for Cell Ranger demultiplexing. The HTO IDs should match the ones used in the features file. The same HTO sample will be used on all multi libraries. This file should contain a unique sample ID for the sample, and the HTO ID(s) associated with that sample. If more than one HTO ID is associated with a sample then a | should be used to separate the tags. More information and examples about the samples section of the multi config file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

The following is an example hto_sample CSV file.

```
sample_id,hashtag_ids
sample1,HTO1
sample2,HTO2|HTO3
```

**Running the Pipeline**

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --hto-sample hto_sample.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --hto-sample hto_sample.csv \
                  --mode slurm
```

#### 3.3.3 Using Seurat Demultiplexing

**Preparing Libraries File**

A CSV file needs to be created that specifies the modalities for each set of FASTQ files associated with each set of libraries to be processed by Cell Ranger. The file contains each Cell Ranger output sample name, FASTQ location information, FASTQ sample name, and library type. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-feature-bc-analysis).

The following an example libraries CSV file.

```
Name,Flowcell,Sample,Type
IL15_LNs,H7CNNBGXG,IL15_LNs,Gene Expression
IL15_LNs,H7CT7BGXG,IL15_LNs_HTO,Antibody Capture
WT,H7CNNBGXG,Control,Gene Expression
WT,H7CT7BGXG,Control,Antibody Capture

```

The Flowcell column is used to identify the path the FASTQ files are located. This is in case there is a situation where the same FASTQ name is used for different modalities across two different runs. In the example above, the name Control was used for the FASTQ files in two different sequencing runs. The files located in the path H7CNNBGXG are associated with a gene expression capture, while the files located in the path H7CT7BGXG are associated with HTO.

The HTO libraries are marked as Antibody Capture in order to treat them as cell surface antibodies, which would then be passed to Seurat to allow it to perform demultiplexing during downstream analysis.

**Preparing Features File**

A CSV file needs to be created that contains information about each of the hashtags that was used on the data. This file should contain a unique ID for the hashtag, a human readable name, sequence, feature type, read, and pattern. More information about the cmo reference file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

If a features files already exists containing surface antibody tags that were used on the samples, then the hashtag information should be included in addition to this surface antibody information.

The following is an example features.csv file:

```
id,name,sequence,feature_type,read,pattern
CITE_CD64,CD64,AGTGGG,Antibody Capture,R2,5PNN(BC)
CITE_CD8,CD8,TCACCGT,Antibody Capture,R2,5PNNN(BC)
HTO-1,HTO-1,GTCAACTCTTTAGCG,Antibody Capture,R2,5P(BC)
HTO-2,HTO-2,TGATGGCCTATTGGG,Antibody Capture,R2,5P(BC)
HTO-3,HTO-3,CTTGCCGCATGTCAT,Antibody Capture,R2,5P(BC)
```

The Seurat downstream script will search for features that start with `HTO-` or `HTO_`. If either is detected it will pull out those features as a separate HTO assay and perform demultiplexing during Seurat QC processing.


**Running the Pipeline**

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
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
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --mode slurm
```

### 3.4 Run with OCM

**Preparing Libraries File**

A CSV file needs to be created that specifies the modalities for each set of FASTQ files associated with each set of libraries to be processed by Cell Ranger. The file contains each Cell Ranger output sample name, FASTQ location information, FASTQ sample name, and library type. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-feature-bc-analysis).

The following an example libraries CSV file.

```
Name,Flowcell,Sample,Type
IL15_LNs,H7CNNBGXG,IL15_LNs,Gene Expression
IL15_LNs,H7CT7BGXG,IL15_LNs_OCM,Multiplexing Capture
WT,H7CNNBGXG,Control,Gene Expression
WT,H7CT7BGXG,Control,Multiplexing Capture

```

The Flowcell column is used to identify the path the FASTQ files are located. This is in case there is a situation where the same FASTQ name is used for different modalities across two different runs. In the example above, the name Control was used for the FASTQ files in two different sequencing runs. The files located in the path H7CNNBGXG are associated with a gene expression capture, while the files located in the path H7CT7BGXG are associated with OCM.


**Preparing OCM Sample File**

A CSV file needs to be created that specifies the sample to OCM information used for Cell Ranger demultiplexing. The OCM IDs should match the ones used by 10x. The same OCM sample will be used on all multi libraries. This file should contain a unique sample ID for the sample, and the OCM ID(s) associated with that sample. If more than one OCM ID is associated with a sample then a | should be used to separate the tags. More information and examples about the samples section of the multi config file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

The following is an example ocm_sample CSV file.

```
sample_id,ocm_barcode_ids
sample1,OB1
sample2,OB2|OB3
```

**Running the Pipeline**

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --ocm-sample ocm_sample.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --ocm-sample ocm_sample.csv \
                  --mode slurm
```

### 3.5 Run for Flex (Fixed RNA)

**Preparing Libraries File**

A CSV file needs to be created that specifies the modalities for each set of FASTQ files associated with each set of libraries to be processed by Cell Ranger. The file contains each Cell Ranger output sample name, FASTQ location information, FASTQ sample name, and library type. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/running-pipelines/cr-feature-bc-analysis).

The following an example libraries CSV file.

```
Name,Flowcell,Sample,Type
IL15_LNs,H7CNNBGXG,IL15_LNs,Gene Expression
WT,H7CNNBGXG,Control,Gene Expression
```

**Selecting Probe Set**

The probe set provided to the pipeline should match the one that was used during library preparation. These files are provided by 10x and can be found at [their downloads page](https://www.10xgenomics.com/support/software/cell-ranger/downloads). 

For those running the pipeline on the NIH's Biowulf, these references can also be found in the cell-seek references section in OpenOmics.

#### 3.5.1 Singleplex

**Running the Pipeline**

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --probe-set Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --probe-set Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv \
                  --mode slurm
```

#### 3.5.2 Multiplex

**Preparing Probe Sample File**

A CSV file needs to be created that specifies the sample to probe barcode information used for Cell Ranger demultiplexing.  The probe barcode IDs should match the ones used by 10x. The same probe barcode sample will be used on all multi libraries. This file should contain a unique sample ID for the sample, and the probe barcode(s) associated with that sample. If more than one probe barcode is associated with a sample then a | should be used to separate the tags. More information and examples about the samples section of the multi config file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

The following is an example probe_sample CSV file.

```
sample_id,probe_barcode_ids
sample1,BC001
sample2,BC002|BC003
```

**Running the Pipeline**

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --probe-set Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv \
                  --probe-sample probesample.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input H7CNNBGXG/*_R?_fastq.gz H7CT7BGXG/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --probe-set Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv \
                  --probe-sample probe_sample.csv \
                  --mode slurm
```


### 3.6 Run with Custom Filters

This has only been properly tested when run on samples that do not undergo Cell Ranger demultiplexing. 

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
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --filter filter.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --filter filter.csv \
                  --mode slurm
```

### 3.7 Adding Sample Level Metadata

This has only been properly tested when run on samples that do not undergo Cell Ranger demultiplexing. 

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
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --metadata metadata.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --metadata metadata.csv \
                  --mode slurm
```

### 3.8 Running while Forcing Cell Call

#### 3.8.1 General Use 

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
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
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
                  --libraries libraries.csv \
                  --forcecells forcecells.csv \
                  --mode slurm
```

#### 3.8.2 Cell Demultiplexing Specific Use

It is possible to force the number of cells that are called by Cell Ranger for the specific CMO/HTO/OCM sample. In this situation, Cell Ranger will override the demultiplexing algorithm and force the pipeline to use the provided number of cells. 

A CSV file needs to be created with the first column containing the name of the library (the Cell Ranger outputted name used in the libraries.csv file) and the second column containing the sample name used for the associated hashtag, while the third column would be the number of cells to force the sample to. Only the samples included in the CSV file will be run while forcing the cell call. Any other samples that are processed will use the default demultiplexing algorithm.

The following is an example of a force cells CSV file.

```
Name,Sample,Cells
Library1,Sample1,3000
Library1,Sample2,5000
```

Based on this file, Sample1 and Sample 2 in Library 1will be run while being forced to have 3000 and 5000 cells respectively.

The following is an example of running the pipeline for a CMO demultiplexing run.

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*_R?_fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --cmo-sample cmo_sample.csv \
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
                  --libraries libraries.csv \
                  --cmo-sample cmo_sample.csv \
                  --forcecells forcecells.csv \
                  --mode slurm
```