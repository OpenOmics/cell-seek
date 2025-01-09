# <code>cell-seek <b>run</b></code>

## 1. About
The `cell-seek` executable is composed of several inter-related sub commands. Please see `cell-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>cell-seek <b>run</b></code> sub command in more detail. The following page lists the options applicable to each of the different pipelines available to select via the <code><b>--pipeline</b></code> flag. With minimal configuration, the **`run`** sub command enables you to start running cell-seek pipeline.

Setting up the cell-seek pipeline is fast and easy! In its most basic form, <code>cell-seek <b>run</b></code> only has *four required inputs*.

## 2. Synopsis
```text
$ cell-seek run [--help] \
      [--dry-run] [--job-name JOB_NAME] [--mode {{slurm,local}}] \
      [--sif-cache SIF_CACHE] [--singularity-cache SINGULARITY_CACHE] \
      [--silent] [--threads THREADS] [--tmp-dir TMP_DIR] \
      [--cellranger {8.0.0, ...} \
      [--aggregate {{mapped, none}}] [--exclude-introns] \
      [--libraries LIBRARIES] [--features FEATURES] \
      [--filter FILTER] [--metadata METADATA] [--create-bam] \
      [--rename RENAME] \
      --input INPUT [INPUT ...] \
      --output OUTPUT \
      --pipeline {gex, ...} \
      --genome {hg38, ...}
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of FastQ (globbing is supported) to analyze via `--input` argument, an output directory to store results via `--output` argument, the version of the pipeline to run via `--pipeline` argument, and the reference genome to use via `--genome` argument.

Use you can always use the `-h` option for information on a specific command.

The following is a breakdown of the required and optional arguments for each of the versions of the pipeline.

### 2.1 GEX

#### 2.1.1 Required Arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input FastQ file(s) or Cell Ranger folder(s).**  
> *type: file(s) or folder(s)*  
>
> FastQ Input: One or more FastQ files can be provided. The pipeline does NOT support single-end data. From the command-line, each input file should separated by a space. Multiple input FastQ files per sample can be provided. Globbing is supported! This makes selecting FastQ files easy. Input FastQ files should always be gzipp-ed.
>
> ***Example:*** `--input .tests/*.R?.fastq.gz`
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
> A custom reference genome can also be provided via a json file. Additional information for creating this json file can be found in [<code>cell-seek <b>genome</b></code>](../genome).
>
> For prebuilt references please select one of the following options: hg38, mm10, hg2024, mm2024
>
> ***Example:*** `--genome hg2024`

---  
  `--cellranger {7.1.0, 7.2.0, 8.0.0, 9.0.0}`
> **The version of Cell Ranger to run.**   
> *type: string*
>   
> This option specifies which version of Cell Ranger to use when running GEX, VDJ, CITE, or MULTI pipelines. Please select one of the following options: 7.1.0, 7.2.0, 8.0.0, 9.0.0
>
> ***Example:*** `--cellranger 7.1.0`

#### 2.1.2 Analysis Options

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
> Sample metadata file. A CSV file containing sample level metadata information that will be included as new metadata columns during QC analysis. The file should contain a header row with Sample as the column name for the sample IDs, and the name of each metadata column that will be added to their associated samples. Each row is then the entries for each sample with the values that will be included as metadata. If no file is provided then no metadata will be added to the samples. If a cell is left blank for a sample then the metadata column for that sample would be an empty string. This flag is currently only applicable when dealing with GEX projects.
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


### 2.2 VDJ

#### 2.2.1 Required Arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input FastQ file(s) or Cell Ranger folder(s).**  
> *type: file(s) or folder(s)*  
>
> FastQ Input: One or more FastQ files can be provided. The pipeline does NOT support single-end data. From the command-line, each input file should separated by a space. Multiple input FastQ files per sample can be provided. Globbing is supported! This makes selecting FastQ files easy. Input FastQ files should always be gzipp-ed.
>
> ***Example:*** `--input .tests/*.R?.fastq.gz`
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
  `--pipeline vdj`
> **The pipeline to run.**   
> *type: string*
>   
> This option selects the version of the pipeline to run. The documentation provided is based on choosing the option for VDJ.
>
> ***Example:*** `--pipeline vdj`

---  
  `--genome {hg38, mm10, custom.json}`
> **Reference genome.**   
> *type: string*
>   
> This option defines the reference genome of the samples. cell-seek does comes bundled with prebuilt reference files for human and mouse samples, e.g. hg38 or mm10. Since there is no 2024 release VDJ reference, if hg2024 or mm2024 is selected the VDJ reference CR 7.1 release will be used for human, and CR 7.0 release will be used for mouse.
>
> A custom reference genome can also be provided via a json file. Additional information for creating this json file can be found in [<code>cell-seek <b>genome</b></code>](../genome).
>
> For prebuilt references please select one of the following options: hg38, mm10
>
> ***Example:*** `--genome hg38`

---  
  `--cellranger {7.1.0, 7.2.0, 8.0.0, 9.0.0}`
> **The version of Cell Ranger to run.**   
> *type: string*
>   
> This option specifies which version of Cell Ranger to use when running GEX, VDJ, CITE, or MULTI pipelines. Please select one of the following options: 7.1.0, 7.2.0, 8.0.0, 9.0.0
>
> ***Example:*** `--cellranger 7.1.0`

#### 2.2.2 Analysis Options

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

### 2.3 CITE

#### 2.3.1 Required Arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input FastQ file(s) or Cell Ranger folder(s).**  
> *type: file(s) or folder(s)*  
>
> FastQ Input: One or more FastQ files can be provided. The pipeline does NOT support single-end data. From the command-line, each input file should separated by a space. Multiple input FastQ files per sample can be provided. Globbing is supported! This makes selecting FastQ files easy. Input FastQ files should always be gzipp-ed.
>
> ***Example:*** `--input .tests/*.R?.fastq.gz`
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
> A custom reference genome can also be provided via a json file. Additional information for creating this json file can be found in [<code>cell-seek <b>genome</b></code>](../genome).
>
> For prebuilt references please select one of the following options: hg38, mm10, hg2024, mm2024
>
> ***Example:*** `--genome hg38`

---  
  `--cellranger {7.1.0, 7.2.0, 8.0.0, 9.0.0}`
> **The version of Cell Ranger to run.**   
> *type: string*
>   
> This option specifies which version of Cell Ranger to use when running GEX, VDJ, CITE, or MULTI pipelines. Please select one of the following options: 7.1.0, 7.2.0, 8.0.0, 9.0.0
>
> ***Example:*** `--cellranger 7.1.0`


#### 2.3.2 Conditionally Required Arguments

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
> - *Flowcell:* The flowcell ID that contains the FASTQ files for this set of data.  
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

#### 2.3.3 Analysis Options

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

### 2.4 MULTI

There are multiple different combinations of library types that may result in the use of Cell Ranger `multi` analysis. Any combination that combines GEX and VDJ data for cell calls, or the use of HTO with the Cell Ranger hashtag caller would need `multi` analysis.

#### 2.4.1 Required Arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

`--input INPUT [INPUT ...]`  
> **Input FastQ file(s) or Cell Ranger folder(s).**  
> *type: file(s) or folder(s)*  
>
> FastQ Input: One or more FastQ files can be provided. The pipeline does NOT support single-end data. From the command-line, each input file should separated by a space. Multiple input FastQ files per sample can be provided. Globbing is supported! This makes selecting FastQ files easy. Input FastQ files should always be gzipp-ed.
>
> ***Example:*** `--input .tests/*.R?.fastq.gz`
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
> A custom reference genome can also be provided via a json file. Additional information for creating this json file can be found in [<code>cell-seek <b>genome</b></code>](../genome).
>
> For prebuilt references please select one of the following options: hg38, mm10, hg2024, mm2024
>
> ***Example:*** `--genome hg38`

---  
  `--cellranger {7.1.0, 7.2.0, 8.0.0, 9.0.0}`
> **The version of Cell Ranger to run.**   
> *type: string*
>   
> This option specifies which version of Cell Ranger to use when running GEX, VDJ, CITE, or MULTI pipelines. Please select one of the following options: 7.1.0, 7.2.0, 8.0.0, 9.0.0
>
> ***Example:*** `--cellranger 7.1.0`

#### 2.4.2 Conditionally Required Arguments

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
> - *Flowcell:* The flowcell ID that contains the FASTQ files for this set of data.  
> - *Sample:* Name that was used when demultiplexing, this should match the FASTQ files.  
> - *Type:* library type for each sample. List of supported options:  
>        * Gene Expression
>        * CRISPR Guide Capture
>        * Antibody Capture
>        * Multiplexing Capture
>        * VDJ
>        * Custom
>
> ***Example:*** `--libraries libraries.csv`

#### 2.4.3 Analysis Options

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
`--cmo-reference CMOREFERENCE`
> **CMO reference file.**   
> *type: file*
>   
> A CMO reference CSV file containing information for processing hashtag data, which is used in multi analysis if custom hashtags will be processed. This file should contain a unique ID for the hashtag, a human readable name, sequence, feature type, read, and pattern. More information about the cmo reference file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

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
> - *feature_type:* Type of the feature. This should always be multiplexing capture.
> - *read:* Specifies which RNA sequencing read contains the Feature Barcode sequence. Must be R1 or R2, but in most cases R2 is the correct read.
> - *pattern:* Specifies how to extract the sequence of the feature barcode from the read.
>
> ***Example:*** `--cmo-reference cmo_reference.csv`

---  
`--cmo-sample CMOSAMPLE`
> **CMO sample file.**   
> *type: file*
>   
> A CMO sample CSV file containing sample to hashtag information used for multi analysis. CMO IDs should match the ones used in the CMO reference file. If no CMO reference is provided then CMO ID should match the ones used by 10x. The same CMO sample will be used on all multi libraries. This file should contain a unique sample ID for the sample, and the CMO IDs associated with that sample. If more than one CMO ID is associated with a sample then a | should be used to separate the tags. More information and examples about the samples section of the multi config file and its requirements can be found on the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/analysis/cr-3p-multi).

> *Here is an example cmo_sample.csv file:*
> ```
> sample_id,cmo_ids
> sample1,CMO301
> sample2,CMO302|CMO303
> ```

> *Where:*  

> - *sample_id:* Unique sample ID for this hashtagged sample. Must not contain whitespace, quote or comma characters. Each sample ID must be unique.
> - *cmo_ids:* Unique CMO ID(s) that the sample is hashtagged with. Must match either entries in the cmo_reference.csv file or 10x CMO IDs.
>
> ***Example:*** `--cmo-sample cmo_sample.csv`

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
> - *Sample:* The sample ID used for the associated hashtag. This will have to match the value used in the CMO sample file or the CMO reference file that is provided as input. If only a CMO reference file is provided, the pipeline default assigns each hashtag with the IDs of HTO_1, HTO_2, etc.
> - *Cells:* The number of cells the sample should be forced to
>
> In this example, the hashtags HTO_1 and HTO_2 in Library 1 will be run while being forced to 3000 and 5000 cells respectively. Any other libraries or samples that are processed will be run without using the force cells flag.
>
> ***Example:*** `--forcecells forcecells.csv`

### 2.5 ATAC

#### 2.5.1 Required Arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input FastQ file(s) or Cell Ranger folder(s).**  
> *type: file(s) or folder(s)*  
>
> FastQ Input: One or more FastQ files can be provided. The pipeline does NOT support single-end data. From the command-line, each input file should separated by a space. Multiple input FastQ files per sample can be provided. Globbing is supported! This makes selecting FastQ files easy. Input FastQ files should always be gzipp-ed.
>
> ***Example:*** `--input .tests/*.R?.fastq.gz`
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
  `--genome {hg38, mm10, custom.json}`
> **Reference genome.**   
> *type: string*
>   
> This option defines the reference genome of the samples. cell-seek does comes bundled with prebuilt reference files for human and mouse samples, e.g. hg38 or mm10.
>
> A custom reference genome can also be provided via a json file. Additional information for creating this json file can be found in [<code>cell-seek <b>genome</b></code>](../genome).
>
> For prebuilt references please select one of the following options: hg38, mm10
>
> ***Example:*** `--genome hg38`


#### 2.5.2 Analysis Options
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

### 2.6 Multiome

#### 2.6.1 Required Arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

`--input INPUT [INPUT ...]`  
> **Input FastQ file(s) or Cell Ranger folder(s).**  
> *type: file(s) or folder(s)*  
>
> FastQ Input: One or more FastQ files can be provided. The pipeline does NOT support single-end data. From the command-line, each input file should separated by a space. Multiple input FastQ files per sample can be provided. Globbing is supported! This makes selecting FastQ files easy. Input FastQ files should always be gzipp-ed.
>
> ***Example:*** `--input .tests/*.R?.fastq.gz`
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
  `--pipeline multiome`
> **The pipeline to run.**   
> *type: string*
>   
> This option selects the version of the pipeline to run. The documentation provided is based on choosing the option for multiome.
>
> ***Example:*** `--pipeline multiome`

---  
  `--genome {hg38, mm10, custom.json}`
> **Reference genome.**   
> *type: string*
>   
> This option defines the reference genome of the samples. cell-seek does comes bundled with prebuilt reference files for human and mouse samples, e.g. hg38 or mm10.
>
> A custom reference genome can also be provided via a json file. Additional information for creating this json file can be found in [<code>cell-seek <b>genome</b></code>](../genome).
>
> For prebuilt references please select one of the following options: hg38, mm10
>
> ***Example:*** `--genome hg38`


#### 2.6.2 Conditionally Required Arguments

The following arguments are only required when FastQ files are used as input. They are not required when Cell Ranger output file is used as input.

`--libraries LIBRARIES`
> **Libraries file.**   
> *type: file*
>   
> A CSV file containing information about each library. It contains each sample's name, flowcell, demultiplexed name, and library type. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/count#library-csv).

> *Here is an example libraries.csv file:*
> ```
> Name,Flowcell,Sample,Type
> IL15_LNs,H7CNNBGXG,IL15_LNs,Gene Expression
> IL15_LNs,H7CT7BGXG,IL15_LNs_BC,Antibody Capture
> ```

> *Where:*

> - *Name:* name of the sample passed to Cell Ranger.  
> - *Flowcell:* The flowcell ID that contains the FASTQ files for this set of data.  
> - *Sample:* Name that was used when demultiplexing, this should match the FASTQ files.  
> - *Type:* library type for each sample. List of supported options:  
>        * Gene Expression
>        * Chromatin Accessibility
>
> ***Example:*** `--libraries libraries.csv`


#### 2.6.3 Analysis Options

The multiome pipeline currently does not have any applicable analysis flags.


### 2.7 Orchestration Options

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

### 2.8 Miscellaneous Options  
Each of the following arguments are optional, and do not need to be provided.

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
>
> Shows command's synopsis, help message, and an example command
>
> ***Example:*** `--help`

## 3. Example

### 3.1 GEX Example
```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*.R?.fastq.gz \
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
./cell-seek run --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --pipeline gex \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --mode slurm
```

### 3.2 VDJ Example
```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --pipeline vdj \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --pipeline vdj \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --mode slurm
```

### 3.3 CITE Example

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*.R?.fastq.gz \
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
./cell-seek run --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --pipeline cite \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --features features.csv \
                  --mode slurm
```

### 3.4 MULTI Example

#### 3.4.1 GEX and VDJ

The following is an example of running `multi` when combining GEX and VDJ runs.

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*.R?.fastq.gz \
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
./cell-seek run --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --mode slurm
```

#### 3.4.2 Including HTO

The following is an example of running `multi` while providing hashtag information.

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*.R?.fastq.gz \
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
./cell-seek run --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multi \
                  --genome hg38 \
                  --cellranger 8.0.0 \
                  --libraries libraries.csv \
                  --cmo-reference cmo_reference.csv \
                  --mode slurm
```

### 3.5 ATAC Example

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --pipeline atac \
                  --genome hg38 \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --pipeline atac \
                  --genome hg38 \
                  --mode slurm
```

### 3.6 Multiome Example

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cell-seek run --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multiome \
                  --genome hg38 \
                  --libraries libraries.csv \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the cell-seek pipeline
# The slurm mode will submit jobs to
# the cluster. It is recommended running
# the pipeline in this mode.
./cell-seek run --input .tests/*.R?.fastq.gz \
                  --output /data/$USER/output \
                  --pipeline multiome \
                  --genome hg38 \
                  --libraries libraries.csv \
                  --mode slurm
```
