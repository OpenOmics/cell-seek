<div align="center">
   
  <h1>baseline 🔬</h1>
  
  **_long pipeline name_**

  [![tests](https://github.com/OpenOmics/baseline/workflows/tests/badge.svg)](https://github.com/OpenOmics/baseline/actions/workflows/main.yaml) [![docs](https://github.com/OpenOmics/baseline/workflows/docs/badge.svg)](https://github.com/OpenOmics/baseline/actions/workflows/docs.yml) [![GitHub issues](https://img.shields.io/github/issues/OpenOmics/baseline?color=brightgreen)](https://github.com/OpenOmics/baseline/issues)  [![GitHub license](https://img.shields.io/github/license/OpenOmics/baseline)](https://github.com/OpenOmics/baseline/blob/main/LICENSE) 
  
  <i>
    This is the home of the pipeline, baseline. Its long-term goals: to accurately ...insert goal, to infer ..insert goal, and to boldly ..insert goal like no pipeline before!
  </i>
</div>

## Overview
Welcome to baseline! Before getting started, we highly recommend reading through [baseline's documentation](https://openomics.github.io/baseline/).

The **`./baseline`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>baseline <b>run</b></code>](https://openomics.github.io/baseline/usage/run/): Run the baseline pipeline with your input files.
 * [<code>baseline <b>unlock</b></code>](https://openomics.github.io/baseline/usage/unlock/): Unlocks a previous runs output directory.
 * [<code>baseline <b>cache</b></code>](https://openomics.github.io/baseline/usage/cache/): Cache remote resources locally, coming soon!

**baseline** is a comprehensive ...insert long description. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](https://openomics.github.io/baseline/usage/run/) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](https://openomics.github.io/baseline/faq/questions/) prior to [opening an issue on Github](https://github.com/OpenOmics/baseline/issues).

## Dependencies
**Requires:** `singularity>=3.5`  `snakemake>=6.0`

At the current moment, the pipeline uses a mixture of enviroment modules and docker images; however, this will be changing soon! In the very near future, the pipeline will only use docker images. With that being said, [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and [singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step of the pipeline will rely on versioned images from [DockerHub](https://hub.docker.com/orgs/nciccbr/repositories). Snakemake uses singularity to pull these images onto the local filesystem prior to job execution, and as so, snakemake and singularity will be the only two dependencies in the future.

## Installation
Please clone this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/baseline.git
# Change your working directory
cd baseline/
# Add dependencies to $PATH
# Biowulf users should run
module load snakemake singularity
# Get usage information
./baseline -h
```

## Contribute 
This site is a living document, created for and by members like you. baseline is maintained by the members of OpenOmics and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository](https://github.com/OpenOmics/baseline).

## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  