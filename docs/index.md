<div align="center">

  <h1 style="font-size: 250%">cell-seek 🔬</h1>

  <b><i>One single-cell pipeline to rule them all!</i></b><br>
  <a href="https://doi.org/10.5281/zenodo.101815241">
      <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.101815241.svg" alt="DOI">
  </a>
  <a href="https://github.com/OpenOmics/cell-seek/releases">
    <img alt="GitHub release" src="https://img.shields.io/github/v/release/OpenOmics/cell-seek?color=blue&include_prereleases">
  </a>
  <a href="https://hub.docker.com/repository/docker/skchronicles/chicyte">
    <img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/skchronicles/chicyte">
  </a><br>
  <a href="https://github.com/OpenOmics/cell-seek/actions/workflows/main.yaml">
    <img alt="tests" src="https://github.com/OpenOmics/cell-seek/workflows/tests/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/cell-seek/actions/workflows/docs.yml">
    <img alt="docs" src="https://github.com/OpenOmics/cell-seek/workflows/docs/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/cell-seek/issues">
    <img alt="GitHub issues" src="https://img.shields.io/github/issues/OpenOmics/cell-seek?color=brightgreen">
  </a>
  <a href="https://github.com/OpenOmics/cell-seek/blob/main/LICENSE">
    <img alt="GitHub license" src="https://img.shields.io/github/license/OpenOmics/cell-seek">
  </a>

  <p>
    This is the home of the pipeline, cell-seek. Its long-term goals: to accurately process and analyze single cell data like no pipeline before!
  </p>

</div>  


## Overview
Welcome to cell-seek's documentation! This guide is the main source of documentation for users that are getting started with [cell-seek](https://github.com/OpenOmics/cell-seek/).

The **`./cell-seek`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions:

 * [<code>cell-seek <b>run</b></code>](usage/run.md): Run the cell-seek pipeline with your input files.
 * [<code>cell-seek <b>genome</b></code>](usage/genome.md): Create a custom genome reference file that can be used with the pipeline.
 * [<code>cell-seek <b>unlock</b></code>](usage/unlock.md): Unlocks a previous runs output directory.
 * [<code>cell-seek <b>shinycell</b></code>](usage/shinycell.md): Build a ShinyCell2 application from a Seurat RDS object.
 * [<code>cell-seek <b>cache</b></code>](usage/cache.md): Cache remote resources locally, coming soon!

**cell-seek** is a comprehensive set of pipelines to perform the initial processing of different single cell technologies. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](usage/run.md) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/OpenOmics/cell-seek/issues).

## Contribute

This site is a living document, created for and by members like you. cell-seek is maintained by the members of NCBR and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/cell-seek).

## Citation

If you use this software, please cite it as below:  

=== "BibTex"

    ```
    @software{Chen_Kuhn_OpenOmics_cell-seek_2023,
      author       = {Vicky Chen and Skyler Kuhn},
      title        = {OpenOmics/cell-seek},
      month        = nov,
      year         = 2023,
      publisher    = {Zenodo},
      doi          = {10.5281/zenodo.10181524},
      url          = {https://doi.org/10.5281/zenodo.10181524}
    }
    ```

=== "APA"

    ```
    Vicky Chen, & Skyler Kuhn. (2023). OpenOmics/cell-seek [Computer software]. https://doi.org/10.5281/zenodo.10181524
    ```

For more citation style options, please visit the pipeline's [Zenodo page](https://doi.org/10.5281/zenodo.10181524).

## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
