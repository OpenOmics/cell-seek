# CellSight

## About

**CellSight** is an interactive web application for exploring single-cell analysis results. It lets you browse clusters, view marker genes, generate plots, and inspect metadata without writing any code.

Within `cell-seek`, a CellSight app is built from a Seurat RDS object using the <code>cell-seek <b>cellsight</b></code> sub command. With minimal configuration, it validates your inputs, creates the required configuration files, and writes the code needed to deploy an interactive single-cell web app.

You can find more information about CellSight at its [GitHub repository](https://github.com/OpenOmics/CellSight).

## Who this is for

- **Building an app?** See the [<code>cell-seek <b>cellsight</b></code> command reference](cellsight_command.md) for the full list of arguments, examples, and the marker file format.
- **Exploring results someone shared with you?** See the [step-by-step guide for PIs and stakeholders](cellsight/PIs.md) for launching a Shiny App on the NIH Biowulf HPC via HPC OnDemand — no command-line experience required.

## What you need

In its most basic form, building an app requires only *3 inputs*:

- A Seurat RDS object file (with a set of valid cluster identity labels).
- A project title.
- An output directory to write the shiny application files.

You can optionally provide a TSV file with marker genes to add a DEG page to the application.

!!! important "Supported Assays"
    `CellSight` does not support all `Seurat` object assays. **Mileage may vary with the Seurat Object you use.**
    **_HTO assays_** and **_Azimuth cell.annotation_** `(prediction.score.celltype.l*)` assays will be removed from your object when using the CellSight pipeline.

    For more information about supported assays and extended format support please consult the [CellSight repository](https://github.com/OpenOmics/CellSight).
