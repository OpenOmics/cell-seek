# ShinyCell2

## About

**ShinyCell2** is an interactive web application for exploring single-cell analysis results. It lets you browse clusters, view marker genes, generate plots, and inspect metadata without writing any code.

Within `cell-seek`, a ShinyCell2 app is built from a Seurat RDS object using the <code>cell-seek <b>shinycell</b></code> sub command. With minimal configuration, it validates your inputs, creates the required configuration files, and writes the code needed to deploy an interactive single-cell web app.

You can find more information about ShinyCell2 at its [GitHub repository](https://github.com/the-ouyang-lab/ShinyCell2) and the OpenOmics fork of ShinyCell2 [here](https://github.com/OpenOmics/ShinyCell2).

## Who this is for

- **Building an app?** See the [<code>cell-seek <b>shinycell</b></code> command reference](shinycell_command.md) for the full list of arguments, examples, and the marker file format.
- **Exploring results someone shared with you?** See the [step-by-step guide for PIs and stakeholders](shinycell/PIs.md) for launching a Shiny App on the NIH Biowulf HPC via HPC OnDemand — no command-line experience required.

## What you need

In its most basic form, building an app requires only *3 inputs*:

- A Seurat RDS object file (with a set of valid cluster identity labels).
- A project title.
- An output directory to write the shiny application files.

You can optionally provide a TSV file with marker genes to add a DEG page to the application.

!!! important "Supported Assays"
    `ShinyCell2` does not support all `Seurat` object assays. **Mileage may vary with the Seurat Object you use.**
    **_HTO assays_** and **_Azimuth cell.annotation_** `(prediction.score.celltype.l*)` assays will be removed from your object when using the ShinyCell2 pipeline.

    These assay types are unsupported currently by the original ShinyCell2.

    For more information about support for assays and extended format support please consult the original ShinyCell2 repository: https://github.com/the-ouyang-lab/ShinyCell2
