# Pipeline output definition

# CellRanger counts, summary report
pipeline_output += expand(
            join(workpath, "{sample}", "outs", "web_summary.html"),
            sample=lib_samples
        )

# Generate Summaries File, summary report
pipeline_output += expand(
            join(workpath, "finalreport", "summaries", "{sample}_web_summary.html"),
            sample=lib_samples
        )

# CellRanger counts intermediate cleanup
pipeline_output += expand(
    join(workpath, "cleanup", "{sample}.samplecleanup"),
    sample=lib_samples
)

# Seurat inital sample QC
pipeline_output += expand(
    join(workpath, "seurat", "{sample}", "seur_cluster.rds"),
    sample=lib_samples
)

# Seurat sample QC reports
pipeline_output += expand(
    join(workpath, "finalreport", "seurat", "{sample}_QC_Report.html"),
    sample=lib_samples
)

if len(lib_samples) > 1:
  pipeline_output += [join(workpath, "finalreport", "seurat", "Summary_QC_Report.html")]


# Get set of input paths
input_paths = [os.path.dirname(p) for p in inputs]
input_paths_set = []
for p in input_paths:
    if not p in input_paths_set:
        input_paths_set.append(p)

# Functions and rules for processing CITE-seq data

# Function defintions

def force_cells(wildcards):
    """
    Wrapper to get the number of forced cells to use for processing if force cells was requested for the sample
    """
    if wildcards.sample in CELLCOUNT_DICT.keys():
        return(f"--force-cells {CELLCOUNT_DICT[wildcards.sample]}")
    else:
        return('')

def count_intron(wildcards):
    """
    Wrapper to decide whether to include introns for counting.
    See config['options']['exclude_introns'] for the encoded value.
    """
    if exclude_introns:
        return('--include-introns false')
    else:
        return('')

def count_bam(wildcards):
    """
    Wrapper to decide whether to create BAM files during Cell Ranger alignment.
    See config['options']['create_bam'] for the encoded value.
    """
    if create_bam:
        if CELLRANGER in ['7.1.0', '7.2.0']:
            return('')
        else:
            return('--create-bam true')
    else:
        if CELLRANGER in ['7.1.0', '7.2.0']:
            return('--no-bam')
        else:
            return('--create-bam false')

def filterFile(wildcards):
    """
    Wrapper to decide whether to provide a filter file for Seurat
    QC analysis.
    See config['options']['filter'] for the encoded value.
    """
    if filter_file == "None":
        return("")
    else:
        return(f"--filterfile {filter_file}")

def filterFileBool(wildcards):
    """
    Wrapper to get if a filter file was provided
    See config['options']['filter'] for the encoded value.
    """
    if filter_file == "None":
        return("TRUE")
    else:
        return("FALSE")

def metadataFile(wildcards):
    """
    Wrapper to decide whether to provide a metadata file for Seurat
    QC analysis.
    See config['options']['metadata'] for the encoded value.
    """
    if METADATA_FILE == "None":
        return("")
    else:
        return(f"--metadata {METADATA_FILE}")


def seuratQCSummarySamples(wildcards):
    """
    Wrapper to return the sample list into an R friendly input format
    """
    return("c('{}')".format("','".join(lib_samples)))


# Rule definitions
rule librariesCSV:
    output:
        expand(join(workpath, "{sample}_libraries.csv"), sample=lib_samples)
    params:
        rname = "libcsv",
        fastq = ",".join(input_paths_set),
        libraries = libraries,
        create_libs = join("workflow", "scripts", "create_library_files.py"),
    container: config["images"]["cite_base"]
    shell:
        """
        python {params.create_libs} \\
            {params.libraries} \\
            {params.fastq}
        """


rule count:
    input:
        lib = join(workpath, "{sample}_libraries.csv"),
        features = features
    output:
        html = join(workpath, "{sample}", "outs", "web_summary.html")
    log:
        err = "run_{sample}_10x_cellranger_count.err",
        log ="run_{sample}_10x_cellranger_count.log"
    params:
        rname = "count",
        batch = "-l nodes=1:ppn=16,mem=96gb",
        prefix = "{sample}",
        #        numcells = lambda wildcards:s2c[wildcards.sample],
        transcriptome = config["references"][genome]["cite_transcriptome"],
        introns = count_intron,
        createbam = count_bam,
        forcecells = force_cells
    envmodules: config["tools"]["cellranger"][CELLRANGER]
    shell:
        """
        # Remove output directory
        # prior to running cellranger
        if [ -d '{params.prefix}' ]; then
            if ! [ -f '{output.html}' ]; then
                rm -rf '{params.prefix}/'
                cellranger count \\
                    --id={params.prefix} \\
                    --transcriptome={params.transcriptome} \\
                    --libraries={input.lib} \\
                    --feature-ref={input.features} {params.introns} {params.createbam} {params.forcecells} \\
                2>{log.err} 1>{log.log}
            fi
        else
            cellranger count \\
                --id={params.prefix} \\
                --transcriptome={params.transcriptome} \\
                --libraries={input.lib} \\
                --feature-ref={input.features} {params.introns} {params.createbam} {params.forcecells} \\
            2>{log.err} 1>{log.log}
        fi
        """


rule summaryFiles:
    input:
        expand(join(workpath, "{sample}", "outs", "web_summary.html"), sample=lib_samples)
    output:
        join(workpath, "finalreport", "metric_summary.xlsx"),
        expand(join(workpath, "finalreport", "summaries", "{sample}_web_summary.html"), sample=lib_samples)
    params:
        rname = "sumfile",
        batch = "-l nodes=1:ppn=1",
        summarize = join("workflow", "scripts", "generateSummaryFiles.py"),
    container: config["images"]["cite_base"]
    shell:
        """
        python {params.summarize}
        """

rule seurat:
    input:
        join(workpath, "{sample}", "outs", "web_summary.html"),
    output:
        rds = join(workpath, "seurat", "{sample}", "seur_cite_cluster.rds")
    log:
        join("seurat", "{sample}", "seurat.log")
    params:
        rname = "seurat",
        sample = "{sample}",
        outdir = join(workpath, "seurat/{sample}"),
        data = join(workpath, "{sample}/outs/filtered_feature_bc_matrix/"),
        rawdata = join(workpath, "{sample}/outs/raw_feature_bc_matrix/"),
        seurat = join("workflow", "scripts", "seurat_adt.R"),
        rlibs = config["r_libs"]["ext"]
    # envmodules: "R/4.1"
    container: config["images"]["cite_base"]
    shell:
        """
        # Add external packages to R .libPath
        export R_LIBS_USER='{params.rlibs}'
        R --no-save --args \\
            {params.outdir} \\
            {params.data} \\
            {params.rawdata} \\
            {params.sample} \\
            {genome} \\
        < {params.seurat} > {log}
        """


rule seurat_rmd_report:
    input:
        join(workpath, "seurat", "{sample}", "seur_cite_cluster.rds"),
        azimuth = join(workpath, "azimuth", "{sample}", "azimuth_prediction.rds")
    output:
        html = join(workpath, "seurat", "{sample}", "{sample}_seurat.html")
    params:
        rname = "seurat_rmd_report",
        sample = "{sample}",
        outdir = join(workpath, "seurat/{sample}"),
        seurat = join("workflow", "scripts", "seurat_adt_plot.Rmd"),
        html = join(workpath, "seurat", "{sample}", "{sample}_seurat.html"),
        celltype = join(workpath, "azimuth", "{sample}")
    # envmodules: "R/4.1"
    container: config["images"]["cite_base"]
    shell:
        """
        R -e "rmarkdown::render('{params.seurat}', params=list(workdir = '{params.outdir}', sample='{params.sample}', celltype='{params.celltype}'), output_file = '{params.html}')"
        """


rule seurat_aggregate:
    input:
        rds = expand(join(workpath, "seurat", "{sample}", "seur_cite_cluster.rds"), sample=lib_samples)
    output:
        rds = join(workpath, "seurat", "SeuratAggregate", "multimode.integrated.rds")
    log:
        join("seurat", "SeuratAggregate", "seurat.log")
    params:
        rname = "seurat_aggregate",
        sample = "aggregate",
        outdir = join(workpath, "seurat", "SeuratAggregate"),
        seurat = join("workflow", "scripts", "seurat_adt_aggregate.R"),
        sample_ref = lib_samples[0]
    # envmodules: "R/4.1"
    container: config["images"]["cite_base"]
    shell:
        """
        R --no-save --args \\
            {params.outdir} \\
            {genome} \\
            {params.sample_ref} \\
            {input.rds} \\
        < {params.seurat} > {log}
        """


rule seurat_aggregate_rmd_report:
    input:
        join(workpath, "seurat", "SeuratAggregate", "multimode.integrated.rds")
    output:
        html = join(workpath, "seurat", "SeuratAggregate", "SeuratAggregate_seurat.html")
    params:
        rname = "seurat_aggregate_rmd_report",
        sample = "Aggregate",
        outdir = join(workpath, "seurat", "SeuratAggregate"),
        seurat = join("workflow", "scripts", "seurat_adt_aggregate_report.Rmd"),
        html = join(workpath, "seurat", "SeuratAggregate", "SeuratAggregate_seurat.html")
    # envmodules: "R/4.1"
    container: config["images"]["cite_base"]
    shell:
        """
        R -e "rmarkdown::render('{params.seurat}', params=list(workdir = '{params.outdir}', sample='{params.sample}'), output_file = '{params.html}')"
        """

rule sampleCleanup:
    input:
        rules.count.output
    output:
        cleanup = touch(join(workpath, "cleanup", "{sample}.samplecleanup"))
    params:
        rname = "sampleCleanup",
        cr_temp = join(workpath, "{sample}", "SC_RNA_COUNTER_CS")
    shell:
        """
        if [ -d '{params.cr_temp}' ]; then
            rm -r {params.cr_temp}
        fi
        """

rule seuratQC:
    input:
        join(workpath, "{sample}", "outs", "web_summary.html")
    output:
        rds = join(workpath, "seurat", "{sample}", "seur_cluster.rds"),
        cell_filter = join(workpath, "seurat", "{sample}", "cell_filter_info.csv"),
        matrix = [
            join(workpath, "seurat", "{sample}", "cite-seq-matrix", "export_HTO_matrix.csv"),
            join(workpath, "seurat", "{sample}", "cite-seq-matrix", "export_ADT_matrix.csv")
        ]
    log:
        join("seurat", "{sample}", "seurat.log")
    params:
        rname = "seuratQC",
        sample = "{sample}",
        outdir = join(workpath, "seurat", "{sample}"),
        data = join(workpath, "{sample}", "outs", "filtered_feature_bc_matrix"),
        seurat = join("workflow", "scripts", "seuratCiteSampleQC.R"),
        filter = filterFile,
        metadata = metadataFile
    envmodules: config["tools"]["rversion"]
    shell:
        """
        unset __RLIBSUSER
        unset R_LIBS_USER

        Rscript --vanilla {params.seurat} \\
            --workdir {params.outdir} \\
            --datapath {params.data} \\
            --sample {params.sample} \\
            --project export \\
            {params.filter} \\
            {params.metadata} \\
            > {log}
        """

rule seuratQCReport:
    input:
        rds = rules.seuratQC.output.rds,
        cell_filter = rules.seuratQC.output.cell_filter
    output:
        report = join(workpath, "seurat", "{sample}", "{sample}_QC_Report.html")
    params:
        rname = "seuratQCReport",
        sample = "{sample}",
        seuratdir = join(workpath, "seurat", "{sample}"),
        filter = filterFileBool,
        tmpdir = tmpdir,
        script = join(workpath, "workflow", "scripts", "seuratCiteSampleQCReport.Rmd")
    envmodules: config["tools"]["rversion"]
    shell:
        """
        unset __RLIBSUSER
        unset R_LIBS_USER

        cd {params.tmpdir}
        cp {params.script} ./{params.sample}.Rmd
        R -e "rmarkdown::render('{params.sample}.Rmd', params=list(seuratdir='{params.seuratdir}', sample='{params.sample}', defaultfilter={params.filter}), output_file='{output.report}')"
        """

rule copySeuratQCReport:
  input:
    report = rules.seuratQCReport.output.report
  output:
    report = join(workpath, "finalreport", "seurat", "{sample}_QC_Report.html")
  params:
    rname = "copySeuratQCReport"
  shell:
    """
    cp {input.report} {output.report}
    """

rule cellFilterSummary:
    input:
        cell_filters = expand(rules.seuratQC.output.cell_filter, sample=lib_samples)
    output:
        cell_filter_summary = join(workpath, "Project_Cell_Filters.csv")
    params:
        rname = "cellFilterSummary",
        seuratdir = join(workpath, "seurat"),
        filename = "cell_filter_info.csv",
        script = join("workflow", "scripts", "cellFilterSummary.R")
    envmodules: config["tools"]["rversion"]
    shell:
        """
        unset __RLIBSUSER
        unset R_LIBS_USER

        Rscript {params.script} --datapath {params.seuratdir} --filename {params.filename} --output {output.cell_filter_summary}
        """

rule seuratQCSummaryReport:
    input:
        rds = expand(rules.seuratQC.output.rds, sample=lib_samples),
        cell_filter = rules.cellFilterSummary.output.cell_filter_summary
    output:
        report = join(workpath, "seurat", "Summary_QC_Report.html")
    params:
        rname = "seuratQCSummaryReport",
        samples = seuratQCSummarySamples,
        seuratdir = join(workpath, "seurat"),
        script = join(workpath, "workflow", "scripts", "seuratCiteSampleQCSummaryReport.Rmd")
    envmodules: config["tools"]["rversion"]
    shell:
        """
        unset __RLIBSUSER
        unset R_LIBS_USER

        R -e "rmarkdown::render('{params.script}', params=list(seuratdir='{params.seuratdir}', samples={params.samples}, cellfilter='{input.cell_filter}'), output_file='{output.report}')"
        """

rule copySeuratQCSummaryReport:
  input:
    report = rules.seuratQCSummaryReport.output.report
  output:
    report = join(workpath, "finalreport", "seurat", "Summary_QC_Report.html")
  params:
    rname = "copySeuratQCSummaryReport"
  shell:
    """
    cp {input.report} {output.report}
    """
