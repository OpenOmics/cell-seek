# Pipeline output definition

# CellRanger counts, summary report
pipeline_output += expand(
    join(workpath, "{sample}", "outs", "web_summary.html"),
    sample=samples
)

# Generate Summaries File, summary report
pipeline_output += expand(
    join(workpath, "finalreport", "summaries", "{sample}_web_summary.html"),
    sample=samples
)

# CellRanger counts intermediate cleanup
pipeline_output += expand(
    join(workpath, "cleanup", "{sample}.samplecleanup"),
    sample=samples
)

#CellRanger aggregate
if aggr != "":
    # CellRanger aggregate analysis
    pipeline_output += [join(workpath, 'aggregate.complete')]

    # CellRanger aggregate intermediate cleanup
    pipeline_output += [join(workpath, "cleanup", "aggregate.aggregatecleanup")]

# Seurat inital sample QC
pipeline_output += expand(
    join(workpath, "seurat", "{sample}", "seur_cluster.rds"),
    sample=samples
)

# Seurat sample QC reports
pipeline_output += expand(
    join(workpath, "seurat", "{sample}", "{sample}_QC_Report.html"),
    sample=samples
)

# Cell Filter Summary File
pipeline_output += [join(workpath, "Project_Cell_Filters.csv")]

# Function definitions
def filterFastq(wildcards):
    return(','.join(set([os.path.dirname(i) for i in input_fastq if len(re.findall(f"{wildcards.sample}_[\w]*R2[\w]*.fastq.gz", i)) > 0])))


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
        return('')
    else:
        return('--no-bam')

def aggr_norm(wildcards):
    """
    Wrapper to decide which normalization method to use for Cell Ranger aggregate.
    If no option was provided while running then aggregate will not be run.
    See config['options']['aggregate'] for the encoded value.
    """
    supported = ['mapped', 'none']
    if aggr.lower() in supported:
        return(aggr.lower())
    else:
        raise Exception(f"\nAn unsupported input was provided for the aggregate flag. Please check the aggregate value under options in the config file to change it so that either the entry is left blank or is one of the following: {supported}\n\nThe currently provided value is: {aggr}\n")

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


rule count:
    output:
        html = join(workpath, "{sample}", "outs", "web_summary.html")
    log:
        err = "run_{sample}_10x_cellranger_count.err",
        log ="run_{sample}_10x_cellranger_count.log"
    params:
        rname = "count",
        batch = "-l nodes=1:ppn=16,mem=96gb",
        prefix = "{sample}",
        transcriptome = config["references"][genome]["gex_transcriptome"],
        excludeintrons = count_intron,
        createbam = count_bam,
        fastqs = filterFastq
    envmodules: config["tools"]["cellranger"]
    shell:
        """
        # Remove output directory
        # prior to running cellranger
        if [ -d '{params.prefix}' ]; then
            rm -rf '{params.prefix}/'
        fi

        cellranger count \\
            --id {params.prefix} \\
            --sample {params.prefix} \\
            --transcriptome {params.transcriptome} \\
            --fastqs {params.fastqs} \\
            {params.excludeintrons} \\
            {params.createbam} \\
        2>{log.err} 1>{log.log}
        """

rule summaryFiles:
    input:
        expand(join(workpath, "{sample}", "outs", "web_summary.html"), sample=samples)
    output:
        join(workpath, "finalreport", "metric_summary.xlsx"),
        expand(join(workpath, "finalreport", "summaries", "{sample}_web_summary.html"), sample=samples)
    params:
        rname = "sumfile",
        batch = "-l nodes=1:ppn=1",
        summarize = join("workflow", "scripts", "generateSummaryFiles.py"),
    container: config["images"]["cite_base"]
    shell:
        """
        python {params.summarize}
        """

rule aggregateCSV:
    input:
        expand(join(workpath, "{sample}", "outs", "web_summary.html"), sample=samples)
    output:
        join(workpath, "AggregatedDatasets.csv")
    params:
        rname = "aggregateCSV",
        script = join("workflow", "scripts", "generateAggregateCSV_GEX.py"),
        analysis = workpath
    shell:
       """
       python {params.script} {params.analysis}
       """

rule aggregate:
    input:
        csv=join(workpath, "AggregatedDatasets.csv")
    output:
        touch(join(workpath, "aggregate.complete"))
    log:
        err="run_10x_aggregate.err",
        log="run_10x_aggregate.log"
    params:
        rname = "aggregate",
        id = "AggregateDatasets",
        norm = aggr_norm
    envmodules: config["tools"]["cellranger"]
    shell:
        """
        cellranger aggr \\
            --id {params.id} \\
            --csv {input.csv} \\
            --normalize={params.norm} \\
        2>{log.err} 1>{log.log}
        """

rule seuratQC:
    input:
        join(workpath, "{sample}", "outs", "web_summary.html")
    output:
        rds = join(workpath, "seurat", "{sample}", "seur_cluster.rds"),
        cell_filter = join(workpath, "seurat", "{sample}", "cell_filter_info.csv")
    log:
        join("seurat", "{sample}", "seurat.log")
    params:
        rname = "seuratQC",
        sample = "{sample}",
        outdir = join(workpath, "seurat", "{sample}"),
        data = join(workpath, "{sample}", "outs", "filtered_feature_bc_matrix"),
        seurat = join("workflow", "scripts", "seuratSampleQC.R"),
        filter = filterFile
    shell:
        """
        module load R/4.3.0
        Rscript {params.seurat} \\
            --workdir {params.outdir} \\
            --datapath {params.data} \\
            --sample {params.sample} \\
            {params.filter} \\
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
	tmpdir = tmpdir,
        script = join(workpath, "workflow", "scripts", "seuratSampleQCReport.Rmd")
    shell:
        """
        module load R/4.3.0
	cd {params.tmpdir}
	cp {params.script} ./{params.sample}.Rmd
        R -e "rmarkdown::render('{params.sample}.Rmd', params=list(seuratdir='{params.seuratdir}', sample='{params.sample}'), output_file='{output.report}')"
        """

rule cellFilterSummary:
    input:
        cell_filters = expand(rules.seuratQC.output.cell_filter, sample=samples)
    output:
        cell_filter_summary = join(workpath, "Project_Cell_Filters.csv")
    params:
        rname = "cellFilterSummary",
        seuratdir = join(workpath, "seurat"),
        filename = "cell_filter_info.csv",
        script = join("workflow", "scripts", "cellFilterSummary.R")
    shell:
        """
        module load R/4.3.0
        Rscript {params.script} --datapath {params.seuratdir} --filename {params.filename} --output {output.cell_filter_summary}
        """

rule sampleCleanup:
    input:
        html = rules.count.output.html
    output:
        cleanup = touch(join(workpath, "cleanup", "{sample}.samplecleanup"))
    params:
        rname = "sampleCleanup",
        cr_temp = join(workpath, "{sample}", "SC_RNA_COUNTER_CS")
    shell:
        """
        rm -r {params.cr_temp}
        """

rule aggregateCleanup:
    input:
        rules.aggregate.output
    output:
        cleanup = touch(join(workpath, "cleanup", "{sample}.aggregatecleanup"))
    params:
        rname = "aggregateCleanup",
        cr_temp = join(workpath, "AggregateDatasets", "SC_RNA_AGGREGATOR_CS")
    shell:
        """
        rm -r {params.cr_temp}
        """
