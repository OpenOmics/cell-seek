# Pipeline output definition

# Single sample libraries files for cellranger count
pipeline_output += expand(
            join(workpath, "{sample}_libraries.csv"),
            sample=lib_samples
        )

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

# Get set of input paths
input_paths = [os.path.dirname(p) for p in inputs]
input_paths_set = []
for p in input_paths:
    if not p in input_paths_set:
        input_paths_set.append(p)

# Functions and rules for processing CITE-seq data

# Function defintions
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
        join(workpath, "{sample}", "outs", "web_summary.html")
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
        createbam = count_bam
    envmodules: config["tools"]["cellranger"][CELLRANGER]
    shell:
        """
        # Remove output directory
        # prior to running cellranger
        if [ -d '{params.prefix}' ]; then
            rm -rf '{params.prefix}/'
        fi

        cellranger count \\
            --id={params.prefix} \\
            --transcriptome={params.transcriptome} \\
            --libraries={input.lib} \\
            --feature-ref={input.features} \\
            {params.introns} \\
            {params.createbam} \\
        2>{log.err} 1>{log.log}
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
        rm -r {params.cr_temp}
        """
