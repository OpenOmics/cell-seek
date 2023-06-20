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

# Function defitions
def count_intron(wildcards):
    """
    Wrapper to decide whether to include introns for counting.
    See config['options']['pre_mrna'] for the encoded value.
    """
    if exclude_introns:
        return('--include-introns=false')
    else:
        return('')


rule count:
    output:
        join(workpath, "{sample}", "outs", "web_summary.html")
    log:
        err = "run_{sample}_10x_cellranger_count.err",
        log ="run_{sample}_10x_cellranger_count.log"
    params:
        rname = "count",
        batch = "-l nodes=1:ppn=16,mem=96gb",
        prefix = "{sample}",
        transcriptome = config["references"][genome]["transcriptome"],
        excludeintrons = count_intron
    envmodules: config["tools"]["cellranger"]
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
            {params.excludeintrons} \\
        2>{log.err} 1>{log.log}
        """

print(expand(join(workpath, "finalreport", "summaries", "{sample}_web_summary.html"), sample=samples))
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
