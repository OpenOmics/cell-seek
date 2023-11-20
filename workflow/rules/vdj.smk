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

def filterFastq(wildcards):
    return(','.join(set([os.path.dirname(i) for i in input_fastq if len(re.findall(f"{wildcards.sample}_[\w]*R2[\w]*.fastq.gz", i)) > 0])))


# Function definitions


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
        reference = config["references"][genome]["vdj_ref"],
        fastqs = filterFastq
    envmodules: config["tools"]["cellranger"]
    shell:
        """
        # Remove output directory
        # prior to running cellranger
        if [ -d '{params.prefix}' ]; then
            rm -rf '{params.prefix}/'
        fi

        cellranger vdj \\
            --id {params.prefix} \\
            --sample {params.prefix} \\
            --reference {params.reference} \\
            --fastqs {params.fastqs} \\
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

rule sampleCleanup:
    input:
        rules.count.output
    output:
        cleanup = touch(join(workpath, "cleanup", "{sample}.samplecleanup"))
    params:
        rname = "sampleCleanup",
        cr_temp = join(workpath, "{sample}", "SC_VDJ_ASSEMBLER_CS")
    shell:
        """
        rm -r {params.cr_temp}
        """
