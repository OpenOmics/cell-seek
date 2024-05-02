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

# Function definitions
def filterFastq(wildcards):
    """
    Wrapper to get a comma separated list of the directories where the FASTQ files associated with the sample are located
    """
    filter_paths = []
    for sample in sample_rename(wildcards).split(','):
        filter_paths += [os.path.dirname(i) for i in input_fastq if len(re.findall(f"{sample}_[\w]*R2[\w]*.fastq.gz", i)) > 0]
    return(','.join(set(filter_paths)))
    #return(','.join(set([os.path.dirname(i) for i in input_fastq if len(re.findall(f"{wildcards.sample}_[\w]*R2[\w]*.fastq.gz", i)) > 0])))

def sample_rename(wildcards):
    """
    Wrapper to get the FASTQ file names to use processing if the sample was requested to be renamed
    """
    if wildcards.sample in RENAME_DICT.values():
        names = [i[0] for i in RENAME_DICT.items() if wildcards.sample == i[1]]
        return(','.join(names))
    else:
        return(wildcards.sample)


rule count:
    output:
        html = join(workpath, "{sample}", "outs", "web_summary.html")
    log:
        err = "run_{sample}_10x_cellranger_count.err",
        log ="run_{sample}_10x_cellranger_count.log"
    params:
        rname = "count",
        batch = "-l nodes=1:ppn=16,mem=96gb",
        id = "{sample}",
        sample = sample_rename,
        reference = config["references"][genome]["atac_ref"],
        fastqs = filterFastq
    envmodules: config["tools"]["cellranger-atac"]
    shell:
        """
        # Remove output directory
        # prior to running cellranger
        if [ -d '{params.id}' ]; then
            if ! [ -f '{output.html}' ]; then
              rm -rf '{params.id}/'
              cellranger-atac count \\
                  --id {params.id} \\
                  --sample {params.sample} \\
                  --reference {params.reference} \\
                  --fastqs {params.fastqs} \\
              2>{log.err} 1>{log.log}
            fi
        else
            cellranger-atac count \\
                --id {params.id} \\
                --sample {params.sample} \\
                --reference {params.reference} \\
                --fastqs {params.fastqs} \\
            2>{log.err} 1>{log.log}
        fi
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
        summarize = join("workflow", "scripts", "generateSummaryFiles_ATAC.py"),
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
        cr_temp = join(workpath, "{sample}", "SC_ATAC_COUNTER_CS")
    shell:
        """
        rm -r {params.cr_temp}
        """
