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
pipeline_output += [join(workpath, "finalreport", "metric_summary.xlsx")]

# CellRanger counts intermediate cleanup
pipeline_output += expand(
    join(workpath, "cleanup", "{sample}.samplecleanup"),
    sample=samples
)


# Get set of input paths
input_paths = [os.path.dirname(p) for p in inputs]
input_paths_set = []
for p in input_paths:
    if not p in input_paths_set:
        input_paths_set.append(p)

# Functions and rules for processing Multiome data

# Function definitions
def count_introns(wildcards):
    """
    Wrapper to decide whether to include introns for counting.
    See config['options']['exclude_introns'] for the encoded value.
    """
    if exclude_introns:
        return('--gex-exclude-introns true')
    else:
        return('')


def count_bam(wildcards):
    """
    Wrapper to decide whether to create BAM files during Cell Ranger alignment - currently unused
    See config['options']['create_bam'] for the encoded value.
    """
    if create_bam:
        return('')
    else:
        return('--no-bam')


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
        lib = join(workpath, "{sample}_libraries.csv")
    output:
        html = join(workpath, "{sample}", "outs", "web_summary.html")
    log:
        err = "run_{sample}_10x_cellranger_count.err",
        log = "run_{sample}_10x_cellranger_count.log"
    params:
        rname = "count",
        prefix = "{sample}",
        reference = config["references"][genome]["arc_ref"],
        introns = count_introns
    envmodules: config["tools"]["cellranger-arc"]
    shell:
        """
        # Remove output directory
        # prior to running cellranger
        if [ -d '{params.prefix}' ]; then
            if ! [ -f '{output.html}' ]; then
                rm -rf '{params.prefix}/'
                cellranger-arc count \\
                    --id={params.prefix} \\
                    --reference={params.reference} \\
                    --libraries={input.lib} {params.introns} \\
                2>{log.err} 1>{log.log}
            fi
        else
            cellranger-arc count \\
                --id={params.prefix} \\
                --reference={params.reference} \\
                --libraries={input.lib} {params.introns} \\
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
        summarize = join("workflow", "scripts", "generateSummaryFiles_MULTIOME.py"),
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
        cr_temp = join(workpath, "{sample}", "SC_ATAC_GEX_COUNTER_CS")
    shell:
        """
        if [ -d '{params.cr_temp}' ]; then
            rm -r {params.cr_temp}
        fi
        """
