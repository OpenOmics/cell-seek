# Pipeline output definition
from textwrap import dedent

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

# scATAC preliminary analysis outputs
## sample specific QC reports
pipeline_output += expand(
    join(workpath, "scATAC_analysis", "{sample}", "{sample}.QC_Report.html"),
    sample=samples
)
## cohort level QC report
pipeline_output += [join(workpath, "scATAC_analysis", "cohort", "cell_filter_info.csv")]

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

def force_cells(wildcards):
    """
    Wrapper to get the number of forced cells to use for processing if force cells was requested for the sample
    """
    if wildcards.sample in CELLCOUNT_DICT.keys():
        return(f"--force-cells {CELLCOUNT_DICT[wildcards.sample]}")
    else:
        return('')

rule count:
    output:
        html = join(workpath, "{sample}", "outs", "web_summary.html"),
        matrix = join(workpath, "{sample}", "outs", "filtered_peak_bc_matrix.h5"),
        fragments = join(workpath, "{sample}", "outs", "fragments.tsv.gz"),
        filterfile = join(workpath, "{sample}", "outs", "singlecell.csv")
    log:
        err = "run_{sample}_10x_cellranger_count.err",
        log ="run_{sample}_10x_cellranger_count.log"
    params:
        rname = "count",
        batch = "-l nodes=1:ppn=16,mem=96gb",
        id = "{sample}",
        sample = sample_rename,
        reference = config["references"][genome]["atac_ref"],
        fastqs = filterFastq,
        forcecells = force_cells
    envmodules: config["tools"]["cellranger-atac"]
    threads: 24
    shell:
        dedent("""
        # Remove output directory
        # prior to running cellranger
        if [ -d '{params.id}' ]; then
            if ! [ -f '{output.html}' ]; then
              rm -rf '{params.id}/'
              cellranger-atac count \\
                  --id {params.id} \\
                  --sample {params.sample} \\
                  --localcores {threads} \\
                  --reference {params.reference} \\
                  --fastqs {params.fastqs} {params.forcecells} \\
              2>{log.err} 1>{log.log}
            fi
        else
            cellranger-atac count \\
                --id {params.id} \\
                --sample {params.sample} \\
                --localcores {threads} \\
                --reference {params.reference} \\
                --fastqs {params.fastqs} {params.forcecells} \\
            2>{log.err} 1>{log.log}
        fi
        """)

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
        html = join(workpath, "{sample}", "outs", "web_summary.html")
    output:
        cleanup = touch(join(workpath, "cleanup", "{sample}.samplecleanup"))
    params:
        rname = "sampleCleanup",
        cr_temp = join(workpath, "{sample}", "SC_ATAC_COUNTER_CS")
    shell:
        """
        if [ -d '{params.cr_temp}' ]; then
            rm -r {params.cr_temp}
        fi
        """


rule prelim_analysis_one:
    input:
        matrix = join(workpath, "{sample}", "outs", "filtered_peak_bc_matrix.h5"),
        fragments = join(workpath, "{sample}", "outs", "fragments.tsv.gz"),
        filterfile = join(workpath, "{sample}", "outs", "singlecell.csv")
    output:
        filter_info = join(workpath, "scATAC_analysis", "{sample}", "cell_filter_info.csv"),
        report = join(workpath, "scATAC_analysis", "{sample}", "{sample}.QC_Report.html")
    params:
        rname = "prelim_analysis_one",
        script = join("/opt", "scripts", "signacSampleQC.R"),
        scriptrmd = join("/opt", "scripts", "signacSampleQC.Rmd"),
        genes = join(config["references"][genome]["atac_ref"], "genes", "genes.gtf.gz"),
        genome = genome,
        outdir = lambda wc: join(workpath, "scATAC_analysis", wc.sample),
        project = lambda wc: f"Preliminary QC Report for Cell-seek Sample {wc.sample} Analysis"
    container: config["images"]["signac_base"]
    shell:
        dedent("""
        Rscript {params.script} \\
            --sample {wildcards.sample} \\
            --featurematrix {input.matrix} \\
            --fragments {input.fragments} \\
            --barcodes {input.filterfile} \\
            --genes {params.genes} \\
            --genome {params.genome} \\
            --project "{params.project}" \\
            -o {params.outdir}
        R -e "rmarkdown::render('{params.scriptrmd}',
                params=list(signacdir='/data/OpenOmics/dev/datasets/input_artifacts/test_cellseek_atac/new_scripts/signacSampleQC_output',
                            thresholds='{output.filter_info}',
                            sample='{wildcards.sample}',
                            defaultfilter=TRUE),
                            output_file='{params.outdir}/{wildcards.sample}/{wildcards.sample}.QC_Report.html')"
        """)



rule prelim_analysis_all:
    input:
        matrix = expand(join(workpath, "{sample}", "outs", "filtered_peak_bc_matrix.h5"), sample=samples),
        fragments = expand(join(workpath, "{sample}", "outs", "fragments.tsv.gz"), sample=samples),
        filterfile = expand(join(workpath, "{sample}", "outs", "singlecell.csv"), sample=samples)
    output:
        filter_info = join(workpath, "scATAC_analysis", "cohort", "cell_filter_info.csv"),
        report = join(workpath, "scATAC_analysis", "cohort", "Cohort_QC_Report.html")
    params:
        rname = "prelim_analysis_all",
        script = join("/opt", "scripts", "signacMultiSampleQC.R"),
        genes = join(config["references"][genome]["atac_ref"], "genes", "genes.gtf.gz"),
        genome = genome,
        sids = ','.join(samples),
        scriptrmd = join("/opt", "scripts", "signacMultiSampleQC.Rmd"),
        outdir = join(workpath, "scATAC_analysis", "cohort"),
        project = "Preliminary QC Report for Cell-seek Multi-Sample Analysis"
    container: config["images"]["signac_base"]
    shell:
        dedent("""
        Rscript {params.script} \\
            --sample {params.sids} \\
            --matrix {input.matrix} \\
            --fragments {input.fragments} \\
            --barcodes {input.filterfile} \\
            --genes {params.genes} \\
            --genome {params.genome} \\
            --project "{params.project}" \\
            -o {params.outdir} 
        R -e "rmarkdown::render('{params.scriptrmd}',
                params=list(signacdir='{params.outdir}',
                            thresholds='{output.filter_info}',
                            sample='{params.sids}',
                            defaultfilter=TRUE),
                            output_file='{params.outdir}/Cohort_QC_Report.html')"
        """)