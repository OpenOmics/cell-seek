wildcard_constraints:
    sample="|".join(lib_samples)

# Pipeline output definition

# Cell Ranger multi output
pipeline_output += expand(
    join(workpath, "{sample}", "outs", "config.csv"),
    sample=lib_samples
)

# Output from summaryFiles
pipeline_output += [join(workpath, "finalreport", "metric_summary.xlsx")]

# CellRanger counts intermediate cleanup
pipeline_output += expand(
    join(workpath, "cleanup", "{sample}.samplecleanup"),
    sample=lib_samples
)


#pipeline_output += [join(workpath, "Project_Cell_Filters.csv")]
pipeline_output += expand(join(workpath, "finalreport", "seurat", "Summary_QC_Report.html"), sample=lib_samples)
pipeline_output += expand(join(workpath, "seurat", "{sample}", "SeuratQCReportCopy_{sample}.complete"), sample=lib_samples)

# Get set of input paths
input_paths = [os.path.dirname(p) for p in inputs]
input_paths_set = []
for p in input_paths:
    if not p in input_paths_set:
        input_paths_set.append(p)

# Functions and rules for processing MULTI data

# Function definitions

def conditional_flags(wildcards):
    """
    Wrapper to decide what flags to provide to the multi config file generator.
    See config['options'] for the encoded values.
    """
    flags = []
    if exclude_introns:
        flags.append('--exclude_introns')

    if create_bam:
        flags.append('--create_bam')

    if features != "None":
        flags.append(f"--feature {features}")

    if CMO_REF != "None":
        flags.append(f"--cmoref {CMO_REF}")

    if CMO_SAMPLE != "None":
        flags.append(f"--cmosample {CMO_SAMPLE}")

    if HTO_SAMPLE != "None":
        flags.append(f"--htosample {HTO_SAMPLE}")

    if OCM_SAMPLE != "None":
        flags.append(f"--ocmsample {OCM_SAMPLE}")

    if PROBE_SAMPLE != "None":
        flags.append(f"--probesample {PROBE_SAMPLE}")

    if PROBE_SET != "None":
        flags.append(f"--probeset {PROBE_SET}")

    if CELLCOUNT_LIBRARIES:
        if wildcards.sample in CELLCOUNT_DICT.keys():
            flags.append(f"--multiplexforcecells {' '.join([','.join(i) for i in CELLCOUNT_DICT[wildcards.sample]])}")
    else:
        if wildcards.sample in CELLCOUNT_DICT.keys():
            flags.append(f"--forcecells {CELLCOUNT_DICT[wildcards.sample]}")

    if libraries != 'None':
      f = open(libraries, 'r')
      for line in f:
          if all([i in line for i in [wildcards.sample, 'VDJ']]):
              flags.append(f'--vdjref {config["references"][genome]["vdj_ref"]}')
              break

    return(' '.join(flags))

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



def get_last_two_path_components(path_str):
    """
    Extracts the last two components of a path using os.path.
    """
    path_str = os.path.split(path_str)[0]
    components = []
    head, tail = os.path.split(path_str)
    while tail:
        components.insert(0, tail)
        head, tail = os.path.split(head)

    if len(components) >= 2:
        return os.path.join(components[-2], components[-1])
    elif len(components) == 1:
        return components[0]
    else:
        return ""

def seuratQCSummarySamples(wildcards):
    """
    Wrapper to return the sample list into an R friendly input format
    """
    lib_samples_expanded = [get_last_two_path_components(i) for i in glob.glob(join(workpath, "seurat", "*", "*", "seur_cluster.rds"))]
    return("c('{}')".format("','".join(lib_samples_expanded)))

def seuratQCSummarySkippedSamples(wildcards):
    """
    Wrapper to return the sample list into an R friendly input format
    """
    lib_samples_expanded = set(get_last_two_path_components(i) for i in glob.glob(join(workpath, "seurat", "*", "*", "seur_cluster.rds")))
    lib_samples_expanded_2 = set(get_last_two_path_components(i) for i in glob.glob(join(workpath, "seurat", "*", "*", "seuratQC.complete")))
    
    return("c('{}')".format("','".join(lib_samples_expanded_2-lib_samples_expanded)))


def seuratQCReportCopy_output(wildcards):
    # note 1: ck_output is the same as OUTDIR

    # note 2: checkpoints will have attributes for each of the checkpoint
    # rules, accessible by name.
    #checkpoint_output = checkpoints.multi.get(**wildcards).output['outs']
    checkpoint_output = checkpoints.multi.get(**wildcards).output['outs']
    if int(CELLRANGER.split('.')[0]) >= 10:
      SMP = glob_wildcards(os.path.join(checkpoint_output, "{count_sample}", "sample_filtered_feature_bc_matrix")).count_sample
    else:
      SMP = glob_wildcards(os.path.join(checkpoint_output, "{count_sample}", "count", "sample_filtered_feature_bc_matrix")).count_sample
    return expand(join(workpath, "seurat", "{sample}", "{count_sample}", "copySeuratQCReport.complete"), **wildcards, count_sample=SMP)

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

rule multiConfig:
    input:
        join(workpath, "{sample}_libraries.csv")
    output:
        join(workpath, "{sample}.csv")
    params:
        #numcells = lambda wildcards:dict2[wildcards.sample],
        rname = "multiconfig",
        flags = conditional_flags,
        transcriptome = config["references"][genome]['gex_transcriptome'],
        cellranger = CELLRANGER,
        multiconfig = join("workflow", "scripts", "write_multiconfig.py")
    shell:
        """
        python {params.multiconfig} \\
            -o {output} \\
            --ref {params.transcriptome} \\
            -l {input} \\
            --cellranger {params.cellranger} \\
            {params.flags}
        """

checkpoint multi:
    input:
        join(workpath, "{sample}.csv")
    output:
        config = join(workpath, "{sample}", "outs", "config.csv"),
        outs = directory(join(workpath, "{sample}", "outs", "per_sample_outs"))
    log:
        err = "run_{sample}_10x_cellranger_count.err",
        log ="run_{sample}_10x_cellranger_count.log"
    params:
        rname = "count",
        batch = "-l nodes=1:ppn=16,mem=96gb",
        prefix = "{sample}",
        transcriptome = config["references"][genome]["gex_transcriptome"],
    envmodules: config["tools"]["cellranger"][CELLRANGER]
    shell:
        """
        # Remove output directory
        # prior to running cellranger
        if [ -d '{params.prefix}' ]; then
            if ! [ -f '{output}' ]; then
                rm -rf '{params.prefix}/'
                cellranger multi \\
                    --id={params.prefix} \\
                    --csv={input} \\
                2>{log.err} 1>{log.log}
            fi
        else
            cellranger multi \\
                --id={params.prefix} \\
                --csv={input} \\
            2>{log.err} 1>{log.log}
        fi
        """

rule summaryFiles:
    input:
        expand(join(workpath, "{sample}", "outs", "config.csv"), sample=lib_samples)
    output:
        join(workpath, "finalreport", "metric_summary.xlsx")
    params:
        rname = "sumfile",
        batch = "-l nodes=1:ppn=1",
        summarize = join("workflow", "scripts", "generateSummaryFiles_MULTI.py"),
    container: config["images"]["cite_base"]
    shell:
        """
        python {params.summarize}
        """

rule sampleCleanup:
    input:
       join(workpath, "{sample}", "outs", "config.csv") 
    output:
        cleanup = touch(join(workpath, "cleanup", "{sample}.samplecleanup"))
    params:
        rname = "sampleCleanup",
        cr_temp = join(workpath, "{sample}", "SC_MULTI_CS")
    shell:
        """
        if [ -d '{params.cr_temp}' ]; then
           rm -r {params.cr_temp}
        fi
        """

rule seuratQC:
    input:
        data = join(workpath, "{sample}", "outs", "per_sample_outs", "{count_sample}", "sample_filtered_feature_bc_matrix") if int(CELLRANGER.split('.')[0]) >= 10 else join(workpath, "{sample}", "outs", "per_sample_outs", "{count_sample}", "count", "sample_filtered_feature_bc_matrix")
    output:
        complete = touch(join(workpath, "seurat", "{sample}", "{count_sample}", "seuratQC.complete"))
#        rds = join(workpath, "seurat", "{sample}", "{count_sample}", "seur_cluster.rds"),
#        cell_filter = join(workpath, "seurat", "{sample}", "{count_sample}", "cell_filter_info.csv"),
#        matrix = [
#            join(workpath, "seurat", "{sample}", "{count_sample}", "cite-seq-matrix", "export_HTO_matrix.csv"),
#            join(workpath, "seurat", "{sample}", "{count_sample}", "cite-seq-matrix", "export_ADT_matrix.csv")
#        ]
    params:
        rname = "seuratQC",
        sample = "{count_sample}",
        outdir = join(workpath, "seurat", "{sample}", "{count_sample}"),
        seurat = join("workflow", "scripts", "seuratCiteSampleQC.R"),
        filter = filterFile,
        metadata = metadataFile,
        log = join(workpath, "seurat", "{sample}", "{count_sample}", "seurat.log")
    envmodules: config["tools"]["rversion"]
    shell:
        """
        unset __RLIBSUSER
        unset R_LIBS_USER

        Rscript --vanilla {params.seurat} \\
            --workdir {params.outdir} \\
            --datapath {input.data} \\
            --sample {params.sample} \\
            --project export \\
            {params.filter} \\
            {params.metadata} \\
            > {params.log}
        """


rule seuratQCReport:
    input:
        complete = join(workpath, "seurat", "{sample}", "{count_sample}", "seuratQC.complete")
        #rds = rules.seuratQC.output.complete,
        #cell_filter = rules.seuratQC.output.complete
    output:
        complete = touch(join(workpath, "seurat", "{sample}", "{count_sample}", "seuratQCReport.complete"))
    params:
        rname = "seuratQCReport",
        sample = "{sample} - {count_sample}",
        seuratdir = join(workpath, "seurat", "{sample}", "{count_sample}"),
        rds = join(workpath, "seurat", "{sample}", "{count_sample}", "seur_cluster.rds"),
        filter = filterFileBool,
        tmpdir = tmpdir,
        report = join(workpath, "seurat", "{sample}", "{count_sample}_QC_Report.html"),
        script = join(workpath, "workflow", "scripts", "seuratCiteSampleQCReport.Rmd")
    envmodules: config["tools"]["rversion"]
    shell:
        """
        unset __RLIBSUSER
        unset R_LIBS_USER

	if [ -f "{params.rds}" ]; then
          cd {params.tmpdir}
          cp {params.script} "./{params.sample}.Rmd"
          R -e "rmarkdown::render('{params.sample}.Rmd', params=list(seuratdir='{params.seuratdir}', sample='{params.sample}', defaultfilter={params.filter}), output_file='{params.report}')"
       fi
        """

rule copySeuratQCReport:
  input:
    complete = rules.seuratQCReport.output.complete
  output:
    complete = touch(join(workpath, "seurat", "{sample}", "{count_sample}", "copySeuratQCReport.complete"))
  params:
    report = rules.seuratQCReport.params.report,
    copied_report = join(workpath, "finalreport", "seurat", "{sample} - {count_sample}_QC_Report.html"),
    rname = "copySeuratQCReport"
  shell:
    """
    DIRECTORY=$(dirname "{params.copied_report}")
    if [ -f "{params.report}" ]; then
      if [[ ! -d $DIRECTORY ]]; then
        mkdir -p $DIRECTORY
      fi
      cp {params.report} '{params.copied_report}'
    fi
    """

rule finishAllSeuratQCReportCopy:
  input:
    report = seuratQCReportCopy_output
  output:
    complete = join(workpath, "seurat", "{sample}", "SeuratQCReportCopy_{sample}.complete")
  params:
    rname = "finishAllSeuratQCReportCopy"
  shell:
    """
    touch {output.complete}
    """

rule cellFilterSummary:
    input:
        expand(join(workpath, "seurat", "{sample}", "SeuratQCReportCopy_{sample}.complete"), sample=lib_samples)
    output:
        cell_filter_summary = join(workpath, "Project_Cell_Filters.csv")
    params:
        rname = "cellFilterSummary",
        seuratdir = join(workpath, "seurat", "*"),
        filename = "cell_filter_info.csv",
        script = join("workflow", "scripts", "cellFilterSummary.R")
    envmodules: config["tools"]["rversion"]
    shell:
        """
        unset __RLIBSUSER
        unset R_LIBS_USER

        Rscript {params.script} --datapath "{params.seuratdir}" --filename {params.filename} --output {output.cell_filter_summary}
        """

rule seuratQCSummaryReport:
    input:
        complete = expand(join(workpath, "seurat", "{sample}", "SeuratQCReportCopy_{sample}.complete"), sample=lib_samples),
        cell_filter = rules.cellFilterSummary.output.cell_filter_summary
    output:
        report = join(workpath, "seurat", "Summary_QC_Report.html")
    params:
        rname = "seuratQCSummaryReport",
        samples = seuratQCSummarySamples,
        test = seuratQCSummarySkippedSamples,
        seuratdir = join(workpath, "seurat"),
        cleanreport = join(workpath, "seurat", "Clean_Summary_QC_Report.html"),
        script = join(workpath, "workflow", "scripts", "seuratCiteSampleQCSummaryReport.Rmd")
    envmodules: config["tools"]["rversion"]
    shell:
        """
        unset __RLIBSUSER
        unset R_LIBS_USER

        R -e "rmarkdown::render('{params.script}', params=list(seuratdir='{params.seuratdir}', samples={params.samples}, cellfilter='{input.cell_filter}', skipped={params.test}), output_file='{output.report}')"
        if [[ "{params.test}" != "c('')" ]]; then
        R -e "rmarkdown::render('{params.script}', params=list(seuratdir='{params.seuratdir}', samples={params.samples}, cellfilter='{input.cell_filter}'), output_file='{params.cleanreport}')"
        fi
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
