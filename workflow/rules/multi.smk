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

    if cmo_ref != "None":
        flags.append(f"--cmoref {cmo_ref}")

    if cmo_sample != "None":
        flags.append(f"--cmosample {cmo_sample}")

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

rule multi:
    input:
        join(workpath, "{sample}.csv")
    output:
        join(workpath, "{sample}", "outs", "config.csv")
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
